#include "simulation/Simulation.h"

#include "io/FileReader.h"
#include "io/VTKWriter.h"
#include "physics/LinkedCellParticleContainer.h"

#include <chrono>  // for benchmarking
#include <memory>
#include <utility>
// process signal handling, TODO: Implement signal handling for all simulation types (?)
#include <atomic>
#include <csignal>

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif  // SPDLOG_ACTIVE_LEVEL
#include "spdlog/spdlog.h"

// Atomic flag for loop control
std::atomic<bool> simulationRunning{true};

void sigint_handler(int signal) {
  if (signal == SIGINT) {
    simulationRunning = false;
  }
}

Simulation::Simulation(SimulationConfig& config)
    : endTime(config.tEnd),
      dt(config.deltaT),
      simulationMode(config.simulationMode),
      writeFrequency(config.writeFrequency),
      outputBasename(std::move(config.outputBasename)),
      particleGenerator(std::move(config.particleGenerator)) {
  // Initialize particle container and force calculation strategy
  switch (config.containerType) {
    case ContainerType::DIRECT:
      particles = std::make_unique<ParticleContainer>();
      if (config.useParallelization) {
        forceCalc =
            std::make_unique<LennardJonesForceParallel>(*particles, *config.epsilon, *config.sigma, *config.cutoff);
      } else {
        forceCalc = std::make_unique<LennardJonesForce>(*particles, *config.epsilon, *config.sigma, *config.cutoff,
                                                        config.gravity.value_or(0.0));
      }
      break;
    case ContainerType::LINKED:
      particles =
          std::make_unique<LinkedCellParticleContainer>(*config.domainSize, *config.cutoff, *config.boundaryTypes);
      forceCalc = std::make_unique<LennardJonesForce>(
          *particles, *config.epsilon, *config.sigma, *config.cutoff,
          config.gravity.value_or(0.0));  // TODO: Does parallelization not work with linked cell containers?
      break;
  }

  // Initialize thermostat
  if (config.thermostatConfig) {
    const ThermostatConfig& thermoConfig = *config.thermostatConfig;

    // Target temperature and initial temperature initialization
    double actualTargetTemp;
    if (thermoConfig.tempTarget) {
      actualTargetTemp = *thermoConfig.tempTarget;
      needToAutoSetTargetTemperature = false;
    } else if (thermoConfig.tempInit) {
      // Fallback: Target temp missing but T_init exists, target = init
      actualTargetTemp = *thermoConfig.tempInit;
      needToAutoSetTargetTemperature = false;
    } else {
      // Fallback: Target temp and initial temp both missing, target = current (calculated).
      actualTargetTemp = 0.0;  // we don't know "Current" yet, so set dummy and flag it.
      needToAutoSetTargetTemperature = true;
    }

    // Create thermostat object
    thermostat =
        std::make_unique<Thermostat>(*particles, thermoConfig.nThermostat, actualTargetTemp,
                                     thermoConfig.tempDelta.value_or(std::numeric_limits<double>::infinity()));

    initialTemperature = thermoConfig.tempInit;  // copy optional, used in simulation setup
  } else {
    thermostat = nullptr;
  }
}

Simulation::~Simulation() = default;

void Simulation::plotParticles(const int iteration) const {
  outputWriter::VTKWriter::plotParticles(*particles, outputBasename, iteration);
}

void Simulation::run() {
  setupSimulation();  // set up particles and objects

  switch (simulationMode) {
    case SimulationMode::BENCHMARK:
      runBenchmark();
      break;
    case SimulationMode::FILE_OUTPUT:
      runFileOutput();
      break;
  }
}

void Simulation::runFileOutput() {
  const int thermoFrequency = thermostat ? thermostat->getUpdateFrequency() : 1;

  double currentTime = 0;
  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (currentTime < endTime) {
    // calculate new x
    forceCalc->calculateX(dt);
    for (auto& p : *particles) {
      p.setOldF(p.getF());  // store f(t_n) for v update
    }
    // calculate new f
    forceCalc->calculateF();
    // calculate new v
    forceCalc->calculateV(dt);

    iteration++;
    if (thermostat && (iteration % thermoFrequency == 0)) {
      thermostat->updateTemperature();
    }
    if (iteration % writeFrequency == 0) {
      plotParticles(iteration);
    }

    currentTime += dt;
  }
}

void Simulation::runBenchmark() {
  // Set up interrupting signal (CTRL + C)
  std::signal(SIGINT, sigint_handler);
  simulationRunning = true;

  using namespace std::chrono;  // used for benchmarking

  // Benchmark begin
  const auto chronoStart = steady_clock::now();

  const int thermoFrequency = thermostat ? thermostat->getUpdateFrequency() : 1;

  double currentTime = 0;
  long iteration = 0;

  // For this loop, we assume current x, current F and current v are known
  while (currentTime < endTime) {
    // Calculate the new positions of the particles
    forceCalc->calculateX(dt);
    for (auto& p : *particles) {
      p.setOldF(p.getF());  // store F(t_n) for v update
    }
    // Calculate new forces and velocities
    forceCalc->calculateF();
    forceCalc->calculateV(dt);

    iteration++;
    if (thermostat && (iteration % thermoFrequency == 0)) {
      thermostat->updateTemperature();
    }

    currentTime += dt;
  }

  const auto chronoEnd = steady_clock::now();
  const auto elapsed = duration_cast<duration<double>>(chronoEnd - chronoStart).count();

  std::signal(SIGINT, SIG_DFL);
  spdlog::set_level(spdlog::level::info);
  if (!simulationRunning) {
    SPDLOG_INFO("Simulation stopped early by user (SIGINT).");
  } else {
    SPDLOG_INFO("Benchmark finished normally.");
  }
  SPDLOG_INFO("Time elapsed: {} s", elapsed);
  SPDLOG_INFO("Total iterations: {}", iteration);
  if (iteration > 0) {
    double timePerIteration = elapsed / static_cast<double>(iteration);
    SPDLOG_INFO("Mean time per iteration: {:.6f} s", timePerIteration);
  }
  spdlog::set_level(spdlog::level::off);
}

void Simulation::setupSimulation() {
  // Generate particles
  particleGenerator->generate(*particles);

  // Set up special thermostat situations (if applicable)
  if (thermostat) {
    // Apply T_init
    if (initialTemperature.has_value()) {
      SPDLOG_INFO("Applying initial thermostat temperature (Brownian motion): {}", *initialTemperature);
      thermostat->initializeTemperature(*initialTemperature);
    }

    // Fix target temperature (if it was missing)
    if (needToAutoSetTargetTemperature) {
      double currentTemp = thermostat->calculateCurrentTemperature();
      thermostat->setTargetTemperature(currentTemp);
      SPDLOG_INFO("Thermostat target set to initial calculated system temperature: {}", currentTemp);
    }
  }
}
