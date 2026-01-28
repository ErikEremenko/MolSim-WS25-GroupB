#include "simulation/Simulation.h"

#include "io/FileReader.h"
#ifdef ENABLE_VTK_OUTPUT
#include "io/VTKWriter.h"
#endif
#include "io/CheckpointWriter.h"
#include "physics/LinkedCellParticleContainer.h"

#include <atomic>
#include <chrono>  // for benchmarking
#include <cmath>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <memory>
#include <utility>

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
      currentTime(config.startTime),
      startTime(config.startTime),
      startIteration(config.startIteration),
      // --- TODO: Refactor these variables, look at Simulation.h for more info ---
      epsilon(1.0),
      sigma(1.0),
      cutoff(config.linkedCellCutoff.value_or(3.0)),
      gravity(0.0),
      // --------------------------------------------------------------------------
      dimensions(config.dimensions.value_or(3)),
      domainSize(config.domainSize.value_or(std::array<double, 3>{0.0, 0.0, 0.0})),
      membraneDimY(config.membraneDimY.value_or(0)),
      simulationMode(config.simulationMode),
      writeFrequency(config.writeFrequency),
      checkpointFrequency((config.simulationMode == SimulationMode::BENCHMARK) ? 0 : config.checkpointFrequency),
      particleGenerator(std::move(config.particleGenerator)),
      outputBasename(std::move(config.outputBasename)) {

  // Initialize boundary strings
  if (config.boundaryTypes) {
    for (size_t i = 0; i < 6; ++i) {
      switch ((*config.boundaryTypes)[i]) {
        case BoundaryType::OUTFLOW:
          boundaryTypeStrings[i] = "OUTFLOW";
          break;
        case BoundaryType::REFLECTIVE:
          boundaryTypeStrings[i] = "REFLECTIVE";
          break;
        case BoundaryType::PERIODIC:
          boundaryTypeStrings[i] = "PERIODIC";
          break;
      }
    }
  } else {
    boundaryTypeStrings.fill("OUTFLOW");
  }

  // Initialize particle container
  switch (config.containerType) {
    case ContainerType::DIRECT:
      particles = std::make_unique<ParticleContainer>();
      SPDLOG_INFO("Using DIRECT particle container (O(n^2) force calculation)");
      break;
    case ContainerType::LINKED:
      if (!config.domainSize || !config.linkedCellCutoff || !config.boundaryTypes) {
        SPDLOG_ERROR("LINKED container requires domainSize, linkedCellCutoff, and boundaryTypes to be set!");
        throw std::runtime_error("Missing required parameters for LINKED container");
      }
      particles = std::make_unique<LinkedCellParticleContainer>(*config.domainSize, *config.linkedCellCutoff,
                                                                *config.boundaryTypes);
      SPDLOG_INFO("Using LINKED cell particle container (cutoff={})", *config.linkedCellCutoff);
      break;
  }

  // Initialize forces
  for (auto& forceConfig : config.forceConfigs) {
    switch (forceConfig.forceType) {
      case ForceType::LENNARD_JONES:
        forces.push_back(std::make_unique<LennardJonesForce>(*particles, *forceConfig.epsilon, *forceConfig.sigma,
                                                             *forceConfig.cutoff));
        SPDLOG_INFO("Initialized LennardJonesForce (epsilon={}, sigma={}, cutoff={})", *forceConfig.epsilon,
                    *forceConfig.sigma, *forceConfig.cutoff);
        break;

      case ForceType::SMOOTHED_LJ:
        forces.push_back(std::make_unique<SmoothedLJForce>(*particles, *forceConfig.epsilon, *forceConfig.sigma,
                                                           *forceConfig.cutoff, *forceConfig.rl));
        SPDLOG_INFO("Initialized SmoothedLJForce (epsilon={}, sigma={}, r_c={}, r_l={})", *forceConfig.epsilon,
                    *forceConfig.sigma, *forceConfig.cutoff, *forceConfig.rl);
        break;

      case ForceType::TRUNCATED_LJ:
        forces.push_back(std::make_unique<TruncatedLJForce>(*particles));
        SPDLOG_INFO("Initialized TruncatedLJForce (repulsive-only, uses per-particle sigma/epsilon)");
        break;

      case ForceType::GLOBAL_GRAVITY: {
        int axis = forceConfig.gravityAxis.value_or(1);  // Default to y-axis
        forces.push_back(std::make_unique<GlobalGravityForce>(*particles, *forceConfig.gravity, axis));
        const char* axisName = (axis == 0) ? "x" : (axis == 1) ? "y" : "z";
        SPDLOG_INFO("Initialized GlobalGravityForce (g={}, axis={})", *forceConfig.gravity, axisName);
        break;
      }

      case ForceType::HARMONIC_MEMBRANE:
        forces.push_back(
            std::make_unique<HarmonicMembraneForce>(*particles, *forceConfig.stiffness, *forceConfig.avgBondLength));
        SPDLOG_INFO("Initialized HarmonicMembraneForce (k={}, r0={})", *forceConfig.stiffness,
                    *forceConfig.avgBondLength);
        break;

      case ForceType::CONSTANT_FORCE:
        forces.push_back(std::make_unique<ConstantForce>(*particles, forceConfig.forceX.value_or(0.0),
                                                         forceConfig.forceY.value_or(0.0),
                                                         forceConfig.forceZ.value_or(0.0), *forceConfig.endTime,
                                                         currentTime, forceConfig.targetIndices, membraneDimY));
        SPDLOG_INFO("Initialized ConstantForce (F=({}, {}, {}), end_time={}, targets={})",
                    forceConfig.forceX.value_or(0.0), forceConfig.forceY.value_or(0.0),
                    forceConfig.forceZ.value_or(0.0), *forceConfig.endTime, forceConfig.targetIndices.size());
        break;
    }
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
    thermostat = std::make_unique<Thermostat>(*particles, thermoConfig.nThermostat, actualTargetTemp,
                                              thermoConfig.tempDelta.value_or(std::numeric_limits<double>::infinity()));

    initialTemperature = thermoConfig.tempInit;  // copy optional, used in simulation setup
  } else {
    thermostat = nullptr;
  }
}

Simulation::~Simulation() = default;

void Simulation::plotParticles(const int iteration) const {
#ifdef ENABLE_VTK_OUTPUT
  outputWriter::VTKWriter::plotParticles(*particles, outputBasename, iteration);
  SPDLOG_DEBUG("Successfully wrote particles to file, iteration={}", iteration);
#else
  SPDLOG_WARN("VTK output disabled, skipping plotParticles for iteration {}", iteration);
#endif
}

void Simulation::writeCheckpoint(const int iteration, const double time) const {
  // Determine output directory based on current working directory
  std::string outputDirectory = "build/output/checkpoints";
  if (std::filesystem::exists("CMakeCache.txt")) {
    outputDirectory = "output/checkpoints";
  }

  outputWriter::CheckpointWriter::writeCheckpoint(
      *particles, outputBasename + "_checkpoint_" + std::to_string(iteration) + ".yaml", iteration, time,
      outputBasename, writeFrequency, checkpointFrequency, endTime, dt, epsilon, sigma, cutoff, gravity, domainSize,
      boundaryTypeStrings, outputDirectory);
}

void Simulation::run() {
  setupSimulation();

  switch (simulationMode) {
    case SimulationMode::BENCHMARK:
      SPDLOG_INFO("Simulation loop starting in benchmark mode...");
      runBenchmark();
      break;
    case SimulationMode::FILE_OUTPUT:
      SPDLOG_INFO("Simulation loop starting in file output mode...");
      runFileOutput();
      break;
  }
}

void Simulation::runFileOutput() {
  const int thermoFrequency = thermostat ? thermostat->getUpdateFrequency() : 1;

  // Use member currentTime instead of local variable
  int iteration = startIteration;

  // For this loop, we assume: current positions, forces and velocities are known
  while (currentTime < endTime) {
    // Calculate the new position of the particles
    ForceCalc::calculateX(*particles, dt);

    // Store the force from the previous time stop for velocity update
    /* TODO: Optimization idea: introduce a 'first force' flag to ForceCalc to avoid
    * iterating an extra time over the particles, thus eliminating the loop below.
    */
    for (auto& p : *particles) {
      p.setOldF(p.getF());
      p.setF({});
    }

    // Calculate the forces acting on the particles
    for (const auto& force : forces) {
      force->calculateF();
    };

    // Calculate the velocities of the particles
    ForceCalc::calculateV(*particles, dt);

    iteration++;
    // Update temperature
    if (thermostat && (iteration % thermoFrequency == 0)) {
      thermostat->updateTemperature();
    }
    // Write state of particles to VTK
    if (writeFrequency > 0 && iteration % writeFrequency == 0) {
      plotParticles(iteration);
    }
    // Write state of particles to checkpoint file
    if (checkpointFrequency > 0 && iteration % checkpointFrequency == 0) {
      writeCheckpoint(iteration, currentTime);
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

  // Use member currentTime instead of local variable
  long iteration = startIteration;

  // For this loop, we assume: current positions, forces and velocities are known
  while (currentTime < endTime) {
    // Calculate the new position of the particles
    ForceCalc::calculateX(*particles, dt);

    // Store the force from the previous time stop for velocity update
    /* TODO: Optimization idea: introduce a 'first force' flag to ForceCalc to avoid
    * iterating an extra time over the particles, thus eliminating the loop below.
    */
    for (auto& p : *particles) {
      p.setOldF(p.getF());
      p.setF({});
    }

    // Calculate the forces acting on the particles
    for (auto& force : forces) {
      force->calculateF();
    };

    // Calculate the velocities of the particles
    ForceCalc::calculateV(*particles, dt);

    iteration++;
    // Update temperature
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

    size_t numParticles = particles->size();
    double mups = (static_cast<double>(iteration) * numParticles) / elapsed;
    SPDLOG_INFO("Molecule-Updates per Second (MUPS): {:.2f}", mups);

    // Write benchmark results to file
    std::string outputDir = "build/output";
    if (std::filesystem::exists("CMakeCache.txt")) {
      outputDir = "output";  // Assume we are in build/
    }

    std::filesystem::create_directories(outputDir);

    std::ofstream benchFile(outputDir + "/benchmark.txt", std::ios_base::app);
    if (benchFile.is_open()) {
      benchFile << "--------------------------------------------------" << std::endl;
      benchFile << "Benchmark Run: " << outputBasename << std::endl;
      benchFile << "Timestamp: " << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << std::endl;
      benchFile << "Particles: " << numParticles << std::endl;
      benchFile << "Time elapsed: " << elapsed << " s" << std::endl;
      benchFile << "Total iterations: " << iteration << std::endl;
      benchFile << "Mean time per iteration: " << timePerIteration << " s" << std::endl;
      benchFile << "MUPS: " << mups << std::endl;
      benchFile.close();
      SPDLOG_INFO("Benchmark results written to {}/benchmark.txt", outputDir);
    } else {
      SPDLOG_WARN("Could not open benchmark output file.");
    }
  }
  spdlog::set_level(spdlog::level::off);
}

void Simulation::setupSimulation() {
  // Generate particles
  particleGenerator->generate(*particles);

  // Precompute constants for all forces that support it (e.g., LennardJonesForce lookup tables)
  for (const auto& force : forces) {
    force->precomputeConstants();
  }

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
