#include "simulation/Simulation.h"

#include "physics/LinkedCellParticleContainer.h"
#include "io/FileReader.h"
#include "io/VTKWriter.h"
//benchmark
#include <chrono>
// process signal handling, TODO: Implement signal handling for all simulation types (?)
#include <atomic>
#include <csignal>

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif  // SPDLOG_ACTIVE_LEVEL
#include "spdlog/spdlog.h"

// Atomic flag for loop control
std::atomic<bool> simulation_running{true};

void sigint_handler(int signal) {
  if (signal == SIGINT) {
    simulation_running = false;
  }
}

namespace {
LinkedCellParticleContainer::BoundaryType parseBoundary(const std::string& s) {
  if (s == "OUTFLOW")
    return LinkedCellParticleContainer::BoundaryType::OUTFLOW;
  if (s == "REFLECTIVE")
    return LinkedCellParticleContainer::BoundaryType::REFLECTIVE;
  if (s == "PERIODIC")
    return LinkedCellParticleContainer::BoundaryType::PERIODIC;
  throw std::runtime_error("Unknown boundary type in YAML: " + s);
}
}  // namespace

// BaseSimulation below
BaseSimulation::BaseSimulation(double end_time, double dt, int write_frequency, const std::string& base_name,
                               SimulationMode simulationMode)
    : end_time(end_time),
      dt(dt),
      write_frequency(write_frequency),
      base_name(base_name),
      simulationMode(simulationMode) {}

BaseSimulation::BaseSimulation(SimulationMode simulationMode)
    : end_time(0.0), dt(0.0), write_frequency(10), base_name("MD_vtk"), simulationMode(simulationMode) {}

BaseSimulation::~BaseSimulation() = default;

void BaseSimulation::plotParticles(const int iteration, const std::string& outputBaseName) const {
  const std::string& out_name(outputBaseName);
  outputWriter::VTKWriter::plotParticles(*particles, out_name, iteration);
}

void BaseSimulation::runFileOutput(int frequency, const std::string& outputBaseName) {
  // TODO: Do the following checks when reading the config file, not here!
  if (frequency < 1) {
    SPDLOG_INFO("Write frequency must be a positive integer, but was given {}", frequency);
    throw std::invalid_argument("Write frequency must be a positive integer");
  } else if (outputBaseName.empty()) {
    SPDLOG_INFO("Write frequency must be a positive integer, but was given an empty string");
    throw std::invalid_argument("Base name must not be an empty string");
  }
  constexpr double start_time = 0;

  double current_time = start_time;
  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
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
    if (iteration % frequency == 0) {
      plotParticles(iteration, outputBaseName);
    }
    current_time += dt;
  }
}

void BaseSimulation::runBenchmark() {
  std::signal(SIGINT, sigint_handler);
  simulation_running = true;
  using namespace std::chrono;
  // Used for benchmarks
  const auto chronoStart = steady_clock::now();

  // Benchmark begin
  constexpr double start_time = 0;
  double current_time = start_time;
  long iteration = 0;

  // For this loop, we assume current x, current F and current v are known
  while (current_time < end_time) {
    // Calculate the new positions of the particles
    forceCalc->calculateX(dt);
    for (auto& p : *particles) {
      p.setOldF(p.getF());  // store F(t_n) for v update
    }
    // Calculate new forces and velocities
    forceCalc->calculateF();
    forceCalc->calculateV(dt);

    current_time += dt;
    iteration++;
  }

  const auto chronoEnd = steady_clock::now();
  const auto elapsed = duration_cast<duration<double>>(chronoEnd - chronoStart).count();

  std::signal(SIGINT, SIG_DFL);
  spdlog::set_level(spdlog::level::info);
  if (!simulation_running) {
    SPDLOG_INFO("Simulation stopped early by user (SIGINT).");
  } else {
    SPDLOG_INFO("Benchmark finished normally.");
  }
  SPDLOG_INFO("Time elapsed: {} s", elapsed);
  SPDLOG_INFO("Total iterations: {}", iteration);
  if (iteration > 0) {
    double time_per_iter = elapsed / static_cast<double>(iteration);
    SPDLOG_INFO("Mean time per iteration: {:.6f} s", time_per_iter);
  }
  spdlog::set_level(spdlog::level::off);
}

void BaseSimulation::run() {
  setupSimulation();

  if (simulationMode == SimulationMode::FILE_OUTPUT) {
    runFileOutput(write_frequency, base_name);
  } else if (simulationMode == SimulationMode::BENCHMARK) {
    runBenchmark();
  }
}

// CollisionSimulation below
CollisionSimulation::CollisionSimulation(std::string inputFilename, double end_time, double dt,
                                         const SimulationMode simulationMode)
    : BaseSimulation(end_time, dt, 10, "MD_vtk", simulationMode), inputFilename(std::move(inputFilename)) {
  particles = std::make_unique<ParticleContainer>();
  constexpr double sigma = 1.0;
  constexpr double cutoffRadius = 2.5 * sigma;
  forceCalc = std::make_unique<LennardJonesForce>(*particles, 5.0, sigma, cutoffRadius, 0);
}

void CollisionSimulation::setupSimulation() {
  CuboidFileReader reader(inputFilename);
  reader.readFile(*particles);
}

// CollisionSimulationParallel below
CollisionSimulationParallel::CollisionSimulationParallel(std::string inputFilename, double end_time, double dt,
                                                         const SimulationMode simulationMode)
    : BaseSimulation(end_time, dt, 10, "MD_vtk", simulationMode), inputFilename(std::move(inputFilename)) {
  particles = std::make_unique<ParticleContainer>();
  constexpr double sigma = 1.0;
  constexpr double cutoffRadius = 2.5 * sigma;
  forceCalc = std::make_unique<LennardJonesForceParallel>(*particles, 5.0, sigma, cutoffRadius);
}

void CollisionSimulationParallel::setupSimulation() {
  CuboidFileReader reader(inputFilename);
  reader.readFile(*particles);
}

// YAMLSimulation below
YAMLSimulation::YAMLSimulation(std::string inputFilename, const SimulationMode simulationMode, const ContainerKind kind,
                               const Parallelization parallelization)
    : BaseSimulation(simulationMode), inputFilename(std::move(inputFilename)), reader(this->inputFilename) {
  this->dt = reader.getDeltaT();
  this->end_time = reader.getTend();
  this->write_frequency = reader.getWriteFrequency();
  this->base_name = reader.getOutputBaseName();

  // imulation parameters for checkpoint generation
  epsilon = reader.getEpsilon();
  sigma = reader.getSigma();
  cutoffRadius = reader.getCutoff();
  domainSize = reader.getDomainSize();
  boundaryTypes = reader.getBoundaryTypesRaw();
  checkpointFrequency = reader.getCheckpointFrequency();

  double gravity = reader.getGravity();

  // Checkpoint resume state (0 if not a checkpoint file)
  startIteration = reader.getCheckpointIteration();
  startTime = reader.getCheckpointTime();

  if (kind == ContainerKind::DIRECT) {
    // Legacy O(n^2) implementation
    particles = std::make_unique<ParticleContainer>();
    forceCalc = std::make_unique<LennardJonesForce>(*particles, epsilon, sigma, cutoffRadius, gravity);
    if (parallelization == Parallelization::ON) {
      // Parallel direct sum LennardJones
      forceCalc = std::make_unique<LennardJonesForceParallel>(*particles, epsilon, sigma, cutoffRadius);
    } else {
      // Serial direct sum LennardJones
      forceCalc = std::make_unique<LennardJonesForce>(*particles, epsilon, sigma, cutoffRadius, gravity);
    }
  } else {
    // Linked cell implementation -> O(n)
    std::array<LinkedCellParticleContainer::BoundaryType, 6> boundaryTypesEnum{};
    for (int i = 0; i < 6; ++i) {
      boundaryTypesEnum[i] = parseBoundary(boundaryTypes[i]);
    }
    particles = std::make_unique<LinkedCellParticleContainer>(domainSize, cutoffRadius, boundaryTypesEnum);
    forceCalc = std::make_unique<LennardJonesForce>(*particles, epsilon, sigma, cutoffRadius, gravity);
  }
}

void YAMLSimulation::setupSimulation() {
  reader.readFile(*particles);
  SPDLOG_INFO("YAML Simulation configured. dt={}, t_end={}, write_frequency={}, base_name={}", dt, end_time,
              write_frequency, base_name);
}

void YAMLSimulation::writeCheckpoint(int iteration, double currentTime) const {
  std::string checkpointFilename = base_name + "_checkpoint_" + std::to_string(iteration) + ".yaml";
  outputWriter::CheckpointWriter::writeCheckpoint(*particles, checkpointFilename, iteration, currentTime, base_name,
                                                  write_frequency, checkpointFrequency, end_time, dt, epsilon, sigma,
                                                  cutoffRadius, domainSize, boundaryTypes);
}

void YAMLSimulation::runFileOutput(int frequency, const std::string& outputBaseName) {
  if (frequency < 1) {
    SPDLOG_INFO("Write frequency must be a positive integer, but was given {}", frequency);
    throw std::invalid_argument("Write frequency must be a positive integer");
  } else if (outputBaseName.empty()) {
    SPDLOG_INFO("Base name must not be an empty string");
    throw std::invalid_argument("Base name must not be an empty string");
  }

  // Use checkpoint values if resuming, otherwise start from 0
  double current_time = startTime;
  int iteration = startIteration;

  if (startIteration > 0) {
    SPDLOG_INFO("Resuming simulation from checkpoint: iteration={}, time={}", startIteration, startTime);
  }

  // For this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    forceCalc->calculateX(dt);
    for (auto& p : *particles) {
      p.setOldF(p.getF());  // Store f(t_n) for v update (Störmer-Verlet)
    }
    forceCalc->calculateF();
    forceCalc->calculateV(dt);

    iteration++;
    current_time += dt;

    if (iteration % frequency == 0) {
      plotParticles(iteration, outputBaseName);
    }
    if (checkpointFrequency > 0 && iteration % checkpointFrequency == 0) {
      writeCheckpoint(iteration, current_time);
    }
  }

  // Write final checkpoint at end of simulation
  if (checkpointFrequency > 0) {
    writeCheckpoint(iteration, current_time);
    SPDLOG_INFO("Final checkpoint written at iteration {}", iteration);
  }
}

// YAMLThermostatSimulation below
void YAMLThermostatSimulation::runBenchmark() {
  constexpr double start_time = 0;

  double current_time = start_time;
  int iteration = 0;
  const int thermostatFrequency = thermostat.getUpdateFrequency();

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {  // TODO: Refactor these loops (also in runBenchmark)
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
    if (iteration % thermostatFrequency == 0) {  // TODO: Optimize this if check
      thermostat.updateTemperature();
    }
    current_time += dt;
  }
}

void YAMLThermostatSimulation::runFileOutput(int frequency, const std::string& outputBaseName) {
  if (frequency < 1) {
    SPDLOG_ERROR("Write frequency must be a positive integer, but was given {}", frequency);
    throw std::invalid_argument("Write frequency must be a positive integer");
  } else if (outputBaseName.empty()) {
    SPDLOG_ERROR("Base name must not be an empty string");
    throw std::invalid_argument("Base name must not be an empty string");
  }

  // Use checkpoint values if resuming, otherwise start from 0
  double current_time = startTime;
  int iteration = startIteration;
  const int thermostatFrequency = thermostat.getUpdateFrequency();

  if (startIteration > 0) {
    SPDLOG_INFO("Resuming thermostat simulation from checkpoint: iteration={}, time={}", startIteration, startTime);
  }

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {  // TODO: Refactor these loops (also in runBenchmark)
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
    current_time += dt;

    if (iteration % frequency == 0) {
      plotParticles(iteration, outputBaseName);
    }
    if (iteration % thermostatFrequency == 0) {  // TODO: Optimize this if check
      thermostat.updateTemperature();
    }
    if (checkpointFrequency > 0 && iteration % checkpointFrequency == 0) {
      writeCheckpoint(iteration, current_time);
    }
  }

  // Write final checkpoint at end of simulation
  if (checkpointFrequency > 0) {
    writeCheckpoint(iteration, current_time);
    SPDLOG_INFO("Final checkpoint written at iteration {}", iteration);
  }
}

YAMLThermostatSimulation::YAMLThermostatSimulation(std::string inputFilename, SimulationMode simulationMode,
                                                   const Thermostat& thermostat, ContainerKind kind,
                                                   Parallelization parallelization)
    : YAMLSimulation(std::move(inputFilename), simulationMode, kind, parallelization), thermostat(thermostat) {}
