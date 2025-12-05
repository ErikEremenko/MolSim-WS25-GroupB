#include "Simulation.h"

#include "LinkedCellParticleContainer.h"
#include "io/FileReader.h"
#include "io/VTKWriter.h"

#include <chrono>

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif  // SPDLOG_ACTIVE_LEVEL
#include "spdlog/spdlog.h"

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

// Simulation run methods
void BaseSimulation::runFileOutput(int frequency, const std::string& outputBaseName) const {
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

void BaseSimulation::runBenchmark() const {
  using namespace std::chrono;
  // used for benchmark
  const auto chronoStart = steady_clock::now();

  // Benchmark begin
  constexpr double start_time = 0;
  double current_time = start_time;

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
  }

  const auto chronoEnd = steady_clock::now();
  const auto elapsed = duration_cast<duration<double>>(chronoEnd - chronoStart);
  spdlog::set_level(spdlog::level::info);
  SPDLOG_INFO("Time elapsed: {} s", elapsed.count());
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

// CollisionSimulation definitions
CollisionSimulation::CollisionSimulation(std::string inputFilename, double end_time, double dt,
                                         const SimulationMode simulationMode)
    : BaseSimulation(end_time, dt, 10, "MD_vtk", simulationMode), inputFilename(std::move(inputFilename)) {
  particles = std::make_unique<ParticleContainer>();
  constexpr double sigma = 1.0;
  constexpr double cutoffRadius = 2.5 * sigma;
  forceCalc = std::make_unique<LennardJonesForce>(*particles, 5.0, sigma, cutoffRadius);
}

void CollisionSimulation::setupSimulation() {
  CuboidFileReader reader(inputFilename);
  reader.readFile(*particles);
}

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

YAMLSimulation::YAMLSimulation(std::string inputFilename, const SimulationMode simulationMode, const ContainerKind kind,
                               const Parallelization parallelization)
    : BaseSimulation(simulationMode), inputFilename(std::move(inputFilename)), reader(this->inputFilename) {
  this->dt = reader.getDeltaT();
  this->end_time = reader.getTend();
  this->write_frequency = reader.getWriteFrequency();
  this->base_name = reader.getOutputBaseName();

  double epsilon = reader.getEpsilon();
  double sigma = reader.getSigma();
  double cutoffRadius = reader.getCutoff();
  const double repulsionDistance = reader.getLJRepulsionDistance();
  const auto domainSize = reader.getDomainSize();
  const auto boundariesRaw = reader.getBoundaryTypesRaw();

  if (kind == ContainerKind::DIRECT) {
    // legacy O(n^2) implementation
    particles = std::make_unique<ParticleContainer>();
    forceCalc = std::make_unique<LennardJonesForce>(*particles, epsilon, sigma, cutoffRadius, repulsionDistance);
    if (parallelization == Parallelization::ON) {
      // parallel direct sum LennardJones
      forceCalc = std::make_unique<LennardJonesForceParallel>(*particles, epsilon, sigma, cutoffRadius);
    } else {
      // serial direct sum LennardJones
      forceCalc = std::make_unique<LennardJonesForce>(*particles, epsilon, sigma, cutoffRadius, repulsionDistance);
    }
  } else {
    // linked cell implementation -> O(n)
    std::array<LinkedCellParticleContainer::BoundaryType, 6> boundaryTypes{};
    for (int i = 0; i < 6; ++i) {
      boundaryTypes[i] = parseBoundary(boundariesRaw[i]);
    }
    particles = std::make_unique<LinkedCellParticleContainer>(domainSize, cutoffRadius, boundaryTypes);

    //  serial LJ with linked cells + reflective boundaries
    forceCalc = std::make_unique<LennardJonesForce>(*particles, epsilon, sigma, cutoffRadius, repulsionDistance);
  }
}

void YAMLSimulation::setupSimulation() {
  reader.readFile(*particles);
  SPDLOG_INFO("YAML Simulation configured. dt={}, t_end={}, write_frequency={}, base_name={}", dt, end_time,
              write_frequency, base_name);
}
