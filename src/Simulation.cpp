#include "Simulation.h"

#include "io/FileReader.h"
#include "io/VTKWriter.h"

#include <chrono>
#include <iostream>

#include "spdlog/spdlog.h"

BaseSimulation::BaseSimulation(double end_time, double dt, SimulationMode simulationMode)
    : end_time(end_time), dt(dt), simulationMode(simulationMode) {}
BaseSimulation::~BaseSimulation() = default;

void BaseSimulation::plotParticles(const int iteration) const {
  const std::string out_name("MD_vtk");
  outputWriter::VTKWriter::plotParticles(*particles, out_name, iteration);
}

// Simulation run methods
void BaseSimulation::runFileOutput() const {
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
    if (iteration % 10 == 0) {
      plotParticles(iteration);
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
    runFileOutput();
  } else if (simulationMode == SimulationMode::BENCHMARK) {
    runBenchmark();
  }
}

// CollisionSimulation definitions
CollisionSimulation::CollisionSimulation(std::string inputFilename, double end_time, double dt,
                                         const SimulationMode simulationMode)
    : BaseSimulation(end_time, dt, simulationMode), inputFilename(std::move(inputFilename)) {
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
    : BaseSimulation(end_time, dt, simulationMode), inputFilename(std::move(inputFilename)) {
  particles = std::make_unique<ParticleContainer>();
  constexpr double sigma = 1.0;
  constexpr double cutoffRadius = 2.5 * sigma;
  forceCalc = std::make_unique<LennardJonesForceParallel>(*particles, 5.0, sigma, cutoffRadius);
}

void CollisionSimulationParallel::setupSimulation() {
  CuboidFileReader reader(inputFilename);
  reader.readFile(*particles);
}