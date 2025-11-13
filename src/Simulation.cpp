#include "Simulation.h"

#include "io/FileReader.h"
#include "io/VTKWriter.h"

#include <chrono>
#include <iostream>

BaseSimulation::BaseSimulation(double end_time, double dt, SimulationMode simulationMode)
  : end_time(end_time),
    dt(dt),
    simulationMode(simulationMode) {

}
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
  auto chronoStart = high_resolution_clock::now();

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

  auto chronoEnd = high_resolution_clock::now();
  auto elapsed = chronoEnd - chronoStart;
  std::cout << "Time elapsed: " << elapsed.count() << " s\n";
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
CollisionSimulation::CollisionSimulation(
  std::string inputFilename, double end_time,
  double dt, SimulationMode simulationMode
  )
    : BaseSimulation(end_time, dt, simulationMode), inputFilename(std::move(inputFilename))
{
  particles = std::make_unique<ParticleContainer>();
  forceCalc = std::make_unique<LennardJonesForce>(*particles, 5.0, 1.0);
}

void CollisionSimulation::setupSimulation() {
  FileCuboidReader reader(inputFilename);
  reader.readFile(*particles);
}
