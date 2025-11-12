#include "Simulation.h"

Simulation::Simulation(char* filename, const double end_time, const double dt,
             ParticleContainer& particles, ForceCalc& calcMethod)
      : filename(filename),
        end_time(end_time),
        dt(dt),
        particles(particles),
        calcMethod(calcMethod) {}
Simulation::~Simulation() = default;

void Simulation::plotParticles(const int iteration) const {
  const std::string out_name("MD_vtk");
  outputWriter::VTKWriter::plotParticles(particles, out_name, iteration);
}

void Simulation::loadParticles() const {
  FileCuboidReader::readFile(particles, filename);
}

void Simulation::run() const {
  constexpr double start_time = 0;

  // Load particles from file

  double current_time = start_time;
  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calcMethod.calculateX(dt);
    for (auto& p : particles) {
      p.setOldF(p.getF());  // store f(t_n) for v update
    }
    // calculate new f
    calcMethod.calculateF();
    // calculate new v
    calcMethod.calculateV(dt);

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    current_time += dt;
  }
}
