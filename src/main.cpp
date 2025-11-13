#include "Simulation.h"

#include <iostream>

int main(const int argc, char* argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  // Read arguments from the command line
  if (argc != 4) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename t_end delta_t" << std::endl;
    return 1;
  }

  ParticleContainer particles;
  LennardJonesForce ljf(particles, 5.0, 1.0);

  const Simulation simulation(argsv[1], std::stod(argsv[2]),
                              std::stod(argsv[3]), particles, ljf,
                              SimulationMode::FILE_OUTPUT);
  simulation.loadParticles();
  simulation.run();

  std::cout << "output written. Terminating..." << std::endl;

  return 0;
}