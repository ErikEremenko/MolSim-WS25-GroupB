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

  CollisionSimulation simulation(argsv[1], std::stod(argsv[2]),
                                 std::stod(argsv[3]),
                                 SimulationMode::FILE_OUTPUT);
  simulation.run();

  std::cout << "output written. Terminating..." << std::endl;

  return 0;
}