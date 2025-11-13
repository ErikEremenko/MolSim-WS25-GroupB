#include "Simulation.h"

#include <iostream>

int main(const int argc, char* argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  // Read arguments from the command line
  if (argc != 5) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename t_end delta_t [file | benchmark]" << std::endl;
    return 1;
  }

  SimulationMode simulation_mode;
  std::string mode_arg = argsv[4];
  if (mode_arg == "file") {
    simulation_mode = SimulationMode::FILE_OUTPUT;
  } else if (mode_arg == "benchmark") {
    simulation_mode = SimulationMode::BENCHMARK;
  } else {
    std::cout << "Invalid simulation mode provided: " << mode_arg << std::endl;
    std::cout << "Please use 'file' or 'benchmark'" << std::endl;
    return 1;
  }
  CollisionSimulation simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), simulation_mode);
  simulation.run();

  // /*
  // std::cout << "output written. Terminating..." << std::endl;
  // */

  return 0;
}