
#include <float.h>

#include <iostream>
#include <list>


#include "utils/ArrayUtils.h"

#include "Simulation.h"

int main(int argc, char* argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  // Read arguments from the command line
  if (argc != 4) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename t_end delta_t" << std::endl;
    return 1;
  }

  ParticleContainer particles;
  StormerVerletMethod svm(particles);

  Simulation simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), particles, svm);
  simulation.loadParticles();
  simulation.run();
  
  std::cout << "output written. Terminating..." << std::endl;
  
  return 0;
}
