
#include <iostream>
#include <list>

#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

/**** forward declaration of the calculation functions ****/

/**
 * calculate the force for all particles
 */
void calculateF();

/**
 * calculate the position for all particles
 */
void calculateX();

/**
 * calculate the position for all particles
 */
void calculateV();

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

constexpr double start_time = 0;
constexpr double end_time = 1000;
constexpr double delta_t = 0.014;

// TODO: what data structure to pick?
std::list<Particle> particles;

int main(int argc, char *argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 2) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename" << std::endl;
  }

  FileReader fileReader;
  fileReader.readFile(particles, argsv[1]);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX();
    // calculate new f
    calculateF();
    // calculate new v
    calculateV();

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

void calculateF() {
  for (auto &p_i : particles) {
    std::array<double, 3> force_i = {0., 0., 0.};
    for (auto &p_j : particles) {
      if (p_i == p_j) {
        continue;
      }

      std::array<double, 3> dist;
      double norm_squared = 0.;
      for (int i = 0; i < 3; i++) {
        dist[i] = p_j.getX()[i] - p_i.getX()[i];
        norm_squared += dist[i] * dist[i];
      }
      if (norm_squared == 0) {
        // avoid division by zero when particles are at the same position
        continue;
      }
      double norm = std::sqrt(norm_squared);
      double norm_cubed = norm_squared * norm;
      for (int i = 0; i < 3; i++) {
        force_i[i] += (p_i.getM() * p_j.getM()) / norm_cubed * (dist[i]);
      }
    }
    p_i.set_f(force_i);
  }
}

void calculateX() {
  for (auto &p : particles) {
    // @TODO: insert calculation of position updates here!
  }
}

void calculateV() {
  for (auto &p : particles) {
    // @TODO: insert calculation of veclocity updates here!
  }
}

void plotParticles(int iteration) {
  std::string out_name("MD_vtk");

  outputWriter::XYZWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
