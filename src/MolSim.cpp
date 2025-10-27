
#include <float.h>

#include <iostream>
#include <list>

#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "utils/ArrayUtils.h"

#define EPSILON 1e-12

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
std::vector<Particle> particles;

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

  for (auto &p : particles) {
    p.set_f({0., 0., 0.});
  }

  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      Particle &p_i = particles[i];
      Particle &p_j = particles[j];

      std::array<double, 3> dist;
      double norm_squared = 0.;
      for (int k= 0; k < 3; ++k) {
        dist[k] = p_j.getX()[k] - p_i.getX()[k];
        norm_squared += dist[k] * dist[k];
      }
      if (norm_squared < EPSILON) {
        // avoid division by zero and numerical explosion when particles are at the same position or very close
        continue;
      }
      double norm = std::sqrt(norm_squared);
      double norm_cubed = norm_squared * norm;
      std::array<double, 3> f_vector;

      for (int k = 0; k < 3; ++k) {
        f_vector[k] = (p_i.getM() * p_j.getM()) / norm_cubed * dist[k];
      }

      // apply forces using Newton's third law

      auto force_i = p_i.getF();
      auto force_j = p_j.getF();

      for (int k = 0; k < 3; ++k) {
        // actio est reactio
        force_i[k] += f_vector[k];
        force_j[k] -= f_vector[k];
      }
      p_i.set_f(force_i);
      p_j.set_f(force_j);
    }
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

void plotParticles(const int iteration) {
  const std::string out_name("MD_vtk");
  outputWriter::XYZWriter writer;
  outputWriter::XYZWriter::plotParticles(particles, out_name, iteration);
}
