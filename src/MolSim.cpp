
#include <float.h>

#include <iostream>
#include <list>

#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
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
void calculateX(double delta_t);

/**
 * calculate the position for all particles
 */
void calculateV(double delta_t);

/**
 * plot the particles to a xyz-file
 */
void plotParticles(int iteration);

std::vector<Particle> particles;

int main(int argc, char* argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 4) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename t_end delta_t" << std::endl;
    return 1;
  }

  char* filename = argsv[1];
  const double end_time = std::stod(argsv[2]);
  const double delta_t = std::stod(argsv[3]);
  constexpr double start_time = 0;

  FileReader fileReader;
  fileReader.readFile(particles, filename);

  double current_time = start_time;

  int iteration = 0;

  // for this loop, we assume: current x, current f and current v are known
  while (current_time < end_time) {
    // calculate new x
    calculateX(delta_t);
    for (auto& p : particles) {
      p.setOldF(p.getF());  // store f(t_n) for v update
    }
    // calculate new f
    calculateF();
    // calculate new v
    calculateV(delta_t);

    iteration++;
    if (iteration % 10 == 0) {
      plotParticles(iteration);
    }
    // std::cout << "Iteration " << iteration << " finished." << std::endl;

    current_time += delta_t;
  }

  std::cout << "output written. Terminating..." << std::endl;
  return 0;
}

void calculateF() {
  for (auto& p : particles) {
    p.setF({});
  }

  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      constexpr double norm_squared_soft_const = 1e-9;
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      if (norm < EPSILON) {
        // avoid division by zero and numerical explosion when particles are at the same position or very close
        continue;
      }
      const double norm_squared_softened =
          norm * norm + norm_squared_soft_const;
      const double norm_cubed_softened = norm_squared_softened * norm;

      const auto F_vector =
          ((p_i.getM() * p_j.getM()) / norm_cubed_softened) * dist;

      // apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
      auto F_i = p_i.getF();
      auto F_j = p_j.getF();
      // actio est reactio
      F_i = F_i + F_vector;
      F_j = F_j - F_vector;
      p_i.setF(F_i);
      p_j.setF(F_j);
    }
  }
}

void calculateX(const double delta_t) {
  for (auto& p : particles) {
    const auto x_curr = p.getX();
    const double m = p.getM();
    std::array<double, 3> F = p.getF();
    std::array<double, 3> v = p.getV();
    const auto a = (1.0 / m) * F;
    const auto x_new = x_curr + delta_t * v + 0.5 * (delta_t * delta_t) * a;
    p.setX(x_new);
  }
}

void calculateV(const double delta_t) {
  for (auto& p : particles) {
    const auto v_curr = p.getV();
    const double m_i = p.getM();
    const auto F = p.getF();
    const auto F_old = p.getOldF();
    const auto v_new = v_curr + (delta_t / (2.0 * m_i)) * (F_old + F);
    p.setV(v_new);
  }
}

void plotParticles(const int iteration) {
  const std::string out_name("MD_vtk");
  /*outputWriter::XYZWriter writer;
  outputWriter::XYZWriter::plotParticles(particles, out_name, iteration);*/
  outputWriter::VTKWriter writer;
  writer.plotParticles(particles, out_name, iteration);
}
