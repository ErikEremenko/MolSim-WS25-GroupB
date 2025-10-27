
#include <float.h>

#include <iostream>
#include <list>

#include "FileReader.h"
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
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

int main(int argc, char *argsv[]) {
  std::cout << "Hello from MolSim for PSE!" << std::endl;
  if (argc != 4) {
    std::cout << "Erroneous programme call! " << std::endl;
    std::cout << "./molsym filename t_end delta_t" << std::endl;
    return 1;
  }

  char *filename = argsv[1];
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
    for (auto &p : particles) {
      p.set_old_f(p.getF()); // store f(t_n) for v update
    }
    // calculate new f
    calculateF();
    // calculate new v
    calculateV(delta_t);

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
      std::array<double, 3> F_vector;

      for (int k = 0; k < 3; ++k) {
        F_vector[k] = (p_i.getM() * p_j.getM()) / norm_cubed * dist[k];
      }

      // apply forces using Newton's third law

      auto force_i = p_i.getF();
      auto force_j = p_j.getF();

      for (int k = 0; k < 3; ++k) {
        // actio est reactio
        force_i[k] += F_vector[k];
        force_j[k] -= F_vector[k];
      }
      p_i.set_f(force_i);
      p_j.set_f(force_j);
    }
  }
}

void calculateX(double delta_t) {
  for (auto &p : particles) {
    std::array<double, 3> x_curr = p.getX();
    std::array<double, 3> x_new;
    double m_i = p.getM();
    std::array<double, 3> F_i = p.getF();
    std::array<double, 3> v_i = p.getV();
    std::array<double, 3> a_i;

    // calculate component-wise acceleration and coordinates
    for (int j = 0; j < 3; ++j) {
      a_i[j] = F_i[j] / (2 * m_i);
    }
    for (int j = 0; j < 3; ++j) {
      x_new[j] = x_curr[j] + delta_t * v_i[j] + (delta_t * delta_t) * a_i[j];
    }
    p.set_x(x_new);
  }
}

void calculateV(double delta_t) {
  for (auto &p : particles) {
    std::array<double, 3> v_curr = p.getV();
    std::array<double, 3> v_new;
    const double m_i = p.getM();
    std::array<double, 3> F_i = p.getF();
    std::array<double, 3> F_i_old = p.getOldF();

    for (int j = 0; j < 3; ++j) {
      v_new[j] = v_curr[j] + delta_t * (F_i_old[j] + F_i[j]) / (2*m_i);
    }
    p.set_v(v_new);
  }
}

void plotParticles(const int iteration) {
  const std::string out_name("MD_vtk");
  /*outputWriter::XYZWriter writer;
  outputWriter::XYZWriter::plotParticles(particles, out_name, iteration);*/
  outputWriter::VTKWriter writer;
  writer.plotParticles(particles, out_name, iteration);

}
