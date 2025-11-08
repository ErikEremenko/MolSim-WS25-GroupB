#include "CalcMethod.h"
#include <omp.h>
#include "utils/ArrayUtils.h"

CalcMethod::~CalcMethod() = default;

void StormerVerletMethod::calculateX(const double dt) {
#pragma omp parallel for
  for (auto& p : particles) {
    const auto x_curr = p.getX();
    const double m = p.getM();
    const auto F = p.getF();
    const auto v = p.getV();
    const auto a = (1.0 / m) * F;
    const auto x_new = x_curr + dt * v + 0.5 * (dt * dt) * a;
    p.setX(x_new);
  }
}

void StormerVerletMethod::calculateV(const double dt) {
#pragma omp parallel for
  for (auto& p : particles) {
    const auto v_curr = p.getV();
    const double m_i = p.getM();
    const auto F = p.getF();
    const auto F_old = p.getOldF();
    const auto v_new = v_curr + (dt / (2.0 * m_i)) * (F_old + F);
    p.setV(v_new);
  }
}

void StormerVerletMethod::calculateGravityF() {
  for (auto& p : particles) {
    p.setF({});
  }

  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      if (norm == 0) {
        // avoid division by zero
        throw std::overflow_error(
            "Calculated a zero norm between particles. This is likely caused "
            "by an incorrect initialization of the Simulation.");
      }
      const double norm3 = norm * norm * norm;

      const auto F_vector = ((p_i.getM() * p_j.getM()) / norm3) * dist;

      // apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
      auto F_i = p_i.getF();
      auto F_j = p_j.getF();
      // actio est reactio
      p_i.setF(F_i + F_vector);
      p_j.setF(F_j - F_vector);
    }
  }
}

void StormerVerletMethod::calculateLennardJonesF(const double epsilon,
                                                 const double sigma) {
  const size_t n_particles = particles.size();
  const double sigma2 = sigma * sigma;
  const double sigma6 = sigma2 * sigma2 * sigma2;

#pragma omp parallel
  {

#pragma omp for
    for (auto& p : particles) {
      p.setF({});
    }
    // accumulation buffer
    std::vector<std::array<double, 3>> F_buf(n_particles, {0.0, 0.0, 0.0});
#pragma omp for schedule(guided)
    for (size_t i = 0; i < n_particles; ++i) {
      for (size_t j = i + 1; j < n_particles; ++j) {
        auto& p_i = particles[i];
        auto& p_j = particles[j];

        const auto dist = p_j.getX() - p_i.getX();
        const double norm = ArrayUtils::L2Norm(dist);
        if (norm == 0) {
          // avoid division by zero
          throw std::overflow_error(
              "Calculated a zero norm between particles. This is likely caused "
              "by an incorrect initialization of the Simulation.");
        }
        const double inv_norm2 = 1.0 / (norm * norm);
        const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

        const double crossing_norm_quot_6 = sigma6 * inv_norm6;
        const double crossing_norm_quot_12 =
            crossing_norm_quot_6 * crossing_norm_quot_6;

        const auto F_vector =
            (24.0 * epsilon) * inv_norm2 *
            (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12) * dist;
        F_buf[i] = F_buf[i] + F_vector;
        F_buf[j] = F_buf[j] - F_vector;
      }
    }

    // load buffered forces into particles
#pragma omp for
    for (size_t i = 0; i < n_particles; ++i) {
      particles[i].setF(particles[i].getF() + F_buf[i]);
    }
  }
}