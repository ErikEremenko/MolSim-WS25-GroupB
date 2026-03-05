#include "physics/ForceCalc.h"

ForceCalc::ForceCalc(ParticleContainer& particles) : particles(particles) {}
ForceCalc::~ForceCalc() = default;

void ForceCalc::calculateX(ParticleContainer& particles, const double dt) {
  const double dt_sq_half = 0.5 * dt * dt;
  const size_t n = particles.size();

  // Parallelize outer loop (each particle update is independent)
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    auto& p = particles[i];
    const auto& x_curr = p.getX();
    const auto& v = p.getV();
    const auto& F = p.getF();
    const double inv_m = 1.0 / p.getM();

    // compute all 3 components (SIMD friendly)
    std::array<double, 3> x_new;
#pragma omp simd
    for (int d = 0; d < 3; ++d) {
      x_new[d] = x_curr[d] + dt * v[d] + dt_sq_half * inv_m * F[d];
    }
    p.setX(x_new);
  }
}

void ForceCalc::calculateV(ParticleContainer& particles, const double dt) {
  const size_t n = particles.size();

  // Parallelize outer loop (each particle update is independent)
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    auto& p = particles[i];
    const auto& v_curr = p.getV();
    const auto& F = p.getF();
    const auto& F_old = p.getOldF();
    const double dt_half_inv_m = dt / (2.0 * p.getM());

    std::array<double, 3> v_new;
#pragma omp simd
    for (int d = 0; d < 3; ++d) {
      v_new[d] = v_curr[d] + dt_half_inv_m * (F_old[d] + F[d]);
    }
    p.setV(v_new);
  }
}

void ForceCalc::precomputeConstants() {}
