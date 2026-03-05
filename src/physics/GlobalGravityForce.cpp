#include "physics/GlobalGravityForce.h"

#include <spdlog/spdlog.h>

GlobalGravityForce::GlobalGravityForce(ParticleContainer& particles, double g, int axis)
    : ForceCalc(particles), gravity(g), axis(axis) {
  if (axis < 0 || axis > 2) {
    SPDLOG_WARN("Invalid gravity axis {}, defaulting to y-axis (1)", axis);
    this->axis = 1;
  }
}

void GlobalGravityForce::calculateF() {
  const size_t n = particles.size();
  const double g = gravity;
  const int ax = axis;

#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    Particle& p = particles[i];
    auto F = p.getF();
    F[ax] += p.getM() * g;
    p.setF(F);
  }
}
