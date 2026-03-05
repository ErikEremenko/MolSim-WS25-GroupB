#include "physics/HarmonicMembraneForce.h"

#include <cmath>

#include "utils/PhysicsConstants.h"

using namespace PhysicsConstants;

HarmonicMembraneForce::HarmonicMembraneForce(ParticleContainer& particles, double k, double r0)
    : ForceCalc(particles), stiffness(k), avgBondLength(r0), diagonalBondLength(SQRT_2 * r0) {}

void HarmonicMembraneForce::calculateF() {
  const size_t numParticles = particles.size();
  const double k = stiffness;
  const double r0 = avgBondLength;
  const double r0_diag = diagonalBondLength;

  for (size_t i = 0; i < numParticles; ++i) {
    Particle& p = particles[i];
    const auto& px = p.getX();
    auto& pf = p.getF();

    // Process direct neighbors (bond length = r0)
    for (const int neighborID : p.getDirectNeighbors()) {
      if (neighborID < 0 || neighborID >= static_cast<int>(numParticles))
        continue;
      const Particle& neighbor = particles[neighborID];
      const auto& nx = neighbor.getX();

      const double dx = nx[0] - px[0];
      const double dy = nx[1] - px[1];
      const double dz = nx[2] - px[2];
      const double distSq = dx * dx + dy * dy + dz * dz;

      if (distSq > 0.0) {
        const double inv_norm = 1.0 / std::sqrt(distSq);
        const double factor = k * (1.0 - r0 * inv_norm);
        pf[0] += factor * dx;
        pf[1] += factor * dy;
        pf[2] += factor * dz;
      }
    }

    // Process diagonal neighbors (bond length = sqrt(2) * r0)
    for (const int neighborID : p.getDiagonalNeighbors()) {
      if (neighborID < 0 || neighborID >= static_cast<int>(numParticles))
        continue;
      const Particle& neighbor = particles[neighborID];
      const auto& nx = neighbor.getX();

      const double dx = nx[0] - px[0];
      const double dy = nx[1] - px[1];
      const double dz = nx[2] - px[2];
      const double distSq = dx * dx + dy * dy + dz * dz;

      if (distSq > 0.0) {
        const double inv_norm = 1.0 / std::sqrt(distSq);
        const double factor = k * (1.0 - r0_diag * inv_norm);
        pf[0] += factor * dx;
        pf[1] += factor * dy;
        pf[2] += factor * dz;
      }
    }
  }
}
