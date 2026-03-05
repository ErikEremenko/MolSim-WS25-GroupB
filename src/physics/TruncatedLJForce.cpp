#include "physics/TruncatedLJForce.h"

#include <cmath>

#include "physics/LinkedCellParticleContainer.h"
#include "utils/PhysicsConstants.h"

using namespace PhysicsConstants;

TruncatedLJForce::TruncatedLJForce(ParticleContainer& particles) : ForceCalc(particles) {}

void TruncatedLJForce::calculateF() {
  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);

  auto applyTruncatedLJ = [](Particle& p_i, Particle& p_j) {
    const auto& xi = p_i.getX();
    const auto& xj = p_j.getX();

    const double dx = xj[0] - xi[0];
    const double dy = xj[1] - xi[1];
    const double dz = xj[2] - xi[2];
    const double distSq = dx * dx + dy * dy + dz * dz;

    // Lorentz-Berthelot mixing rules
    const double sigma_ij = (p_i.getSigma() + p_j.getSigma()) * 0.5;
    const double sigma2 = sigma_ij * sigma_ij;
    // Cutoff at 2^(1/6) * sigma (repulsion only): cutoffSq = 2^(1/3) * sigma^2

    if (const double cutoffSq = TWO_POW_1_3 * sigma2; distSq > 0.0 && distSq < cutoffSq) {
      const double epsilon_ij = std::sqrt(p_i.getEpsilon() * p_j.getEpsilon());
      const double sigma6 = sigma2 * sigma2 * sigma2;
      const double sigma12 = sigma6 * sigma6;

      const double inv_distSq = 1.0 / distSq;
      const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
      const double term = 48.0 * epsilon_ij * sigma12 * inv_distSq * inv_distSq3 * (0.5 / sigma6 - inv_distSq3);

      const double fx = term * dx;
      const double fy = term * dy;
      const double fz = term * dz;

      auto& fi = p_i.getF();
      auto& fj = p_j.getF();
      fi[0] += fx;
      fi[1] += fy;
      fi[2] += fz;
      fj[0] -= fx;
      fj[1] -= fy;
      fj[2] -= fz;
    }
  };

  if (lc) {
    lc->handleOutflowBoundaries();
    lc->iteratePairsTemplate(applyTruncatedLJ);
  } else {
    const size_t n_particles = particles.size();
    for (size_t i = 0; i < n_particles; ++i) {
      for (size_t j = i + 1; j < n_particles; ++j) {
        applyTruncatedLJ(particles[i], particles[j]);
      }
    }
  }
}
