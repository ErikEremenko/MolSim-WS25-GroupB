#include "physics/SmoothedLJForce.h"

#include <spdlog/spdlog.h>
#include <cmath>

#include "physics/LinkedCellParticleContainer.h"

SmoothedLJForce::SmoothedLJForce(ParticleContainer& particles, const double epsilon, const double sigma,
                                 double cutoffRadius, double smoothingRadius)
    : ForceCalc(particles),
      epsilon(epsilon),
      sigma(sigma),
      cutoffRadius(cutoffRadius),
      smoothingRadius(smoothingRadius),
      cutoffRadiusSq(cutoffRadius * cutoffRadius),
      smoothingRadiusSq(smoothingRadius * smoothingRadius),
      rcMinusRlCubed((cutoffRadius - smoothingRadius) * (cutoffRadius - smoothingRadius) *
                     (cutoffRadius - smoothingRadius)) {
  if (smoothingRadius >= cutoffRadius) {
    SPDLOG_ERROR("SmoothedLJForce: r_l ({}) must be less than r_c ({})", smoothingRadius, cutoffRadius);
    throw std::invalid_argument("Smoothing radius must be less than the cutoff radius!");
  }
}

void SmoothedLJForce::calculateF() {
  if (tableWidth == 0 || pairEpsilon.empty()) {
    precomputeConstants();
  }
  if (dynamic_cast<LinkedCellParticleContainer*>(&particles)) {
    calculateFLinkedCell();
  } else {
    calculateFDirectSum();
  }
}

void SmoothedLJForce::calculateFDirectSum() const {
  // Reset forces
  for (auto& p : particles) {
    p.setF({});
  }

  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    for (size_t j = i + 1; j < n_particles; ++j) {
      applySmoothedPairForce(particles[i], particles[j]);
    }
  }
}

void SmoothedLJForce::calculateFLinkedCell() {
  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);
  if (!lc) {
    throw std::runtime_error("SmoothedLJForce::calculateFLinkedCell requires LinkedCellParticleContainer");
  }
  lc->handleOutflowBoundaries();

  // Cache member variables for better optimization (avoid repeated this-> access)
  const double cutoffSq = cutoffRadiusSq;
  const double smoothingSq = smoothingRadiusSq;
  const double smoothingR = smoothingRadius;
  const double cutoffR = cutoffRadius;
  const double rcRlCubed = rcMinusRlCubed;
  const double* __restrict__ eps = pairEpsilon.data();
  const double* __restrict__ sig6 = pairSigma6.data();
  const double* __restrict__ sig12 = pairSigma12.data();
  const int tw = tableWidth;

  // Keep force calculation in the lambda for inlining.
  lc->iteratePairsTemplate(
      [cutoffSq, smoothingSq, smoothingR, cutoffR, rcRlCubed, eps, sig6, sig12, tw](Particle& p_i, Particle& p_j) {
        const auto& xi = p_i.getX();
        const auto& xj = p_j.getX();

        const double dx = xj[0] - xi[0];
        const double dy = xj[1] - xi[1];
        const double dz = xj[2] - xi[2];
        const double distSq = dx * dx + dy * dy + dz * dz;

        if (distSq >= cutoffSq) {
          return;
        }

        // Use precomputed lookup tables
        const int idx = p_i.getType() * tw + p_j.getType();
        const double epsilon_ij = eps[idx];
        const double sigma6 = sig6[idx];
        const double sigma12 = sig12[idx];

        const double inv_distSq = 1.0 / distSq;
        const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
        const double sigma6_d6 = sigma6 * inv_distSq3;
        const double sigma12_d12 = sigma12 * inv_distSq3 * inv_distSq3;

        const double U_LJ = 4.0 * epsilon_ij * (sigma12_d12 - sigma6_d6);
        const double F_LJ_scalar = 24.0 * epsilon_ij * inv_distSq * (sigma6_d6 - 2.0 * sigma12_d12);

        double fx = 0.0;
        double fy = 0.0;
        double fz = 0.0;

        if (distSq <= smoothingSq) {
          fx = F_LJ_scalar * dx;
          fy = F_LJ_scalar * dy;
          fz = F_LJ_scalar * dz;
        } else {
          const double d = std::sqrt(distSq);
          const double dMinusRl = d - smoothingR;
          const double dMinusRl2 = dMinusRl * dMinusRl;
          const double S = 1.0 - dMinusRl2 * (3.0 * cutoffR - smoothingR - 2.0 * d) / rcRlCubed;
          const double dS_dd = -6.0 * dMinusRl * (cutoffR - d) / rcRlCubed;

          const double F_total_scalar = S * F_LJ_scalar + U_LJ * dS_dd / d;
          fx = F_total_scalar * dx;
          fy = F_total_scalar * dy;
          fz = F_total_scalar * dz;
        }

        auto& fi = p_i.getF();
        auto& fj = p_j.getF();
        fi[0] += fx;
        fi[1] += fy;
        fi[2] += fz;
        fj[0] -= fx;
        fj[1] -= fy;
        fj[2] -= fz;
      });

  // Apply boundary conditions
  applyPeriodicBoundaries(lc);
}

void SmoothedLJForce::precomputeConstants() {
  // Determine the highest type number for a particle
  int maxTypeNr = 0;
  for (const auto& p : particles) {
    if (p.getType() > maxTypeNr) {
      maxTypeNr = p.getType();
    }
  }
  tableWidth = maxTypeNr + 1;

  const int tableSize = tableWidth * tableWidth;
  pairEpsilon.clear();
  pairSigma6.clear();
  pairSigma12.clear();
  pairEpsilon.resize(tableSize, -1.0);
  pairSigma6.resize(tableSize, -1.0);
  pairSigma12.resize(tableSize, -1.0);

  // Collect unique types and their sigma/epsilon
  std::vector<int> uniqueTypes;
  std::vector<double> typeSigma(tableWidth, -1.0);
  std::vector<double> typeEpsilon(tableWidth, -1.0);

  for (const auto& p : particles) {
    const int t = p.getType();
    if (typeSigma[t] < 0) {
      uniqueTypes.push_back(t);
      typeSigma[t] = p.getSigma();
      typeEpsilon[t] = p.getEpsilon();
    }
  }

  // Compute lookup tables for all type pairs
  for (const int ti : uniqueTypes) {
    for (const int tj : uniqueTypes) {
      const int idx = ti * tableWidth + tj;
      if (pairEpsilon[idx] < 0) {
        // Lorentz-Berthelot mixing rules
        const double sigma_ij = (typeSigma[ti] + typeSigma[tj]) * 0.5;
        const double epsilon_ij = std::sqrt(typeEpsilon[ti] * typeEpsilon[tj]);

        const double sigma2 = sigma_ij * sigma_ij;
        const double sigma6 = sigma2 * sigma2 * sigma2;
        const double sigma12 = sigma6 * sigma6;

        pairEpsilon[idx] = epsilon_ij;
        pairSigma6[idx] = sigma6;
        pairSigma12[idx] = sigma12;

        // Symmetric
        const int idxSym = tj * tableWidth + ti;
        pairEpsilon[idxSym] = epsilon_ij;
        pairSigma6[idxSym] = sigma6;
        pairSigma12[idxSym] = sigma12;
      }
    }
  }

  SPDLOG_DEBUG("SmoothedLJForce: Precomputed constants for {} unique particle types", uniqueTypes.size());
}

void SmoothedLJForce::calcFPeriodicPair(Particle* p1, Particle* p2) const {
  applySmoothedPairForce(*p1, *p2);
}

void SmoothedLJForce::applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const {
  applyPeriodicBoundariesImpl(lc, [this](Particle* p1, Particle* p2) { calcFPeriodicPair(p1, p2); });
}

// applySmoothedPairForce is now defined inline in SmoothedLJForce.h for proper inlining
