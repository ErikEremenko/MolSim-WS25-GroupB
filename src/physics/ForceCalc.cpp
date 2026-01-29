#include "physics/ForceCalc.h"

#include <omp.h>
#include <spdlog/spdlog.h>
#include <cmath>

#include "physics/LinkedCellParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "utils/PhysicsConstants.h"

using namespace PhysicsConstants;

ForceCalc::ForceCalc(ParticleContainer& particles) : particles(particles) {}
ForceCalc::~ForceCalc() = default;

void ForceCalc::calculateX(ParticleContainer& particles, const double dt) {
  const double dt_sq_half = 0.5 * dt * dt;
  const size_t n = particles.size();

  // Manual loop (better SIMD optimization potential)
  for (size_t i = 0; i < n; ++i) {
    auto& p = particles[i];
    const auto& x_curr = p.getX();
    const auto& v = p.getV();
    const auto& F = p.getF();
    const double inv_m = 1.0 / p.getM();

    // compute all 3 components (SIMD friendly)
    std::array<double, 3> x_new{};
#pragma omp simd
    for (int d = 0; d < 3; ++d) {
      x_new[d] = x_curr[d] + dt * v[d] + dt_sq_half * inv_m * F[d];
    }
    p.setX(x_new);
  }
}

void ForceCalc::calculateV(ParticleContainer& particles, const double dt) {
  const size_t n = particles.size();

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

void GravityForce::calculateF() {
  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // Index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      const double norm = ArrayUtils::L2Norm(dist);

      const double norm3 = norm * norm * norm;

      const auto F_vector = ((p_i.getM() * p_j.getM()) / norm3) * dist;

      // Apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
      auto F_i = p_i.getF();
      auto F_j = p_j.getF();
      // Actio est reactio
      p_i.setF(F_i + F_vector);
      p_j.setF(F_j - F_vector);
    }
  }
}

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

  for (size_t i = 0; i < n; ++i) {
    Particle& p = particles[i];
    auto& F = p.getF();
    F[ax] += p.getM() * g;
  }
}

LennardJonesForce::LennardJonesForce(ParticleContainer& particles, const double epsilon, const double sigma,
                                     const double cutoffRadius)
    : ForceCalc(particles),
      epsilon(epsilon),
      sigma(sigma),
      cutoffRadius(cutoffRadius),
      repulsionDistance(TWO_POW_1_6 * sigma),
      cutoffRadiusSq(cutoffRadius * cutoffRadius) {}

void LennardJonesForce::calculateF() {
  if (dynamic_cast<LinkedCellParticleContainer*>(&particles)) {
    calculateFLinkedCell();
  } else {
    calculateFDirectSum();
  }
}

void LennardJonesForce::calculateFDirectSum() {
  for (auto& p : particles) {
    p.setF({});
  }
  const double sigma2 = sigma * sigma;
  const double sigma6 = sigma2 * sigma2 * sigma2;
  const double cutoffRadiusSqLocal = cutoffRadius * cutoffRadius;

  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // Index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      applyLJPairForceGlobal(particles[i], particles[j], sigma6, epsilon, cutoffRadiusSqLocal);
    }
  }
}

LennardJonesForceParallel::LennardJonesForceParallel(ParticleContainer& particles, const double epsilon,
                                                     const double sigma, const double cutoffRadius)
    : ForceCalc(particles), epsilon(epsilon), sigma(sigma), cutoffRadius(cutoffRadius) {}

void LennardJonesForceParallel::calculateF() {
#pragma omp parallel for
  for (auto& p : particles) {
    p.setF({});
  }
  const double sigma2 = sigma * sigma;
  const double sigma6 = sigma2 * sigma2 * sigma2;
  const double cutoffRadiusSqLocal = cutoffRadius * cutoffRadius;
  const size_t n_particles = particles.size();

#pragma omp parallel for schedule(guided)
  for (size_t i = 0; i < n_particles; ++i) {
    // Index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      // Use squared distance to avoid sqrt (optimization)
      const double normSq = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

      // Skip if beyond cutoff (no zero check for performance - see class documentation)
      if (normSq >= cutoffRadiusSqLocal) {
        continue;
      }

      const double inv_norm2 = 1.0 / normSq;
      const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

      const double crossing_norm_quot_6 = sigma6 * inv_norm6;
      const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

      const auto F_vector = (24.0 * epsilon) * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12) * dist;

      auto& Fi = p_i.getF();
      auto& Fj = p_j.getF();
      // Apply forces using Newton's third law and atomic operations
#pragma omp atomic
      Fi[0] += F_vector[0];
#pragma omp atomic
      Fi[1] += F_vector[1];
#pragma omp atomic
      Fi[2] += F_vector[2];

#pragma omp atomic
      Fj[0] -= F_vector[0];
#pragma omp atomic
      Fj[1] -= F_vector[1];
#pragma omp atomic
      Fj[2] -= F_vector[2];
    }
  }
}

void LennardJonesForce::calculateFLinkedCell() {
  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);
  if (!lc) {
    throw std::runtime_error("LennardJonesForce::calculateFLinkedCell requires LinkedCellParticleContainer");
  }
  lc->handleOutflowBoundaries();

  // Use template version to eliminate std::function overhead (~85% overhead on VTune)
  // Cache lookup table pointers for better optimization
  const double* __restrict__ lut1 = pairLookupTable1.data();
  const double* __restrict__ lut2 = pairLookupTable2.data();
  const int tw = tableWidth;
  const double cutoffSq = cutoffRadiusSq;

  // CRITICAL: Force calculation MUST be directly in the lambda for proper inlining.
  // Calling a separate function (even with always_inline) prevents compiler optimizations.
  lc->iteratePairsTemplate([lut1, lut2, tw, cutoffSq](Particle& p_i, Particle& p_j) {
    // Get positions (getX() force-inlined)
    const auto& xi = p_i.getX();
    const auto& xj = p_j.getX();

    // Compute distance vector and squared distance
    const double dx = xj[0] - xi[0];
    const double dy = xj[1] - xi[1];
    const double dz = xj[2] - xi[2];
    const double distSq = dx * dx + dy * dy + dz * dz;

    if (distSq > cutoffSq)
      return;

    // Compute force using lookup tables
    const double inv_distSq = 1.0 / distSq;
    const int idx = p_i.getType() * tw + p_j.getType();
    const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
    const double term = lut1[idx] * inv_distSq * inv_distSq3 * (lut2[idx] - inv_distSq3);

    const double fx = term * dx;
    const double fy = term * dy;
    const double fz = term * dz;

    // Update forces (getF() force-inlined)
    auto& fi = p_i.getF();
    auto& fj = p_j.getF();
    fi[0] += fx;
    fi[1] += fy;
    fi[2] += fz;
    fj[0] -= fx;
    fj[1] -= fy;
    fj[2] -= fz;
  });

  applyReflectiveBoundaries(lc);
  applyPeriodicBoundaries(lc);
}

void LennardJonesForce::applyReflectiveBoundaries(const LinkedCellParticleContainer* lc) const {
  // Ghost particle mirrored across wall at dimension d has ghostDist with only
  // one non-zero component: ghostDist[d] = 2*(wallPos - x[d])
  //   normSq = ghostDist[d]^2 = 4 * (wallPos - x[d])^2 = 4 * distToWall^2
  // -> simpler norm calculation

  const auto& domainOrigin = lc->domain_origin();
  const auto& domainDims = lc->domain_dims();
  const auto& boundaryTypes = lc->boundary_types();

  // Precompute boundary position
  const double minX = domainOrigin[0], maxX = domainOrigin[0] + domainDims[0];
  const double minY = domainOrigin[1], maxY = domainOrigin[1] + domainDims[1];
  const double minZ = domainOrigin[2], maxZ = domainOrigin[2] + domainDims[2];

  // Cache which boundaries are reflective
  const bool refXMin = boundaryTypes[0] == BoundaryType::REFLECTIVE;
  const bool refXMax = boundaryTypes[1] == BoundaryType::REFLECTIVE;
  const bool refYMin = boundaryTypes[2] == BoundaryType::REFLECTIVE;
  const bool refYMax = boundaryTypes[3] == BoundaryType::REFLECTIVE;
  const bool refZMin = boundaryTypes[4] == BoundaryType::REFLECTIVE;
  const bool refZMax = boundaryTypes[5] == BoundaryType::REFLECTIVE;

  for (auto& p : particles) {
    const auto& x = p.getX();
    auto F_total = p.getF();

    // Use precomputed per-type lookup tables
    const int pType = p.getType();
    const double repDistSq = repulsionDistanceSqLookup[pType];  // (2^(1/6) * sigma)^2
    const double sigma6 = sigma6Lookup[pType];
    const double epsilon24 = epsilon24Lookup[pType];  // 24 * epsilon

    // Lambda for computing 1D ghost force contribution (inlined by compiler)
    // normSq = 4 * distToWall^2, ghostDistD = ±2 * distToWall
    auto applyGhostForce1D = [&](const int d, const double distToWall, const double sign) {
      const double normSq = 4.0 * distToWall * distToWall;
      if (normSq < repDistSq && normSq > 0.0) {
        const double inv_norm2 = 1.0 / normSq;
        const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;
        const double s6_inv6 = sigma6 * inv_norm6;
        const double s12_inv12 = s6_inv6 * s6_inv6;
        // ghostDist[d] = sign * 2.0 * distToWall (negative for min wall, pos for max wall)
        const double ghostDistD = sign * 2.0 * distToWall;
        F_total[d] += epsilon24 * inv_norm2 * (s6_inv6 - 2.0 * s12_inv12) * ghostDistD;
      }
    };

    if (refXMin) {
      if (const double distToWall = x[0] - minX; distToWall > 0.0 && distToWall < cutoffRadius)
        applyGhostForce1D(0, distToWall, -1.0);
    }
    if (refXMax) {
      if (const double distToWall = maxX - x[0]; distToWall > 0.0 && distToWall < cutoffRadius)
        applyGhostForce1D(0, distToWall, 1.0);
    }

    if (refYMin) {
      if (const double distToWall = x[1] - minY; distToWall > 0.0 && distToWall < cutoffRadius)
        applyGhostForce1D(1, distToWall, -1.0);
    }
    if (refYMax) {
      if (const double distToWall = maxY - x[1]; distToWall > 0.0 && distToWall < cutoffRadius)
        applyGhostForce1D(1, distToWall, 1.0);
    }

    if (refZMin) {
      if (const double distToWall = x[2] - minZ; distToWall > 0.0 && distToWall < cutoffRadius)
        applyGhostForce1D(2, distToWall, -1.0);
    }
    if (refZMax) {
      if (const double distToWall = maxZ - x[2]; distToWall > 0.0 && distToWall < cutoffRadius)
        applyGhostForce1D(2, distToWall, 1.0);
    }

    p.setF(F_total);
  }
}

void LennardJonesForce::calcFPeriodicBoundary(Particle* p1, Particle* p2) const {
  applyLJPairForceLookup(*p1, *p2, pairLookupTable1.data(), pairLookupTable2.data(), tableWidth, cutoffRadiusSq);
}

void LennardJonesForce::applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const {
  applyPeriodicBoundariesImpl(lc, [this](Particle* p1, Particle* p2) { calcFPeriodicBoundary(p1, p2); });
}

void LennardJonesForce::applyLJPairForceGlobal(Particle& p_i, Particle& p_j, const double sigma6, const double epsilon,
                                               const double cutoffSq) {
  const auto dist = p_j.getX() - p_i.getX();
  const double normSq = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

  if (normSq >= cutoffSq) {
    return;
  }

  const double inv_norm2 = 1.0 / normSq;
  const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

  const double crossing_norm_quot_6 = sigma6 * inv_norm6;
  const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

  const auto F_vector = (24.0 * epsilon) * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12) * dist;

  auto F_i = p_i.getF();
  auto F_j = p_j.getF();
  p_i.setF(F_i + F_vector);
  p_j.setF(F_j - F_vector);
}

// applyLJPairForceLookup is now defined inline in ForceCalc.h for proper inlining

void LennardJonesForce::precomputeConstants() {
  // Determine the highest type number for a particle
  int maxTypeNr = 0;
  for (const auto& p : particles) {
    if (p.getType() > maxTypeNr)
      maxTypeNr = p.getType();
  }
  tableWidth = maxTypeNr + 1;

  // Pre-allocate lookup tables with proper size (avoids repeated push_back)
  const int tableSize = tableWidth * tableWidth;
  pairLookupTable1.clear();
  pairLookupTable2.clear();
  pairLookupTable1.resize(tableSize, -1.0);
  pairLookupTable2.resize(tableSize, -1.0);

  // Pre-allocate lookup tables for reflective boundaries
  repulsionDistanceSqLookup.clear();
  sigma6Lookup.clear();
  epsilon24Lookup.clear();
  repulsionDistanceSqLookup.resize(tableWidth, -1.0);
  sigma6Lookup.resize(tableWidth, -1.0);
  epsilon24Lookup.resize(tableWidth, -1.0);

  // Build a set of unique (type_i, type_j) pairs to avoid O(n^2) iteration
  // Collect unique types first, then iterate over type pairs
  std::vector<int> uniqueTypes;
  std::vector<double> typeSigma(tableWidth, -1.0);
  std::vector<double> typeEpsilon(tableWidth, -1.0);

  for (const auto& p : particles) {
    const int t = p.getType();
    if (typeSigma[t] < 0) {
      // First particle of this type -> record its sigma/epsilon
      uniqueTypes.push_back(t);
      typeSigma[t] = p.getSigma();
      typeEpsilon[t] = p.getEpsilon();

      // Precompute per-type constants for reflective boundaries
      const double pSigma = p.getSigma();
      const double pEpsilon = p.getEpsilon();
      const double sigma2 = pSigma * pSigma;
      const double sigma6 = sigma2 * sigma2 * sigma2;

      repulsionDistanceSqLookup[t] = TWO_POW_1_3 * sigma2;
      sigma6Lookup[t] = sigma6;
      epsilon24Lookup[t] = 24.0 * pEpsilon;
    }
  }

  // For every pair of unique types, compute and store lookup values
  for (const int ti : uniqueTypes) {
    for (const int tj : uniqueTypes) {
      const int idx = ti * tableWidth + tj;
      if (pairLookupTable1[idx] < 0) {
        // Apply Lorentz-Berthelot mixing rules
        const double sigma_ij = (typeSigma[ti] + typeSigma[tj]) * 0.5;
        const double epsilon_ij = std::sqrt(typeEpsilon[ti] * typeEpsilon[tj]);

        // Precompute: 48 * epsilon_ij * sigma_ij^12
        const double sigma2 = sigma_ij * sigma_ij;
        const double sigma4 = sigma2 * sigma2;
        const double sigma6 = sigma4 * sigma2;
        const double sigma12 = sigma6 * sigma6;
        const double eps_48_sigma_12 = 48.0 * epsilon_ij * sigma12;

        // Precompute: 1 / (2 * sigma_ij^6)
        const double sigma_m6_2 = 0.5 / sigma6;

        pairLookupTable1[idx] = eps_48_sigma_12;
        pairLookupTable2[idx] = sigma_m6_2;

        // Symmetric: save to (tj, ti) as well
        const int idxSym = tj * tableWidth + ti;
        pairLookupTable1[idxSym] = eps_48_sigma_12;
        pairLookupTable2[idxSym] = sigma_m6_2;
      }
    }
  }

  SPDLOG_DEBUG("Precomputed constants for {} unique particle types ({} pair combinations)", uniqueTypes.size(),
               uniqueTypes.size() * uniqueTypes.size());
}

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

  // CRITICAL: Force calculation MUST be directly in the lambda for proper inlining.
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

// applySmoothedPairForce is now defined inline in ForceCalc.h for proper inlining

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

ConstantForce::ConstantForce(ParticleContainer& particles, double fx, double fy, double fz, const double endTime,
                             double& currentTime, std::vector<std::pair<int, int>> targetIndices,
                             const int membraneDimY)
    : ForceCalc(particles),
      force({fx, fy, fz}),
      endTime(endTime),
      currentTime(currentTime),
      targetIndices(std::move(targetIndices)),
      membraneDimY(membraneDimY) {}

void ConstantForce::calculateF() {
  // Only apply force before endTime
  if (currentTime >= endTime) {
    return;
  }

  // Apply constant force to target particles identified by their membrane grid indices
  for (const auto& [gridX, gridY] : targetIndices) {
    // Calculate particle index from grid coordinates (matching generateMembrane layout)
    // The particle at (gridX, gridY) has index gridX * yDim + gridY

    if (const int particleIdx = gridX * membraneDimY + gridY;
        particleIdx >= 0 && particleIdx < static_cast<int>(particles.size())) {
      Particle& p = particles[particleIdx];
      p.setF(p.getF() + force);
    }
  }
}