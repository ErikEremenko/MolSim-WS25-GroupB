#include "physics/LennardJonesForce.h"

#include <omp.h>
#include <spdlog/spdlog.h>
#include <cmath>

#include "physics/LinkedCellParticleContainer.h"
#include "utils/ArrayUtils.h"
#include "utils/PhysicsConstants.h"

using namespace PhysicsConstants;

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
#ifdef _OPENMP
    if (useParallel) {
      if (parallelStrategy == ParallelStrategy::COLORING) {
        calculateFLinkedCellParallel1();
      } else {
        calculateFLinkedCellParallel2();
      }
    } else {
      calculateFLinkedCell();
    }
#else
    calculateFLinkedCell();
#endif
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
    // Load positions into local variables (helps register allocation)
    const double xi0 = p_i.getX()[0];
    const double xi1 = p_i.getX()[1];
    const double xi2 = p_i.getX()[2];
    const double xj0 = p_j.getX()[0];
    const double xj1 = p_j.getX()[1];
    const double xj2 = p_j.getX()[2];

    // Compute distance vector and squared distance
    const double dx = xj0 - xi0;
    const double dy = xj1 - xi1;
    const double dz = xj2 - xi2;
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
  // For reflective walls, the ghost distance is only along one axis.
  // That lets us compute normSq as 4 * distToWall * distToWall.

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

  // Parallelize (each particle's reflective boundary force is independent)
  const size_t n = particles.size();
#pragma omp parallel for schedule(static)
  for (size_t i = 0; i < n; ++i) {
    auto& p = particles[i];
    const auto& x = p.getX();
    auto F_total = p.getF();

    // Use precomputed per-type lookup tables
    const int pType = p.getType();
    const double repDistSq = repulsionDistanceSqLookup[pType];
    const double sigma6 = sigma6Lookup[pType];
    const double epsilon24 = epsilon24Lookup[pType];

    // 1D ghost force helper
    auto applyGhostForce1D = [&](const int d, const double distToWall, const double sign) {
      const double normSq = 4.0 * distToWall * distToWall;
      if (normSq < repDistSq && normSq > 0.0) {
        const double inv_norm2 = 1.0 / normSq;
        const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;
        const double s6_inv6 = sigma6 * inv_norm6;
        const double s12_inv12 = s6_inv6 * s6_inv6;
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

void LennardJonesForce::calcFPeriodicBoundaryAtomic(Particle* p1, const std::array<double, 3>& p2_shifted_pos,
                                                    Particle* p2) const {
  // Version of calcFPeriodicBoundary that uses atomics for thread-safe force updates
  // Takes the pre-shifted position to avoid modifying the position of p2
  const std::array<double, 3> dist = {p2_shifted_pos[0] - p1->getX()[0], p2_shifted_pos[1] - p1->getX()[1],
                                      p2_shifted_pos[2] - p1->getX()[2]};

  double term = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

  if (term >= cutoffRadiusSq)
    return;

  term = 1.0 / term;

  const int idx = p1->getType() * tableWidth + p2->getType();
  term = pairLookupTable1[idx] * term * term * term * term * (pairLookupTable2[idx] - term * term * term);

  const double fx = term * dist[0];
  const double fy = term * dist[1];
  const double fz = term * dist[2];

  auto& f1 = p1->getF();
  auto& f2 = p2->getF();
#pragma omp atomic
  f1[0] += fx;
#pragma omp atomic
  f1[1] += fy;
#pragma omp atomic
  f1[2] += fz;
#pragma omp atomic
  f2[0] -= fx;
#pragma omp atomic
  f2[1] -= fy;
#pragma omp atomic
  f2[2] -= fz;
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

// applyLJPairForceLookup is now defined inline in LennardJonesForce.h for proper inlining

void LennardJonesForce::applyPeriodicBoundariesParallel(LinkedCellParticleContainer* lc) const {
  // Parallelized version of applyPeriodicBoundaries using atomics for thread-safe force updates
  const auto domainDims = lc->domain_dims();
  const auto boundaryTypes = lc->boundary_types();
  const auto numCells = lc->num_cells();

  // X-periodic: left wall (x=1) <-> right wall (x=numCells[0]-2)
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
        auto& cell1 = lc->cell_at(1, c1y, c1z);
        for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
          for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
            if (c2y < 1 || c2y > numCells[1] - 2 || c2z < 1 || c2z > numCells[2] - 2)
              continue;

            auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, c2z);
            for (auto& p1 : cell1) {
              for (auto& p2 : cell2) {
                std::array<double, 3> p2_shifted = p2->getX();
                p2_shifted[0] -= domainDims[0];
                calcFPeriodicBoundaryAtomic(p1, p2_shifted, p2);
              }
            }
          }
        }
      }
    }
  }

  // Y-periodic: bottom wall (y=1) <-> top wall (y=numCells[1]-2)
  if (boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
        auto& cell1 = lc->cell_at(c1x, 1, c1z);
        for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {
          for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
            if (c2x < 1 || c2x > numCells[0] - 2 || c2z < 1 || c2z > numCells[2] - 2)
              continue;

            auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, c2z);
            for (auto& p1 : cell1) {
              for (auto& p2 : cell2) {
                std::array<double, 3> p2_shifted = p2->getX();
                p2_shifted[1] -= domainDims[1];
                calcFPeriodicBoundaryAtomic(p1, p2_shifted, p2);
              }
            }
          }
        }
      }
    }
  }

  // Z-periodic: front wall (z=1) <-> back wall (z=numCells[2]-2)
  if (boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
#pragma omp parallel for collapse(2) schedule(static)
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
        auto& cell1 = lc->cell_at(c1x, c1y, 1);
        for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {
          for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
            if (c2x < 1 || c2x > numCells[0] - 2 || c2y < 1 || c2y > numCells[1] - 2)
              continue;

            auto& cell2 = lc->cell_at(c2x, c2y, numCells[2] - 2);
            for (auto& p1 : cell1) {
              for (auto& p2 : cell2) {
                std::array<double, 3> p2_shifted = p2->getX();
                p2_shifted[2] -= domainDims[2];
                calcFPeriodicBoundaryAtomic(p1, p2_shifted, p2);
              }
            }
          }
        }
      }
    }
  }

  // Edge and corner interactions remain sequential (small number of cells)
  // X+Y periodic edges
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC) {
    for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
      auto& cell1 = lc->cell_at(1, 1, c1z);
      for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
        if (c2z < 1 || c2z > numCells[0] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, c2z);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[1] - domainDims[1], 1);
            calcFPeriodicBoundary(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
          }
      }
    }
    for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
      auto& cell1 = lc->cell_at(1, numCells[1] - 2, c1z);
      for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
        if (c2z < 1 || c2z > numCells[0] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, 1, c2z);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
            calcFPeriodicBoundary(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[1] - domainDims[1], 1);
          }
      }
    }
  }

  // X+Z periodic edges
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      auto& cell1 = lc->cell_at(1, c1y, 1);
      for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
        if (c2y < 1 || c2y > numCells[1] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, numCells[2] - 2);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
            calcFPeriodicBoundary(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
          }
      }
    }
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      auto& cell1 = lc->cell_at(1, c1y, numCells[2] - 2);
      for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
        if (c2y < 1 || c2y > numCells[1] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, 1);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
            calcFPeriodicBoundary(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
          }
      }
    }
  }

  // Y+Z periodic edges
  if (boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      auto& cell1 = lc->cell_at(c1x, 1, 1);
      for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {
        if (c2x < 1 || c2x > numCells[0] - 2)
          continue;
        auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, numCells[2] - 2);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[1] - domainDims[1], 1);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
            calcFPeriodicBoundary(p1, p2);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
          }
      }
    }
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      auto& cell1 = lc->cell_at(c1x, 1, numCells[2] - 2);
      for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {
        if (c2x < 1 || c2x > numCells[0] - 2)
          continue;
        auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, 1);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[1] - domainDims[1], 1);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
            calcFPeriodicBoundary(p1, p2);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
          }
      }
    }
  }

  // Corner interactions (fully periodic 3D)
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {

    auto& cell1 = lc->cell_at(1, 1, 1);
    auto& cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, numCells[2] - 2);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
        calcFPeriodicBoundary(p1, p2);
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
      }

    cell1 = lc->cell_at(numCells[0] - 2, 1, 1);
    cell2 = lc->cell_at(1, numCells[1] - 2, numCells[2] - 2);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
        calcFPeriodicBoundary(p1, p2);
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
      }

    cell1 = lc->cell_at(1, 1, numCells[2] - 2);
    cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, 1);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
        calcFPeriodicBoundary(p1, p2);
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
      }

    cell1 = lc->cell_at(numCells[0] - 2, 1, numCells[2] - 2);
    cell2 = lc->cell_at(1, numCells[1] - 2, 1);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
        calcFPeriodicBoundary(p1, p2);
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
      }
  }
}

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

void LennardJonesForce::calculateFLinkedCellParallel1() {
  // Strategy 1: Domain decomposition with static scheduling (coloring approach)
  // Cells in the same phase are spaced so their write regions do not overlap.
  // This avoids races without atomics.
  // 2D uses 3x2 phases, 3D uses 3x2x3 phases.

  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);
  if (!lc) {
    throw std::runtime_error("LennardJonesForce::calculateFLinkedCell requires LinkedCellParticleContainer");
  }
  lc->handleOutflowBoundaries();

  const auto numCells = lc->num_cells();
  const int nx = numCells[0];
  const int ny = numCells[1];
  const int nz = numCells[2];

  const bool is2D = (nz <= 3);

  if (is2D) {
    // 2D path: 3x2 coloring (6 phases)
    for (int px = 0; px < 3; px++)
      for (int py = 0; py < 2; py++) {
#pragma omp parallel for collapse(2) schedule(static)
        for (int lx = 1 + px; lx < nx - 1; lx += 3)
          for (int ly = 1 + py; ly < ny - 1; ly += 2) {
            for (int lz = 1; lz < nz - 1; lz++) {

              auto& cell1 = lc->cell_at(lx, ly, lz);

              // Intra-cell pairs
              for (size_t i = 0; i < cell1.size(); ++i)
                for (size_t j = i + 1; j < cell1.size(); ++j) {
                  calcFPeriodicBoundary(cell1[i], cell1[j]);
                }

              // Inter-cell pairs (13-neighbor half-shell)
              for (int dx = -1; dx < 2; dx++)
                for (int dz = -1; dz < 2; dz++)
                  for (int dy = 0; dy < 2; dy++) {
                    if (dy == 0 && (dz == -1 || (dz == 0 && dx <= 0)))
                      continue;
                    int c2x = lx + dx, c2y = ly + dy, c2z = lz + dz;
                    if (c2x < 1 || c2x >= nx - 1 || c2y < 1 || c2y >= ny - 1 || c2z < 1 || c2z >= nz - 1)
                      continue;

                    auto& cell2 = lc->cell_at(c2x, c2y, c2z);
                    for (auto& p1 : cell1)
                      for (auto& p2 : cell2) {
                        calcFPeriodicBoundary(p1, p2);
                      }
                  }
            }
          }
      }
  } else {
    // 3D path: 3x2x3 coloring (18 phases)
    for (int px = 0; px < 3; px++)
      for (int py = 0; py < 2; py++)
        for (int pz = 0; pz < 3; pz++) {
#pragma omp parallel for collapse(3) schedule(static)
          for (int lx = 1 + px; lx < nx - 1; lx += 3)
            for (int ly = 1 + py; ly < ny - 1; ly += 2)
              for (int lz = 1 + pz; lz < nz - 1; lz += 3) {

                auto& cell1 = lc->cell_at(lx, ly, lz);

                // Intra-cell pairs
                for (size_t i = 0; i < cell1.size(); ++i)
                  for (size_t j = i + 1; j < cell1.size(); ++j) {
                    calcFPeriodicBoundary(cell1[i], cell1[j]);
                  }

                // Inter-cell pairs (13-neighbor half-shell)
                for (int dx = -1; dx < 2; dx++)
                  for (int dz = -1; dz < 2; dz++)
                    for (int dy = 0; dy < 2; dy++) {
                      if (dy == 0 && (dz == -1 || (dz == 0 && dx <= 0)))
                        continue;
                      int c2x = lx + dx, c2y = ly + dy, c2z = lz + dz;
                      if (c2x < 1 || c2x >= nx - 1 || c2y < 1 || c2y >= ny - 1 || c2z < 1 || c2z >= nz - 1)
                        continue;

                      auto& cell2 = lc->cell_at(c2x, c2y, c2z);
                      for (auto& p1 : cell1)
                        for (auto& p2 : cell2) {
                          calcFPeriodicBoundary(p1, p2);
                        }
                    }
              }
        }
  }

  applyReflectiveBoundaries(lc);
  applyPeriodicBoundariesParallel(lc);
}

void LennardJonesForce::calculateFLinkedCellParallel2() {
  // Strategy 2: Task-based parallelization with atomics
  // - Creates one task per cell for dynamic load balancing
  // - Uses atomic operations for force updates (race condition safe)
  // - Work-stealing provides better balance for inhomogeneous particle distributions

  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);
  if (!lc) {
    throw std::runtime_error("LennardJonesForce::calculateFLinkedCellParallel2 requires LinkedCellParticleContainer");
  }
  lc->handleOutflowBoundaries();

  const auto numCells = lc->num_cells();
  const int nx = numCells[0];
  const int ny = numCells[1];
  const int nz = numCells[2];

  // Cache lookup table data for lambda capture
  const double* __restrict__ lut1 = pairLookupTable1.data();
  const double* __restrict__ lut2 = pairLookupTable2.data();
  const int tw = tableWidth;
  const double cutoffSq = cutoffRadiusSq;

#pragma omp parallel
  {
#pragma omp single
    {
      // Spawning one task per cell (OpenMP runtime handles work distribution)
      for (int lx = 1; lx < nx - 1; lx++) {
        for (int ly = 1; ly < ny - 1; ly++) {
          for (int lz = 1; lz < nz - 1; lz++) {
#pragma omp task firstprivate(lx, ly, lz)
            {
              auto& cell1 = lc->cell_at(lx, ly, lz);

              // Intra-cell pairs
              for (size_t i = 0; i < cell1.size(); ++i) {
                for (size_t j = i + 1; j < cell1.size(); ++j) {
                  Particle* p1 = cell1[i];
                  Particle* p2 = cell1[j];

                  const auto& x1 = p1->getX();
                  const auto& x2 = p2->getX();
                  const double dx = x2[0] - x1[0];
                  const double dy = x2[1] - x1[1];
                  const double dz = x2[2] - x1[2];
                  const double distSq = dx * dx + dy * dy + dz * dz;

                  if (distSq >= cutoffSq)
                    continue;

                  const double inv_distSq = 1.0 / distSq;
                  const int idx = p1->getType() * tw + p2->getType();
                  const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
                  const double term = lut1[idx] * inv_distSq * inv_distSq3 * (lut2[idx] - inv_distSq3);

                  const double fx = term * dx;
                  const double fy = term * dy;
                  const double fz = term * dz;

                  // Intra-cell: must still use atomics because a neighboring cell's
                  // inter-cell task may concurrently update the same particles
                  auto& f1 = p1->getF();
                  auto& f2 = p2->getF();
#pragma omp atomic
                  f1[0] += fx;
#pragma omp atomic
                  f1[1] += fy;
#pragma omp atomic
                  f1[2] += fz;
#pragma omp atomic
                  f2[0] -= fx;
#pragma omp atomic
                  f2[1] -= fy;
#pragma omp atomic
                  f2[2] -= fz;
                }
              }

              // Inter-cell pairs
              // Only process forward neighbors to avoid double-counting
              for (int ddx = -1; ddx <= 1; ddx++) {
                for (int ddz = -1; ddz <= 1; ddz++) {
                  for (int ddy = 0; ddy <= 1; ddy++) {
                    // Skipping self and backward neighbors
                    if (ddy == 0 && (ddz < 0 || (ddz == 0 && ddx <= 0)))
                      continue;

                    int c2x = lx + ddx, c2y = ly + ddy, c2z = lz + ddz;
                    if (c2x < 1 || c2x >= nx - 1 || c2y < 1 || c2y >= ny - 1 || c2z < 1 || c2z >= nz - 1)
                      continue;

                    auto& cell2 = lc->cell_at(c2x, c2y, c2z);

                    for (auto& p1 : cell1) {
                      for (auto& p2 : cell2) {
                        const auto& x1 = p1->getX();
                        const auto& x2 = p2->getX();
                        const double dx_p = x2[0] - x1[0];
                        const double dy_p = x2[1] - x1[1];
                        const double dz_p = x2[2] - x1[2];
                        const double distSq = dx_p * dx_p + dy_p * dy_p + dz_p * dz_p;

                        if (distSq >= cutoffSq)
                          continue;

                        const double inv_distSq = 1.0 / distSq;
                        const int idx = p1->getType() * tw + p2->getType();
                        const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
                        const double term = lut1[idx] * inv_distSq * inv_distSq3 * (lut2[idx] - inv_distSq3);

                        const double fx = term * dx_p;
                        const double fy = term * dy_p;
                        const double fz = term * dz_p;

                        // Inter-cell: using atomics for thread safety
                        auto& f1 = p1->getF();
                        auto& f2 = p2->getF();
#pragma omp atomic
                        f1[0] += fx;
#pragma omp atomic
                        f1[1] += fy;
#pragma omp atomic
                        f1[2] += fz;
#pragma omp atomic
                        f2[0] -= fx;
#pragma omp atomic
                        f2[1] -= fy;
#pragma omp atomic
                        f2[2] -= fz;
                      }
                    }
                  }
                }
              }
            }  // end task
          }
        }
      }
    }  // end single (implicit taskwait at end of single)
  }  // end parallel

  applyReflectiveBoundaries(lc);
  applyPeriodicBoundariesParallel(lc);
}
