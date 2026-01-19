#include "physics/ForceCalc.h"

#include <math.h>  // TODO: Replace with <cmath>, never use C headers in C++!
#include <spdlog/spdlog.h>

#include "physics/LinkedCellParticleContainer.h"
#include "simulation/SimulationConfig.h"  // For boundary types, TODO: Refactor boundary types into another file
#include "utils/ArrayUtils.h"

ForceCalc::~ForceCalc() = default;

void ForceCalc::calculateX(const double dt) {
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

void ForceCalc::calculateV(const double dt) {
  for (auto& p : particles) {
    const auto v_curr = p.getV();
    const double m_i = p.getM();
    const auto F = p.getF();
    const auto F_old = p.getOldF();
    const auto v_new = v_curr + (dt / (2.0 * m_i)) * (F_old + F);
    p.setV(v_new);
  }
}

void ForceCalc::precomputeConstants() {}

void GravityForce::calculateF() {
  for (auto& p : particles) {
    p.setF({});
  }

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

LennardJonesForce::LennardJonesForce(ParticleContainer& particles, const double epsilon, const double sigma,
                                     const double cutoffRadius, const double gravity)
    : ForceCalc(particles),
      epsilon(epsilon),
      sigma(sigma),
      cutoffRadius(cutoffRadius),
      repulsionDistance(std::pow(2.0, 1.0 / 6.0) * sigma),
      gravity(gravity),
      cutoffRadiusSq (cutoffRadius * cutoffRadius) {}

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

      // Apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
      auto F_i = p_i.getF();
      auto F_j = p_j.getF();
      // Actio est reactio
      p_i.setF(F_i + F_vector);
      p_j.setF(F_j - F_vector);
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
  for (auto& p : particles) {
    p.setF({0, p.getM() * gravity, 0});
  }

  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);
  if (!lc) {
    throw std::runtime_error("LennardJonesForce::calculateFLinkedCell requires LinkedCellParticleContainer");
  }
  lc->handleOutflowBoundaries();

  // Linked Cells iteration with N3L
  lc->iteratePairs([&](Particle& p_i, Particle& p_j) {
    const auto dist = p_j.getX() - p_i.getX();
    double term = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

    if (term > cutoffRadiusSq)
      return;

    term = 1 / term;

    const int idx = p_i.getType() * tableWidth + p_j.getType();
    term = pairLookupTable1[idx] * term * term * term * term * (pairLookupTable2[idx] - term * term * term);

    const auto F_vec = term * dist;
    p_i.setF(p_i.getF() + F_vec);
    p_j.setF(p_j.getF() - F_vec);
  });

  applyReflectiveBoundaries(lc);
  applyPeriodicBoundaries(lc);
}

void LennardJonesForce::applyReflectiveBoundaries(const LinkedCellParticleContainer* lc) const {

  const auto domainOrigin = lc->domain_origin();
  const auto domainDims = lc->domain_dims();
  const auto boundaryTypes = lc->boundary_types();

  for (auto& p : particles) {
    const auto x = p.getX();
    auto F_total = p.getF();

    // Use precomputed per-type lookup tables instead of computing pow() per particle
    const int pType = p.getType();
    const double p_repulsionDistance = repulsionDistanceLookup[pType];
    const double sigma6 = sigma6Lookup[pType];
    const double p_epsilon = epsilonLookup[pType];

    // Helper computes ghost particle repulsion force for a reflective wall
    auto computeGhostForce = [&](int d, double wallPos) -> std::array<double, 3> {
      auto ghostX = x;
      ghostX[d] = 2.0 * wallPos - x[d];  // Ghost particle is mirrored across the wall

      const auto ghostDist = ghostX - x;
      // Compute squared norm directly to avoid sqrt
      const double normSq = ghostDist[0] * ghostDist[0] + ghostDist[1] * ghostDist[1] + ghostDist[2] * ghostDist[2];
      const double repDistSq = p_repulsionDistance * p_repulsionDistance;

      if (normSq < repDistSq && normSq > 0.) {
        const double inv_norm2 = 1.0 / normSq;
        const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;
        const double crossing_norm_quot_6 = sigma6 * inv_norm6;
        const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

        return (24.0 * p_epsilon * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12)) * ghostDist;
      }
      return {0., 0., 0.};
    };

    for (int d = 0; d < 3; ++d) {
      const double minD = domainOrigin[d];
      const double maxD = domainOrigin[d] + domainDims[d];

      // Handling two opposite boundaries per dimension d -> 6 faces
      if (boundaryTypes[2 * d] == BoundaryType::REFLECTIVE) {
        if (const double distToWall = x[d] - minD; distToWall > 0. && distToWall < cutoffRadius) {
          F_total = F_total + computeGhostForce(d, minD);
        }
      }
      if (boundaryTypes[2 * d + 1] == BoundaryType::REFLECTIVE) {
        if (const double distToWall = maxD - x[d]; distToWall > 0. && distToWall < cutoffRadius) {
          F_total = F_total + computeGhostForce(d, maxD);
        }
      }
    }
    p.setF(F_total);
  }
}

void LennardJonesForce::calcFPeriodicBoundary(Particle* p1, Particle* p2) const {
  // Use optimized lookup tables for mixed sigma/epsilon values
  std::array<double, 3> dist = {p2->getX()[0] - p1->getX()[0], p2->getX()[1] - p1->getX()[1],
                                p2->getX()[2] - p1->getX()[2]};
  
  // Use squared distance to avoid sqrt
  double term = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];

  if (term >= cutoffRadiusSq)
    return;

  term = 1.0 / term;

  // Use precomputed lookup tables (same as calculateFLinkedCell)
  const int idx = p1->getType() * tableWidth + p2->getType();
  term = pairLookupTable1[idx] * term * term * term * term * (pairLookupTable2[idx] - term * term * term);

  const auto F_vec = term * dist;
  p1->setF(p1->getF() + F_vec);
  p2->setF(p2->getF() - F_vec);
}

void LennardJonesForce::applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const {
  const auto domainDims = lc->domain_dims();
  const auto boundaryTypes = lc->boundary_types();
  const auto numCells = lc->num_cells();

  // Face interactions
  // Handle particle pairs across opposite faces of the domain
  // X-periodic: left wall (x=1) <-> right wall (x=numCells[0]-2)
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC) {
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++)
      for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
        auto& cell1 = lc->cell_at(1, c1y, c1z);
        // Iterate through all neighboring cells of cell1 on opposing wall
        for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++)
          for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
            if (c2y < 1 || c2y > numCells[1] - 2 || c2z < 1 || c2z > numCells[2] - 2)
              continue;

            auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, c2z);
            // Iterate through all particle pairs between cells and apply force function
            for (auto& p1 : cell1)
              for (auto& p2 : cell2) {
                p2->setX(p2->getX()[0] - domainDims[0], 0);
                calcFPeriodicBoundary(p1, p2);
                p2->setX(p2->getX()[0] + domainDims[0], 0);
              }
          }
      }
  }

  // Y-periodic: bottom wall (y=1) <-> top wall (y=numCells[1]-2)
  if (boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++)
      for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
        auto& cell1 = lc->cell_at(c1x, 1, c1z);
        // Iterate through all neighboring cells of cell1 on opposing wall
        for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++)
          for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {

            if (c2x < 1 || c2x > numCells[0] - 2 || c2z < 1 || c2z > numCells[2] - 2)
              continue;

            auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, c2z);
            // Iterate through all particle pairs between cells and apply force function
            for (auto& p1 : cell1)
              for (auto& p2 : cell2) {
                p2->setX(p2->getX()[1] - domainDims[1], 1);
                calcFPeriodicBoundary(p1, p2);
                p2->setX(p2->getX()[1] + domainDims[1], 1);
              }
          }
      }
  }

  // Z-periodic: front wall (z=1) <-> back wall (z=numCells[2]-2)
  if (boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++)
      for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
        auto& cell1 = lc->cell_at(c1x, c1y, 1);
        // Iterate through all neighboring cells of cell1 on opposing wall
        for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++)
          for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {

            if (c2x < 1 || c2x > numCells[0] - 2 || c2y < 1 || c2y > numCells[1] - 2)
              continue;

            auto& cell2 = lc->cell_at(c2x, c2y, numCells[2] - 2);
            // Iterate through all particle pairs between cells and apply force function
            for (auto& p1 : cell1)
              for (auto& p2 : cell2) {
                p2->setX(p2->getX()[2] - domainDims[2], 2);
                calcFPeriodicBoundary(p1, p2);
                p2->setX(p2->getX()[2] + domainDims[2], 2);
              }
          }
      }
  }

  // Edge interaction: handle particle pairs along edges (two periodic dimensions)
  // X+y periodic: bottom-left edge <-> top-right edge, top-left edge <-> bottom-right edge
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC) {
    for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
      auto& cell1 = lc->cell_at(1, 1, c1z);
      // Iterate through all neighboring cells of cell1 on opposing side
      for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {

        if (c2z < 1 || c2z > numCells[0] - 2)
          continue;

        auto& cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, c2z);
        // Iterate through all particle pairs between cells and apply force function
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
    // Iterate through all boundary cells of top left edge
    for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
      auto& cell1 = lc->cell_at(1, numCells[1] - 2, c1z);
      // Iterate through all neighboring cells of cell1 on opposing side
      for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {

        if (c2z < 1 || c2z > numCells[0] - 2)
          continue;

        auto& cell2 = lc->cell_at(numCells[0] - 2, 1, c2z);
        // Iterate through all particle pairs between cells
        // Move particle temporarily and apply force function
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

  // X+z periodic: front-left edge <-> back-right edge, back-left edge <-> front-right edge
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      auto& cell1 = lc->cell_at(1, c1y, 1);
      // Iterate through all neighboring cells of cell1 on opposing side
      for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {

        if (c2y < 1 || c2y > numCells[1] - 2)
          continue;

        auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, numCells[2] - 2);
        // Iterate through all particle pairs between cells and apply force function
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
    // Iterate through all boundary cells of back left edge
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      auto& cell1 = lc->cell_at(1, c1y, numCells[2] - 2);
      // Iterate through all neighboring cells of cell1 on opposing side
      for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {

        if (c2y < 1 || c2y > numCells[1] - 2)
          continue;

        auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, 1);
        // Iterate through all particle pairs between cells and apply force function
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

  // Y+z periodic: bottom-front edge <-> top-back edge, bottom-back edge <-> top-front edge
  if (boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      auto& cell1 = lc->cell_at(c1x, 1, 1);
      // Iterate through all neighboring cells of cell1 on opposing side
      for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {

        if (c2x < 1 || c2x > numCells[0] - 2)
          continue;

        auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, numCells[2] - 2);
        // Iterate through all particle pairs between cells and apply force function
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
    // Iterate through all boundary cells of bottom back edge
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      auto& cell1 = lc->cell_at(c1x, 1, numCells[2] - 2);
      // Iterate through all neighboring cells of cell1 on opposing side
      for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {

        if (c2x < 1 || c2x > numCells[0] - 2)
          continue;

        auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, 1);
        // Iterate through all particle pairs between cells and apply force function
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

  // Corner interactions
  // Handle particle pairs between all 8 corners in fully periodic 3d domain
  // 4 corner pairs: (0,0,0)<->(1,1,1), (1,0,0)<->(0,1,1), (0,0,1)<->(1,1,0), (1,0,1)<->(0,1,0)
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

  // Pre-acclocate lookup tables for reflective boundaries
  repulsionDistanceLookup.clear();
  sigma6Lookup.clear();
  epsilonLookup.clear();
  repulsionDistanceLookup.resize(tableWidth, -1.0);
  sigma6Lookup.resize(tableWidth, -1.0);
  epsilonLookup.resize(tableWidth, -1.0);

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
      const double sigma = p.getSigma();
      const double sigma2 = sigma * sigma;
      const double sigma6 = sigma2 * sigma2 * sigma2;
      // Precomputed 2^(1/6) 
      repulsionDistanceLookup[t] = 1.1224620483093729814335330496791795162324111106139867534404095458 * sigma;
      sigma6Lookup[t] = sigma6;
      epsilonLookup[t] = p.getEpsilon();
    }
  }

  // For every pair of unique types, compute and store lookup values
  for (const int ti : uniqueTypes) {
    for (const int tj : uniqueTypes) {
      const int idx = ti * tableWidth + tj;
      if (pairLookupTable1[idx] < 0) {
        // Apply Lorentz-Berthelot mixing rules
        const double sigma_ij = (typeSigma[ti] + typeSigma[tj]) / 2.0;
        const double epsilon_ij = std::sqrt(typeEpsilon[ti] * typeEpsilon[tj]);

        // Precompute: 48 * epsilon_ij * sigma_ij^12
        const double sigma2 = sigma_ij * sigma_ij;
        const double sigma4 = sigma2 * sigma2;
        const double sigma6 = sigma4 * sigma2;
        const double sigma12 = sigma6 * sigma6;
        const double eps_48_sigma_12 = 48.0 * epsilon_ij * sigma12;

        // Precompute: 1 / (2 * sigma_ij^6)
        const double sigma_m6_2 = 1.0 / (2.0 * sigma6);

        pairLookupTable1[idx] = eps_48_sigma_12;
        pairLookupTable2[idx] = sigma_m6_2;

        // Symmetric: save to (tj, ti) as well
        const int idxSym = tj * tableWidth + ti;
        pairLookupTable1[idxSym] = eps_48_sigma_12;
        pairLookupTable2[idxSym] = sigma_m6_2;
      }
    }
  }

  SPDLOG_DEBUG("Precomputed constants for {} unique particle types ({} pair combinations)",
               uniqueTypes.size(), uniqueTypes.size() * uniqueTypes.size());
}
