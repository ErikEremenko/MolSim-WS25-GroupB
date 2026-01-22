#include "physics/ForceCalc.h"

#include <cmath>
#include <spdlog/spdlog.h>

#include "physics/LinkedCellParticleContainer.h"
#include "utils/ArrayUtils.h"

constexpr double sqrtTwo = 1.4142135623730951;  // value of std::sqrt(2)

ForceCalc::ForceCalc(ParticleContainer& particles) : particles(particles) {}
ForceCalc::~ForceCalc() = default;

void ForceCalc::calculateX(ParticleContainer& particles, const double dt) {
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

void ForceCalc::calculateV(ParticleContainer& particles, const double dt) {
  for (auto& p : particles) {
    const auto v_curr = p.getV();
    const double m_i = p.getM();
    const auto F = p.getF();
    const auto F_old = p.getOldF();
    const auto v_new = v_curr + (dt / (2.0 * m_i)) * (F_old + F);
    p.setV(v_new);
  }
}

void GravityForce::calculateF() {
  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // Index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      if (norm == 0.) {
        // Avoid division by zero
        SPDLOG_ERROR(
            "Calculated a zero norm between particles. This is likely caused "
            "by an incorrect initialization of the Simulation.");
        throw std::overflow_error(
            "Calculated a zero norm between particles. This is likely caused "
            "by an incorrect initialization of the Simulation.");
      }
      const double norm3 = norm * norm * norm;

      const auto F_vector = ((p_i.getM() * p_j.getM()) / norm3) * dist;

      // Apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
      // Actio est reactio
      p_i.addF(F_vector);
      p_j.addF(-F_vector);
    }
  }
}

LennardJonesForce::LennardJonesForce(
  ParticleContainer& particles, const double epsilon, const double sigma, const double cutoffRadius, const bool isTruncated = false
  ) : ForceCalc(particles),
      epsilon(epsilon),
      sigma(sigma),
      cutoffRadius(cutoffRadius),
      repulsionDistance(std::pow(2.0, 1.0 / 6.0) * sigma),
      isTruncated(isTruncated) {}

void LennardJonesForce::calculateF() {
  if (dynamic_cast<LinkedCellParticleContainer*>(&particles)) {
    calculateFLinkedCell();
  } else {
    calculateFDirectSum();
  }
}

void LennardJonesForce::calculateFDirectSum() {
  const double sigma2 = sigma * sigma;
  const double sigma6 = sigma2 * sigma2 * sigma2;

  const size_t n_particles = particles.size();
  for (size_t i = 0; i < n_particles; ++i) {
    // Index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      if (norm == 0) {
        // Avoid division by zero
        SPDLOG_ERROR(
            "Calculated a zero norm between particles. This is likely caused by an incorrect initialization of the "
            "Simulation.");
        throw std::overflow_error(
            "Calculated a zero norm between particles. This is likely caused by an incorrect initialization of the "
            "Simulation.");
      } else if (norm >= cutoffRadius || (isTruncated && norm >= repulsionDistance)) {
        continue;
      }
      const double inv_norm2 = 1.0 / (norm * norm);
      const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

      const double crossing_norm_quot_6 = sigma6 * inv_norm6;
      const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

      const auto F_vector = (24.0 * epsilon) * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12) * dist;

      // Apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
      // Actio est reactio
      p_i.addF(F_vector);
      p_j.addF(-F_vector);
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
  const size_t n_particles = particles.size();

#pragma omp parallel for schedule(guided)
  for (size_t i = 0; i < n_particles; ++i) {
    // Index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      const auto dist = p_j.getX() - p_i.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      if (norm == 0) {
        // Avoid division by zero
        SPDLOG_ERROR(
            "Calculated a zero norm between particles. This is likely caused "
            "by an incorrect initialization of the Simulation.");
        throw std::overflow_error(
            "Calculated a zero norm between particles. This is likely caused "
            "by an incorrect initialization of the Simulation.");
      } else if (norm >= cutoffRadius) {
        continue;
      }
      const double inv_norm2 = 1.0 / (norm * norm);
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

  // Linked Cells iteration with N3L
  lc->iteratePairs([&](Particle& p_i, Particle& p_j) {
    const auto dist = p_j.getX() - p_i.getX();
    const double norm = ArrayUtils::L2Norm(dist);

    if (norm == 0) {
      // avoid division by zero
      SPDLOG_ERROR(
          "Calculated a zero norm between particles. This is likely caused "
          "by an incorrect initialization of the Simulation.");
      throw std::overflow_error(
          "Calculated a zero norm between particles. This is likely caused "
          "by an incorrect initialization of the Simulation.");
    } else if (norm >= cutoffRadius || (isTruncated && norm >= repulsionDistance)) {
      return;
    }

    // Get the particle's sigma/epsilon and apply mixing rule
    const auto sigma_ij = (p_i.getSigma() + p_j.getSigma()) / 2;
    const auto epsilon_ij = std::sqrt(p_i.getEpsilon() * p_j.getEpsilon());

    const double sigma2 = sigma_ij * sigma_ij;
    const double sigma6 = sigma2 * sigma2 * sigma2;
    const double inv_norm2 = 1.0 / (norm * norm);
    const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

    const double crossing_norm_quot_6 = sigma6 * inv_norm6;
    const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

    const auto F_vec = (24.0 * epsilon_ij * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12)) * dist;
    p_i.addF(F_vec);
    p_j.addF(-F_vec);
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
    std::array<double, 3> F_sum = {0.0, 0.0, 0.0};

    // Use particle's own sigma/epsilon
    const double p_sigma = p.getSigma();
    const double p_epsilon = p.getEpsilon();
    const double p_repulsionDistance = std::pow(2.0, 1.0 / 6.0) * p_sigma;
    const double sigma2 = p_sigma * p_sigma;
    const double sigma6 = sigma2 * sigma2 * sigma2;

    // Helper computes ghost particle repulsion force for a reflective wall
    auto computeGhostForce = [&](int d, double wallPos) -> std::array<double, 3> {
      auto ghostX = x;
      ghostX[d] = 2.0 * wallPos - x[d];  // Ghost particle is mirrored across the wall

      const auto ghostDist = ghostX - x;
      const double norm = ArrayUtils::L2Norm(ghostDist);

      if (norm < p_repulsionDistance && norm > 0.) {  // Avoid division by zero
        const double inv_norm2 = 1.0 / (norm * norm);
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
          F_sum = F_sum + computeGhostForce(d, minD);
        }
      }
      if (boundaryTypes[2 * d + 1] == BoundaryType::REFLECTIVE) {
        if (const double distToWall = maxD - x[d]; distToWall > 0. && distToWall < cutoffRadius) {
          F_sum = F_sum + computeGhostForce(d, maxD);
        }
      }
    }
    p.addF(F_sum);
  }
}

void LennardJonesForce::calcFPeriodicBoundary(Particle* p1, Particle* p2) const {

  const double sigma2 = sigma * sigma;
  const double sigma6 = sigma2 * sigma2 * sigma2;

  std::array<double, 3> dist = {p2->getX()[0] - p1->getX()[0], p2->getX()[1] - p1->getX()[1],
                                p2->getX()[2] - p1->getX()[2]};
  const double norm = ArrayUtils::L2Norm(dist);

  if (norm == 0) {
    // avoid division by zero
    SPDLOG_ERROR(
        "Calculated a zero norm between particles. This is likely caused "
        "by an incorrect initialization of the Simulation.");
    throw std::overflow_error(
        "Calculated a zero norm between particles. This is likely caused "
        "by an incorrect initialization of the Simulation.");
  } else if (norm >= cutoffRadius)
    return;

  const double inv_norm2 = 1.0 / (norm * norm);
  const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

  const double crossing_norm_quot_6 = sigma6 * inv_norm6;
  const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

  const auto F_vec = (24.0 * epsilon * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12)) * dist;
  p1->addF(F_vec);
  p2->addF(-F_vec);
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

ConstantAccelerationForce::ConstantAccelerationForce(
  ParticleContainer& particles, const double accX, const double accY, const double accZ
  ) : ForceCalc(particles), accX(accX), accY(accY), accZ(accZ) {}

void ConstantAccelerationForce::calculateF() {
  for (auto& p : particles) {
    const double m = p.getM();
    p.addF({m*accX, m*accY, m*accZ});
  }
}

MembraneBondForce::MembraneBondForce(ParticleContainer& particles, const double stiffnessConstant, const double bondLength)
  : ForceCalc(particles), stiffnessConstant(stiffnessConstant), bondLength(bondLength) {}

void MembraneBondForce::calculateF() {
  for (auto& p : particles) {
    for (const int neighborId : p.getDirectNeighbors()) {
      Particle& neighbor = particles[neighborId];
      const auto dist = p.getX() - neighbor.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      const auto force = stiffnessConstant * (norm - bondLength) / norm * dist;
      p.addF(force);
    }

    for (const int neighborId : p.getDiagonalNeighbors()) {
      Particle& neighbor = particles[neighborId];
      const auto dist = p.getX() - neighbor.getX();
      const double norm = ArrayUtils::L2Norm(dist);
      const auto force = stiffnessConstant * (norm - sqrtTwo*bondLength) / norm * dist;
      p.addF(force);
    }
  }
}

MembraneConstantForce::MembraneConstantForce(ParticleContainer& particles, const double dt)
    : ForceCalc(particles), particlesInitialized(false), affectedParticles(), currentTime(0), dt(dt) {}

void MembraneConstantForce::calculateF() {
  if (currentTime >= forceEndTime) return;
  if (!particlesInitialized) {
    /* This has to be done here as we have lazy particle generation and the particle container is empty
     * when the force calculation strategies are being initialized.
     */
    // Get the affected particles - (17/24), (17/25), (18/24) and (18/25)
    // We use the formula: x*yDim + y for indexing, as that is how ParticleGenerator initializes the particles
    affectedParticles[0] = &particles[17*yDim + 24];
    affectedParticles[1] = &particles[17*yDim + 25];
    affectedParticles[2] = &particles[18*yDim + 24];
    affectedParticles[3] = &particles[18*yDim + 25];

    particlesInitialized = true;
  }

  for (auto p : affectedParticles) {
    p->addF({0, 0, forceValue});
  }

  currentTime += dt;
}
