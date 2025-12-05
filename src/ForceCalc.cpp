#include "ForceCalc.h"
#include <math.h>
#include <spdlog/spdlog.h>
#include "utils/ArrayUtils.h"

#include "LinkedCellParticleContainer.h"

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

void GravityForce::calculateF() {
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
      if (norm == 0.) {
        // avoid division by zero
        SPDLOG_ERROR(
            "Calculated a zero norm between particles. This is likely caused "
            "by an incorrect initialization of the Simulation.");
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

LennardJonesForce::LennardJonesForce(ParticleContainer& particles, const double epsilon, const double sigma,
                                     const double cutoffRadius)
    : LennardJonesForce(
          particles, epsilon, sigma, cutoffRadius, std::pow(2.0, 1.0 / 6.0) * sigma  // default repulsion distance
      ) {}
LennardJonesForce::LennardJonesForce(ParticleContainer& particles, const double epsilon, const double sigma,
                                     const double cutoffRadius, const double repulsionDistance)
    : ForceCalc(particles),
      epsilon(epsilon),
      sigma(sigma),
      cutoffRadius(cutoffRadius),
      repulsionDistance(repulsionDistance) {}

void LennardJonesForce::calculateF() {
  if (dynamic_cast<LinkedCellParticleContainer*>(&particles)) {
    calculateFLinkedCell();
  } else {  // O(n^2) implementation
    for (auto& p : particles) {
      p.setF({});
    }
    const double sigma2 = sigma * sigma;
    const double sigma6 = sigma2 * sigma2 * sigma2;

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

        const auto F_vector =
            (24.0 * epsilon) * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12) * dist;

        // apply forces using Newton's third law (O(n^2) -> O(((n^2)/2))
        auto F_i = p_i.getF();
        auto F_j = p_j.getF();
        // actio est reactio
        p_i.setF(F_i + F_vector);
        p_j.setF(F_j - F_vector);
      }
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
    // index offset for Newton's third law
    for (size_t j = i + 1; j < n_particles; ++j) {
      auto& p_i = particles[i];
      auto& p_j = particles[j];

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
      // apply forces using Newton's third law and atomic operations
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
    p.setF({});
  }

  auto* lc = dynamic_cast<LinkedCellParticleContainer*>(&particles);
  if (!lc) {
    throw std::runtime_error("LennardJonesForce::calculateFLinkedCell requires LinkedCellParticleContainer");
  }
  lc->updateCells();

  const double sigma2 = sigma * sigma;
  const double sigma6 = sigma2 * sigma2 * sigma2;

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
    } else if (norm >= cutoffRadius)
      return;

    const double inv_norm2 = 1.0 / (norm * norm);
    const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;

    const double crossing_norm_quot_6 = sigma6 * inv_norm6;
    const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

    const auto F_vec = (24.0 * epsilon * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12)) * dist;
    p_i.setF(p_i.getF() + F_vec);
    p_j.setF(p_j.getF() - F_vec);
  });

  const auto domainOrigin = lc->domain_origin();
  const auto domainDims = lc->domain_dims();
  const auto boundaryTypes = lc->boundary_types();

  for (auto& p : particles) {
    const auto x = p.getX();
    auto F_total = p.getF();

    for (int d = 0; d < 3; ++d) {
      const double minD = domainOrigin[d];
      const double maxD = domainOrigin[d] + domainDims[d];

      // handling two opposite boundaries per dimension d -> 6 faces
      if (boundaryTypes[2 * d] == LinkedCellParticleContainer::BoundaryType::REFLECTIVE) {
        if (const double distToWall = x[d] - minD; distToWall > 0. && distToWall < cutoffRadius) {
          // use ghost / virtual particle for repulsion
          auto ghostX = x;
          ghostX[d] = 2.0 * minD - x[d];  // ghost particle is mirrored across the wall

          const auto ghostDist = ghostX - x;

          if (const double norm = ArrayUtils::L2Norm(ghostDist);
              norm < repulsionDistance && norm > 0) {  // avoid division by zero
            const double inv_norm2 = 1.0 / (norm * norm);
            const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;
            const double crossing_norm_quot_6 = sigma6 * inv_norm6;
            const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

            const auto F_wall =
                (24.0 * epsilon * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12)) * ghostDist;
            F_total = F_total + F_wall;
          }
        }
      }
      if (boundaryTypes[2 * d + 1] == LinkedCellParticleContainer::BoundaryType::REFLECTIVE) {
        if (const double distToWall = maxD - x[d]; distToWall > 0. && distToWall < cutoffRadius) {
          // use ghost / virtual particle for repulsion
          auto ghostX = x;
          ghostX[d] = 2.0 * maxD - x[d];  // ghost particle is mirrored across the wall

          const auto ghostDist = ghostX - x;

          if (const double norm = ArrayUtils::L2Norm(ghostDist);
              norm < repulsionDistance && norm > 0.) {  // avoid division by zero
            const double inv_norm2 = 1.0 / (norm * norm);
            const double inv_norm6 = inv_norm2 * inv_norm2 * inv_norm2;
            const double crossing_norm_quot_6 = sigma6 * inv_norm6;
            const double crossing_norm_quot_12 = crossing_norm_quot_6 * crossing_norm_quot_6;

            const auto F_wall =
                (24.0 * epsilon * inv_norm2 * (crossing_norm_quot_6 - 2.0 * crossing_norm_quot_12)) * ghostDist;
            F_total = F_total + F_wall;
          }
        }
      }
    }
    p.setF(F_total);
  }
}
