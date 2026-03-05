#pragma once

#include "physics/ForceCalc.h"

/**
 * @class TruncatedLJForce
 * @brief Repulsive-only Lennard-Jones potential, truncated at 2^(1/6)·sigma
 * Used for membrane simulations to prevent self-penetration without attraction.
 * Uses Lorentz-Berthelot mixing rules for mixed particle types.
 */
class TruncatedLJForce final : public ForceCalc {
 public:
  explicit TruncatedLJForce(ParticleContainer& particles);

  void calculateF() override;
};
