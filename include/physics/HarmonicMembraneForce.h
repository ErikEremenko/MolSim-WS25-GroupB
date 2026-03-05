#pragma once

#include "physics/ForceCalc.h"

/**
 * @class HarmonicMembraneForce
 * @brief Harmonic potential for membrane bonds between neighboring particles
 * Models interactions between direct and diagonal neighbors in a 2D membrane.
 */
class HarmonicMembraneForce final : public ForceCalc {
 private:
  double stiffness;           ///< Stiffness constant k
  double avgBondLength;       ///< Average bond length r0 for direct neighbors
  double diagonalBondLength;  ///< Bond length for diagonal neighbors: sqrt(2) * r0

 public:
  /**
   * @param particles ParticleContainer with membrane particles (must have neighbor info set)
   * @param k Stiffness constant
   * @param r0 Average bond length for direct neighbors
   */
  HarmonicMembraneForce(ParticleContainer& particles, double k, double r0);

  void calculateF() override;
};
