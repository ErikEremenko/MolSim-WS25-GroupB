#pragma once

#include "physics/ForceCalc.h"

/**
 * @class GlobalGravityForce
 * @brief Applies a constant gravitational acceleration to all particles along a specified axis
 */
class GlobalGravityForce final : public ForceCalc {
 private:
  double gravity;  ///< Gravitational acceleration magnitude
  int axis;        ///< Axis: 0=x, 1=y, 2=z

 public:
  /**
   * @param particles ParticleContainer that stores the particles
   * @param g Gravity acceleration value (e.g., -9.81 for downward gravity)
   * @param axis Axis direction: 0=x, 1=y (default), 2=z
   */
  GlobalGravityForce(ParticleContainer& particles, double g, int axis = 1);

  void calculateF() override;
};
