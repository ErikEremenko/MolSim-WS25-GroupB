#pragma once

#include "physics/ForceCalc.h"

#include <array>
#include <utility>
#include <vector>

/**
 * @class ConstantForce
 * @brief Applies a constant force to specific particles (identified by membrane x/y indices)
 * The force is only applied until a specified end time ("pulling" membrane particles)
 */
class ConstantForce final : public ForceCalc {
 private:
  std::array<double, 3> force;
  double endTime;
  double& currentTime;                             ///< Reference to simulation's current time
  std::vector<std::pair<int, int>> targetIndices;  ///< x/y indices of target particles
  int membraneDimY;                                ///< Y-dimension of membrane grid for index calculation

 public:
  /**
   * @param particles ParticleContainer
   * @param fx Force in x direction
   * @param fy Force in y direction
   * @param fz Force in z direction
   * @param endTime Time after which force stops being applied
   * @param currentTime Reference to the simulation's current time variable
   * @param targetIndices Vector of (x, y) index pairs identifying which particles to pull
   * @param membraneDimY Y-dimension of the membrane for calculating particle indices
   */
  ConstantForce(ParticleContainer& particles, double fx, double fy, double fz, double endTime, double& currentTime,
                std::vector<std::pair<int, int>> targetIndices, int membraneDimY);

  void calculateF() override;
};
