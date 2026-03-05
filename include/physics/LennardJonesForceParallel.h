#pragma once

#include "physics/ForceCalc.h"

/**
 * @class LennardJonesForceParallel
 * @brief Models the Lennard-Jones potential with parallelization
 */
class LennardJonesForceParallel final : public ForceCalc {
 private:
  const double epsilon, sigma, cutoffRadius;

 public:
  LennardJonesForceParallel(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius);
  void calculateF() override;
};
