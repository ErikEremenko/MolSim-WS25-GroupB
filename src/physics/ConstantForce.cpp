#include "physics/ConstantForce.h"

#include <utility>

#include "utils/ArrayUtils.h"

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
