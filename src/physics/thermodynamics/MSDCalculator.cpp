#include "physics/thermodynamics/MSDCalculator.h"

#include "utils/ArrayUtils.h"

MSDCalculator::MSDCalculator(ParticleContainer& particles) : particles(particles) {}

double MSDCalculator::calculateDiffusion() {
  double displacement_squared_sum = 0.0;

  for (auto& p : particles) {
    displacement_squared_sum += ArrayUtils::squaredL2Norm(p.getDisplacement());
    p.resetDisplacement();
  }

  return displacement_squared_sum / particles.size();
}
