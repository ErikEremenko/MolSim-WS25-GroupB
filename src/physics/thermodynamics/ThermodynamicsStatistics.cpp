#include "physics/thermodynamics/ThermodynamicsStatistics.h"

ThermodynamicsStatistics::ThermodynamicsStatistics(ParticleContainer& particles, const int updateFrequency = 1000)
  : msd(particles), rdf(particles), updateFrequency(updateFrequency) {}

int ThermodynamicsStatistics::getUpdateFrequency() const {
  return updateFrequency;
}

void ThermodynamicsStatistics::updateStatistics() {
  double diffusion = msd.calculateDiffusion();
  rdf.calculateDistribution();

  // TODO: Add file output
}
