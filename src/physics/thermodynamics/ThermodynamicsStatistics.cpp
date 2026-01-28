#include "physics/thermodynamics/ThermodynamicsStatistics.h"

ThermodynamicsStatistics::ThermodynamicsStatistics(
  ParticleContainer& particles, const std::array<double, 3>& domainDims, const int updateFrequency = 1000
  ) : msd(particles), rdf(particles, domainDims), updateFrequency(updateFrequency) {}

int ThermodynamicsStatistics::getUpdateFrequency() const {
  return updateFrequency;
}

void ThermodynamicsStatistics::updateStatistics() {
  double diffusion = msd.calculateDiffusion();

  rdf.calculateDistribution();
  auto& rdfBins = rdf.getIntervals();

  // TODO: Add file output
}
