#include "physics/thermodynamics/ThermodynamicsStatistics.h"

ThermodynamicsStatistics::ThermodynamicsStatistics(
  ParticleContainer& particles, const std::array<double, 3>& domainDims
  ) : msd(particles), rdf(particles, domainDims) {}

void ThermodynamicsStatistics::updateStatistics() {
  double diffusion = msd.calculateDiffusion();

  rdf.calculateDistribution();
  auto& rdfBins = rdf.getIntervals();

  // TODO: Add file output
}
