#include "physics/thermodynamics/ThermodynamicsStatistics.h"

ThermodynamicsStatistics::ThermodynamicsStatistics(ParticleContainer& particles,
                                                   const std::array<double, 3>& domainDims, const std::string& baseName)
    : msd(particles), rdf(particles, domainDims), writer(baseName) {}

void ThermodynamicsStatistics::updateStatistics() {
  const auto diffusion = msd.calculateDiffusion();
  const auto& rdfDensities = rdf.calculateDistribution();

  writer.write(diffusion, rdfDensities);
}
