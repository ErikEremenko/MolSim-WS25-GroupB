#include "physics/thermodynamics/ThermodynamicsStatistics.h"

ThermodynamicsStatistics::ThermodynamicsStatistics(ParticleContainer& particles, Thermostat* thermostat,
                                                   const std::array<double, 3>& domainDims, const std::string& baseName)
    : msd(particles), rdf(particles, domainDims), thermostat(thermostat), writer(baseName) {}

void ThermodynamicsStatistics::updateStatistics() {
  // Calculate thermodynamical statistics and temperature
  const auto diffusion = msd.calculateDiffusion();
  const auto& rdfDensities = rdf.calculateDistribution();

  double temperature = 0.0;
  if (thermostat) {  // this is allowed to be nullptr but will practically never be
    temperature = thermostat->calculateCurrentTemperature();
  }

  const int simulationIterations = ++writeIteration * updateFrequency;
  writer.write(simulationIterations, temperature, diffusion, rdfDensities);
}
