#include "Thermostat.h"

#include <algorithm>
#include <cmath>

#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"

double Thermostat::calculateCurrentTemperature() {
  double totalKineticEnergyTimesTwo = 0.0;
  for (const auto& particle : particles) {
    totalKineticEnergyTimesTwo += particle.getM() * ArrayUtils::squaredL2Norm(particle.getV());
  }

  return totalKineticEnergyTimesTwo / (3 * particles.size());  // number of dimensions = 3 in our case
}

Thermostat::~Thermostat() = default;

Thermostat::Thermostat(ParticleContainer& particles, int nThermostat, double tempTarget, double tempDelta)
    : particles(particles), nThermostat(nThermostat), tempTarget(tempTarget), tempDelta(tempDelta) {}

void Thermostat::initializeTemperature(double tempInit) {
  for (auto& particle : particles) {
    const std::array<double, 3> v = maxwellBoltzmannDistributedVelocity(std::sqrt(tempInit / particle.getM()), 3);
    particle.setV(v);
  }
}

void Thermostat::setTemperature(double tempNew) {
  const double tempCurrent = calculateCurrentTemperature();
  const double beta = std::sqrt(tempNew / tempCurrent);  // scaling factor
  for (auto& particle : particles) {
    particle.setV(beta * particle.getV());
  }
}

void Thermostat::updateTemperature() {
  const double tempCurrent = calculateCurrentTemperature();
  const double tempNew = tempCurrent + std::clamp(tempTarget - tempCurrent, -tempDelta, tempDelta);

  const double beta = std::sqrt(tempNew / tempCurrent);  // scaling factor
  for (auto& particle : particles) {
    particle.setV(beta * particle.getV());
  }
}

// Getters
int Thermostat::getUpdateFrequency() const {
  return nThermostat;
}
