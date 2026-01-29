#include "physics/thermodynamics/Thermostat.h"

#include <algorithm>
#include <cmath>

#include "utils/ArrayUtils.h"
#include "utils/MaxwellBoltzmannDistribution.h"

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG
#endif  // SPDLOG_ACTIVE_LEVEL
#include <spdlog/spdlog.h>

double Thermostat::calculateCurrentTemperature() {
  double totalKineticEnergyTimesTwo = 0.0;
  for (const auto& particle : particles) {
    totalKineticEnergyTimesTwo += particle.getM() * ArrayUtils::squaredL2Norm(particle.getV());
  }

  return totalKineticEnergyTimesTwo / (dimensions * particles.size());
}

Thermostat::~Thermostat() = default;

Thermostat::Thermostat(ParticleContainer& particles, int nThermostat, double tempTarget, double tempDelta,
                       int dimensions)
    : particles(particles),
      nThermostat(nThermostat),
      tempTarget(tempTarget),
      tempDelta(tempDelta),
      dimensions(dimensions) {}

void Thermostat::initializeTemperature(double tempInit) {
  for (auto& particle : particles) {
    const std::array<double, 3> v =
        maxwellBoltzmannDistributedVelocity(std::sqrt(tempInit / particle.getM()), dimensions);
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

  SPDLOG_DEBUG("Updated temperature from {} to {}", tempCurrent, tempNew);
}

// Getters and setters
int Thermostat::getUpdateFrequency() const {
  return nThermostat;
}

void Thermostat::setTargetTemperature(const double newTempTarget) {
  tempTarget = newTempTarget;
}
