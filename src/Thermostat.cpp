#include "Thermostat.h"

#include <utils/ArrayUtils.h>
#include <cmath>

Thermostat::~Thermostat() = default;

Thermostat::Thermostat(float tempInit, int nThermostat, float tempTarget, float maxTempDiff)
    : tempCurrent(tempInit), tempTarget(tempTarget), nThermostat(nThermostat) {
  // TODO: No 'valid argument' checks done here, consider adding exception throwing
  tempDiff = (tempTarget - tempInit) / nThermostat;
  if (tempDiff >= 0) {
    heating = true;
  } else {
    heating = false;
    tempDiff = -tempDiff;
  }

  if (tempDiff >= maxTempDiff) {
    tempDiff = maxTempDiff;
  }

  if (tempDiff == 0) {
    holding = true;
  }
}

Thermostat::Thermostat(float tempInit, int nThermostat)
    : Thermostat(tempInit, nThermostat, tempInit, std::numeric_limits<float>::infinity()) {}

inline void Thermostat::updateTemperature(ParticleContainer& particles) {
  // Check and skip if we already reached the target temperature
  if ((heating && tempCurrent >= tempTarget) || (!heating && tempCurrent <= tempTarget)) {
    return;
  }

  float tempNew = tempCurrent + tempDiff;
  const float beta = std::sqrt(tempNew / tempCurrent);  // scaling factor

  for (auto& particle : particles) {
    particle.setV(beta * particle.getV());
  }

  tempCurrent = tempNew;
}

// Getters
int Thermostat::getNThermostat() const {
  return nThermostat;
}
