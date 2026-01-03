#pragma once

/**
 * @class Thermostat
 * Regulates the temperature of a given system.
 */

#include "ParticleContainer.h"

class Thermostat {
private:
  float tempCurrent;
  float tempTarget;
  const int nThermostat;
  float tempDiff;

  bool heating;
  bool holding;
public:
  Thermostat(float tempInit, int nThermostat, float tempTarget, float tempDiff);
  Thermostat(float tempInit, int nThermostat);

  inline void updateTemperature(ParticleContainer& particles);

  ~Thermostat();

  [[nodiscard]] int getNThermostat() const;
};
