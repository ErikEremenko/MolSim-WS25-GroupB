#pragma once

/**
 * @class Thermostat
 * Regulates the temperature of a given system.
 */

#include <limits>

#include "ParticleContainer.h"

class Thermostat {
private:
  /**
   * @brief Particles controlled by this thermostat.
   */
  ParticleContainer& particles;
  /**
   * @brief The number of time steps after which the thermostat is periodically applied.
   */
  const int nThermostat;

  /**
   * @brief Target temperature to reach
   */
  double tempTarget;
  /**
   * @brief Maximal absolute temperature change allowed for one application of the thermostat.
   */
  double tempDelta;

  double calculateCurrentTemperature();
public:
  ~Thermostat();
  Thermostat(ParticleContainer& particles, int nThermostat, double tempTarget, double tempDelta = std::numeric_limits<double>::infinity());

  /**
   * @brief Initializes the temperature of the system using Maxwell-Boltzmann velocities in Anderson thermostat style.
   * Does not preserve previous velocities.
   *
   * @param tempInit The initial temperature of the system
   */
  void initializeTemperature(double tempInit);
  /**
   * @brief Sets the temperature of the system using gradual scaling.
   * Preserves previous velocity directions using a scaling factor.
   *
   * @param tempNew New temperature of the system
   */
  void setTemperature(double tempNew);
  /**
   * @brief Gradually scales the . This function should be called every @nThermoStat iterations.
   */
  void updateTemperature();  // TODO: Make this function inline

  /**
   * @return The update frequency of the thermostat, see @ref nThermostat.
   */
  [[nodiscard]] int getUpdateFrequency() const;
};
