#pragma once

/**
 * @class Thermostat
 * Regulates the temperature of a given system.
 */

#include <limits>

#include "physics/ParticleContainer.h"

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

 public:
  ~Thermostat();
  Thermostat(ParticleContainer& particles, int nThermostat, double tempTarget = 1.0,
             double tempDelta = std::numeric_limits<double>::infinity());

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
   * @brief A single application of the thermostat using gradual temperature scaling.
   * This function should be called every @nThermoStat iterations.
   */
  void updateTemperature();  // TODO: Make this function inline
  /**
   * Calculates the current temperature of the system by calculating the kinetic energy of the particles.
   * @return The current temperature of the system
   */
  double calculateCurrentTemperature();

  // Getters and setters
  /**
   * @brief Gets the frequency of thermostat applications
   * @return nThermostat member
   */
  [[nodiscard]] int getUpdateFrequency() const;

  /**
   * @brief Sets the target temperature of the system.
   * @note Used in @ref Simulation to set the target temperature to the current calculated temperature
   * if no initial temperature or target temperature was specified.
   * @param newTempTarget
   */
  void setTargetTemperature(double newTempTarget);
};
