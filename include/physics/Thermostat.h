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
  int nThermostat;

  /**
   * @brief Target temperature to reach
   */
  double tempTarget;
  /**
   * @brief Maximal absolute temperature change allowed for one application of the thermostat.
   */
  double tempDelta;
  /**
   * @brief Number of spatial dimensions (e.g., 2 for 2D, 3 for 3D)
   */
  int dimensions;

 public:
  ~Thermostat();
  /**
   * @brief Construct a new Thermostat object
   *
   * @param particles The particle container to control
   * @param nThermostat Number of steps between thermostat applications
   * @param tempTarget Target temperature
   * @param tempDelta Maximum temperature change per application
   * @param dimensions Number of spatial dimensions (default 3)
   */
  Thermostat(ParticleContainer& particles, int nThermostat, double tempTarget = 1.0,
             double tempDelta = std::numeric_limits<double>::infinity(), int dimensions = 3);

  /**
   * @brief Initializes the temperature of the system using Maxwell-Boltzmann velocities in Anderson thermostat style.
   * Does not preserve previous velocities.
   *
   * @param tempInit The initial temperature of the system
   */
  void initializeTemperature(double tempInit);
  /**
   * @brief Update the temperature of the system towards the target temperature
   */
  void updateTemperature();
  /**
   * @brief Directly set the temperature of the system (scaling velocities)
   *
   * @param tempNew The new temperature to set
   */
  void setTemperature(double tempNew);
  /**
   * @brief Calculate the current temperature of the system based on kinetic energy
   *
   * @return double Current temperature
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
