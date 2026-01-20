#pragma once

#include <yaml-cpp/yaml.h>

#include <string>
#include <optional>

#include "simulation/SimulationConfig.h"

/**
 * @class YAMLFileReader
 * @brief Reads a YAML file and extracts a simulation configuration.
 *
 */
class YAMLFileReader {  // TODO: Add 'container' key and its parsing
 public:
  /**
   * @brief Loads the YAML file immediately.
   * @param filename Path to the YAML file.
   */
  explicit YAMLFileReader(const std::string& filename);

  /**
   * @brief Parses the simulation configuration from the loaded YAML nodes.
   * @return Simulation configuration
   */
  SimulationConfig getConfig();

 private:
  /**
   * @brief Stores loaded YAML structure.
   */
  YAML::Node config;

  /**
   * @brief Validate configuration keys.
   */
  void checkRequiredKeys() const;

  // Getters for simulation parameters
  /**
 * @brief Gets the base name for output files.
 * @return The base name string used for output file generation.
 */
  std::string getOutputBaseName() const;

  /**
   * @brief Gets the frequency at which simulation output is written.
   * @return The write frequency in number of iterations.
   */
  int getWriteFrequency() const;

  /**
   * @brief Gets the frequency at which checkpoint files are written.
   * @return The checkpoint frequency in number of iterations.
   */
  int getCheckpointFrequency() const;

  /**
   * @brief Gets the simulation end time.
   * @return The total simulation time in seconds.
   */
  double getTEnd() const;

  /**
   * @brief Gets the simulation time step size.
   * @return The time step in seconds.
   */
  double getDeltaT() const;

  /**
   * @brief Gets the number of dimensions for the simulation (2 or 3).
   * @return The number of dimensions.
   */
  std::optional<int> getDimensions() const;

  /**
   * @brief Gets the simulation domain size.
   * @return Array containing domain dimensions [x, y, z].
   */
  std::array<double, 3> getDomainSize() const;

  /**
   * @brief Gets the boundary condition types as raw strings.
   * @return Array of boundary type strings in order: [x_min, x_max, y_min, y_max, z_min, z_max].
   */
  std::array<std::string, 6> getBoundaryTypesRaw() const;

  /**
   * @brief Parses the 'thermostat' YAML block.
   * @return Optional ThermostatConfig. Returns nullopt if block is missing.
   */
  std::optional<ThermostatConfig> getThermostatConfig() const;

  // Checkpoint-related getters
  /** @brief Returns true if this file contains checkpoint particle data */
  bool isCheckpoint() const;
  /** @brief Returns the iteration number from checkpoint (0 if not a checkpoint) */
  int getCheckpointIteration() const;
  /** @brief Returns the simulation time from checkpoint (0.0 if not a checkpoint) */
  double getCheckpointTime() const;

  /**
   * @brief Reads the force type from the given string.
   * @return The parsed force type
   */
  static ForceType getForceType(const std::string& str);
};
