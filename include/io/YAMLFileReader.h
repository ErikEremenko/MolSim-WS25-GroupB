#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include "utils/MaxwellBoltzmannDistribution.h"
#include "simulation/SimulationConfig.h"

/**
 * @class YAMLFileReader
 * @brief Reads a YAML file and extracts a simulation configuration.
 *
 */
class YAMLFileReader {
 public:
  /**
   * @brief Constructor loads the YAML file immediately.
   * @param filename Path to the YAML input file.
   */
  explicit YAMLFileReader(std::string filename);

  /**
   * @brief Parses the simulation configuration from the loaded YAML nodes.
   * @return Simulation configuration
   */
  SimulationConfig getConfig();

  // Getters for simulation parameters
  std::string getOutputBaseName() const;
  int getWriteFrequency() const;
  int getCheckpointFrequency() const;
  double getTEnd() const;
  double getDeltaT() const;
  double getEpsilon() const;
  double getSigma() const;
  double getCutoff() const;
  double getGravity() const;
  std::array<double, 3> getDomainSize() const;
  std::array<std::string, 6> getBoundaryTypesRaw() const;

  // Checkpoint-related getters
  /** @brief Returns true if this file contains checkpoint particle data */
  bool isCheckpoint() const;
  /** @brief Returns the iteration number from checkpoint (0 if not a checkpoint) */
  int getCheckpointIteration() const;
  /** @brief Returns the simulation time from checkpoint (0.0 if not a checkpoint) */
  double getCheckpointTime() const;

 private:
  /**
   * @brief Stores loaded YAML structure.
   */
  YAML::Node config;

  /**
   * @brief Path to the YAML file.
   */
  std::string filename;

  /**
   * @brief Validate configuration keys.
   */
  void checkRequiredKeys() const;
};
