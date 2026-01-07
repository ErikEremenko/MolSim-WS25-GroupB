#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include "io/FileReader.h"
#include "physics/ParticleContainer.h"
#include "utils/MaxwellBoltzmannDistribution.h"

class YAMLFileReader {
 public:
  /**
   * @brief Constructor loads the YAML file immediately.
   * @param filename Path to the YAML input file.
   */
  explicit YAMLFileReader(std::string filename);

  /**
   * @brief Parses cuboids section and queues orders to the particle generator.
   */
  void readFile();

  // Getters for simulation parameters
  std::string getOutputBaseName() const;
  int getWriteFrequency() const;
  int getCheckpointFrequency() const;
  double getTend() const;
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
  // Stores loaded YAML structure
  YAML::Node config;

  std::string filename;

  // used to validate configuration keys
  void checkRequiredKeys() const;
};
