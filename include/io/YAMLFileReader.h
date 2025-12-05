#pragma once
#include <yaml-cpp/yaml.h>
#include <string>
#include "ParticleContainer.h"
#include "io/FileReader.h"
#include "utils/MaxwellBoltzmannDistribution.h"

class YAMLFileReader final : public BaseFileReader {
 public:
  /**
   * @brief Constructor loads the YAML file immediately.
   * @param filename Path to the YAML input file.
   */
  explicit YAMLFileReader(std::string filename);

  /**
   * @brief Parses cuboids section and populates the container.
   * @param particles The container to add particles to.
   */
  void readFile(ParticleContainer& particles) override;

  // Getters for simulation parameters
  std::string getOutputBaseName() const;
  int getWriteFrequency() const;
  double getTend() const;
  double getDeltaT() const;
  double getEpsilon() const;
  double getSigma() const;
  double getCutoff() const;
  std::array<double,3> getDomainSize() const;
  std::array<std::string,6> getBoundaryTypesRaw() const;
  double getLJRepulsionDistance() const;

 private:
  YAML::Node config;  // stores loaded YAML structure

  // used to validate configuration keys
  void checkRequiredKeys() const;
};
