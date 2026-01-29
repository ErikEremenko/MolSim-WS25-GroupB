#pragma once

#include <array>

#include "Thermostat.h"
#include "io/ThermodynamicsWriter.h"
#include "physics/ParticleContainer.h"
#include "physics/thermodynamics/MSDCalculator.h"
#include "physics/thermodynamics/RDFCalculator.h"

// TODO: Write docstrings for this class
class ThermodynamicsStatistics {
 private:
  MSDCalculator msd;
  RDFCalculator rdf;
  Thermostat* thermostat;

  ThermodynamicsWriter writer;
  int writeIteration = 0;

 public:
  static constexpr int updateFrequency = 1000;

  ThermodynamicsStatistics(ParticleContainer& particles, Thermostat* thermostat,
    const std::array<double, 3>& domainDims, const std::string& baseName);

  void updateStatistics();
};
