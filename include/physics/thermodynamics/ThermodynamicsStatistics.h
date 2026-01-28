#pragma once

#include <array>

#include "physics/ParticleContainer.h"
#include "physics/thermodynamics/MSDCalculator.h"
#include "physics/thermodynamics/RDFCalculator.h"

// TODO: Write docstrings for this class
class ThermodynamicsStatistics {
  private:
    MSDCalculator msd;
    RDFCalculator rdf;
  public:
    static constexpr int updateFrequency = 1000;

    ThermodynamicsStatistics(ParticleContainer& particles, const std::array<double, 3>& domainDims);

    void updateStatistics();
};
