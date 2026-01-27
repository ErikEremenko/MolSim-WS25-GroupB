#pragma once

#include "physics/ParticleContainer.h"
#include "physics/thermodynamics/MSDCalculator.h"
#include "physics/thermodynamics/RDFCalculator.h"

// TODO: Write docstrings for this class
class ThermodynamicsStatistics {
  private:
    MSDCalculator msd;
    RDFCalculator rdf;
    const int updateFrequency;
  public:
    ThermodynamicsStatistics(ParticleContainer& particles, int updateFrequency);

    void updateStatistics();

   [[nodiscard]] int getUpdateFrequency() const;
};
