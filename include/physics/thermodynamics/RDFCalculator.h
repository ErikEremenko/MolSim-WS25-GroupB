#pragma once

#include "physics/ParticleContainer.h"

// TODO: Write docstrings for this class
class RDFCalculator {
  private:
    ParticleContainer& particles;
  public:
    explicit RDFCalculator(ParticleContainer& particles);
    void calculateDistribution();
};
