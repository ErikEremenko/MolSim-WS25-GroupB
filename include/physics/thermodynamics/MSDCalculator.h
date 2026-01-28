#pragma once

#include "physics/ParticleContainer.h"

// TODO: Write docstrings for this class
class MSDCalculator {
  private:
    ParticleContainer& particles;
  public:
    explicit MSDCalculator(ParticleContainer& particles);
    double calculateDiffusion();
};
