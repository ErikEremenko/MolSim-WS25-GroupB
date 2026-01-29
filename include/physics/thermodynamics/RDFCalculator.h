#pragma once

#include "physics/ParticleContainer.h"

// TODO: Write docstrings for this class
class RDFCalculator {
 public:
  static constexpr int intervalCount = 50;  // define here due to other constexpr's in the private section
 private:
  static constexpr double sampleWidth = 1.0;                   // this is delta*r in the slide
  static constexpr double invSampleWidth = 1.0 / sampleWidth;  // pre-compute to avoid division in O(N^2) loop

  static constexpr double maxIntervalValue = intervalCount * sampleWidth;

  ParticleContainer& particles;
  std::array<double, 3> domainDims;

  std::array<int, intervalCount> particleCounts{};
  std::array<double, intervalCount> densities{};

  void resetParticleCounts();
 public:
  explicit RDFCalculator(ParticleContainer& particles, const std::array<double, 3>& domainDims);
  const std::array<double, intervalCount>& calculateDistribution();
};
