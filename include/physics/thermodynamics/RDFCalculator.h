#pragma once

#include "physics/ParticleContainer.h"

// TODO: Write docstrings for this class
class RDFCalculator {
  private:
    static constexpr double sampleWidth = 1.0;  // this is delta*r in the slide
    static constexpr double invSampleWidth = 1.0 / sampleWidth;  // pre-compute to avoid division in O(N^2) loop
    static constexpr int intervalCount = 50;
    static constexpr double maxIntervalValue = intervalCount*sampleWidth;

    ParticleContainer& particles;
    std::array<double, 3> domainDims;

    std::array<int, intervalCount> intervals{};

    void resetIntervals();
  public:
    explicit RDFCalculator(ParticleContainer& particles, const std::array<double, 3>& domainDims);
    void calculateDistribution();

    [[nodiscard]] std::array<int, intervalCount>& getIntervals();
};
