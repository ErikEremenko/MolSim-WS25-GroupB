#include "physics/thermodynamics/RDFCalculator.h"

#include "utils/ArrayUtils.h"

#include <cstddef>  // for size_t

RDFCalculator::RDFCalculator(ParticleContainer& particles, const std::array<double, 3>& domainDims)
    : particles(particles), domainDims(domainDims) {}

const std::array<int, RDFCalculator::intervalCount>& RDFCalculator::calculateDistribution() {
  resetIntervals();

  const std::size_t n_particles = particles.size();
  for (std::size_t i = 0; i < n_particles; ++i) {
    for (std::size_t j = i + 1; j < n_particles; ++j) {  // Iterate over particle pairs
      auto& p_i = particles[i];
      auto& p_j = particles[j];

      // Calculate distance using minimum image convention for periodic boundaries
      const auto& x_i = p_i.getX();
      const auto& x_j = p_j.getX();

      std::array<double, 3> distanceVector{x_i[0] - x_j[0], x_i[1] - x_j[1], x_i[2] - x_j[2]};

      for (int dim = 0; dim < 3; dim++) {
        // Assume all boundaries are periodic
        const double halfDim = domainDims[dim] * 0.5;
        if (distanceVector[dim] > halfDim)
          distanceVector[dim] -= domainDims[dim];
        else if (distanceVector[dim] < -halfDim)
          distanceVector[dim] += domainDims[dim];
      }

      const double distance = ArrayUtils::L2Norm(distanceVector);

      // Find which interval to put this pair in
      if (distance >= maxIntervalValue)
        continue;  // pair distance exceeds our RDF range, skip

      const auto intervalIndex = static_cast<std::size_t>(distance * invSampleWidth);
      intervals[intervalIndex] += 2;  // add both particles
    }
  }

  return intervals;
}

void RDFCalculator::resetIntervals() {
  for (int& interval : intervals)
    interval = 0;
}
