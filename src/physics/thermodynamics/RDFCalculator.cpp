#include "physics/thermodynamics/RDFCalculator.h"

#include "utils/ArrayUtils.h"

#include <cmath>
#include <cstddef>  // for size_t

// Source - https://stackoverflow.com/a/49778398
// Posted by Ron, modified by community. See post 'Timeline' for change history
// Retrieved 2026-01-29, License - CC BY-SA 4.0
constexpr double pi = 3.14159265358979323846;

RDFCalculator::RDFCalculator(ParticleContainer& particles, const std::array<double, 3>& domainDims)
    : particles(particles), domainDims(domainDims) {}

const std::array<double, RDFCalculator::intervalCount>& RDFCalculator::calculateDistribution() {
  resetParticleCounts();

  // Calculate particle counts
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
      particleCounts[intervalIndex] += 2;  // add both particles
    }
  }

  // Calculate densities
  for (std::size_t i = 0; i < intervalCount; ++i) {
    // Calculate the volume of the spherical shell (difference of two spheres)
    const double r_inner = static_cast<double>(i) * sampleWidth;
    const double r_outer = r_inner + sampleWidth;
    const double shellVolume = (4.0 / 3.0) * pi * ((r_outer * r_outer * r_outer) - (r_inner * r_inner * r_inner));

    densities[i] = static_cast<double>(particleCounts[i]) / shellVolume;
  }

  return densities;
}

void RDFCalculator::resetParticleCounts() {
  for (int& particleCount : particleCounts)
    particleCount = 0;
}
