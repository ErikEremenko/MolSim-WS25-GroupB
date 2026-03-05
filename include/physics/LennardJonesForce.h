#pragma once

#include "physics/ForceCalc.h"
#include "simulation/SimulationConfig.h"

#include <vector>

/**
 * @class LennardJonesForce
 * @brief Models the Lennard-Jones potential
 *
 * @warning This class assumes that each particle type uniquely maps to a (sigma, epsilon) pair.
 * Particles with the same type MUST have identical sigma and epsilon values.
 * Violating this assumption leads to incorrect force calculations!
 *
 * @warning For performance, zero-distance checks between particles are NOT performed in
 * calculateFLinkedCell(). Ensure particles are properly initialized with distinct positions.
 * Overlapping particles will cause division by zero and NaN/inf propagation!
 */
class LennardJonesForce final : public ForceCalc {
 private:
  /// Global Lennard-Jones epsilon parameter
  const double epsilon;
  /// Global Lennard-Jones sigma parameter
  const double sigma;
  /// Cutoff radius beyond which particle interactions are ignored
  const double cutoffRadius;
  /// Distance below which particles repel (2^(1/6) * sigma)
  const double repulsionDistance;

  /// Squared cutoff radius for optimized distance comparisons
  const double cutoffRadiusSq;

  /**
   * @brief A lookup table for every pair of particle types. Used for storing precomputed information about every pair of particle types.
   */
  std::vector<double> pairLookupTable1;
  /**
   * @brief Lookup table 2 for precomputed pair constants: 1 / (2 * sigma_ij^6)
   */
  std::vector<double> pairLookupTable2;
  /**
   * @brief Lookup table for precomputed squared repulsion distances per type: (2^(1/6) * sigma)^2
   * Used for optimized distance comparisons without sqrt
   */
  std::vector<double> repulsionDistanceSqLookup;
  /**
   * @brief Lookup table for precomputed sigma^6 per type (for reflective boundaries)
   */
  std::vector<double> sigma6Lookup;
  /**
   * @brief Lookup table for 24 * epsilon per type (precomputed multiplier for reflective boundaries)
   */
  std::vector<double> epsilon24Lookup;

  /// Width of the lookup table (max type + 1)
  int tableWidth;

 public:
  /**
   * @param particles ParticleContainer that stores the particles used by the calculation method
   * @param epsilon Epsilon in the Lennard-Jones potential formula
   * @param sigma Sigma in the Lennard-Jones potential formula
   * @param cutoffRadius Distance beyond which interactions between the particles are not calculated (ignored)
   */
  LennardJonesForce(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius);

  void calculateF() override;
  void calculateFDirectSum();
  void calculateFLinkedCell();

  /**
   * @brief Strategy 1: stride-coloring (domain decomposition)
   *
   * 2D: 3x2 coloring in X/Y (6 phases).
   * 3D: 3x2x3 coloring in X/Y/Z (18 phases).
   *
   * Uses the full 13-neighbor half-shell.
   * Phase spacing keeps write regions separate.
   *
   * Pros: no atomics, deterministic.
   * Cons: 6-18 barriers, can imbalance on inhomogeneous data.
   */
  void calculateFLinkedCellParallel1();

  /**
   * @brief Strategy 2: task-based with atomics
   *
   * One task per cell.
   * Atomics protect concurrent updates.
   *
   * Pros: fewer barriers, adapts to uneven distributions.
   * Cons: atomic + task overhead.
   */
  void calculateFLinkedCellParallel2();

  void precomputeConstants() override;

  /**
   * @brief Set the parallel strategy for force calculation
   * @param strategy COLORING or TASKBASED
   */
  void setParallelStrategy(ParallelStrategy strategy) { parallelStrategy = strategy; }

  /**
   * @brief Enable or disable parallelization
   */
  void setUseParallel(bool enable) { useParallel = enable; }

 private:
  static void applyLJPairForceGlobal(Particle& p_i, Particle& p_j, double sigma6, double epsilon, double cutoffSq);

  /**
   * @brief Compute and apply LJ pair force using lookup tables (must be defined in header for inlining)
   * 
   * This function is critical for performance ~80% of CPU time) and has to be be inlined.
   * Moving the definition to the .cpp file prevents inlining across translation units.
   */
  __attribute__((always_inline)) static inline void applyLJPairForceLookup(Particle& p_i, Particle& p_j,
                                                                           const double* __restrict__ lut1,
                                                                           const double* __restrict__ lut2,
                                                                           const int tw, const double cutoffSq) {
    const auto& xi = p_i.getX();
    const auto& xj = p_j.getX();

    const double dx = xj[0] - xi[0];
    const double dy = xj[1] - xi[1];
    const double dz = xj[2] - xi[2];
    const double distSq = dx * dx + dy * dy + dz * dz;

    if (distSq > cutoffSq) {
      return;
    }

    const double inv_distSq = 1.0 / distSq;
    const int idx = p_i.getType() * tw + p_j.getType();
    const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
    const double term = lut1[idx] * inv_distSq * inv_distSq3 * (lut2[idx] - inv_distSq3);

    const double fx = term * dx;
    const double fy = term * dy;
    const double fz = term * dz;

    auto& fi = p_i.getF();
    auto& fj = p_j.getF();
    fi[0] += fx;
    fi[1] += fy;
    fi[2] += fz;
    fj[0] -= fx;
    fj[1] -= fy;
    fj[2] -= fz;
  }

  ParallelStrategy parallelStrategy = ParallelStrategy::COLORING;
  bool useParallel = false;

  void applyReflectiveBoundaries(const class LinkedCellParticleContainer* lc) const;
  void calcFPeriodicBoundary(Particle* p1, Particle* p2) const;
  void calcFPeriodicBoundaryAtomic(Particle* p1, const std::array<double, 3>& p2_shifted_pos, Particle* p2) const;
  void applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const;
  void applyPeriodicBoundariesParallel(LinkedCellParticleContainer* lc) const;
};
