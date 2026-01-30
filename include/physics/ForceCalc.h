#pragma once

#include "physics/LinkedCellParticleContainer.h"
#include "physics/ParticleContainer.h"
#include "simulation/SimulationConfig.h"

#include <cmath>
#include <utility>
#include <vector>

/**
 * @class ForceCalc
 * @brief Abstract base class used for implementing different force calculation strategies
 *
 * @section perf_requirements Performance Requirements
 * For optimal performance with millions/billions of particle operations:
 *
 * @warning **Zero-Distance Checks**: Force calculation methods in derived classes may omit
 * zero-distance checks for performance. Ensure particles are initialized with distinct positions.
 * Overlapping particles will cause division by zero, resulting in NaN/inf values that silently
 * propagate through the simulation.
 *
 * @warning **Particle Types**: The LennardJonesForce class requires that particle types
 * uniquely identify (sigma, epsilon) pairs. Use the automatic type assignment in YAMLFileReader
 * or ensure manual type assignments are consistent.
 */
class ForceCalc {
 protected:
  /**
   * @brief Particles to apply the force calculations on
   */
  ParticleContainer& particles;

 public:
  /**
  * @brief Constructor
  * @param particles ParticleContainer that stores the particles used by the calculation method
  */
  explicit ForceCalc(ParticleContainer& particles);
  virtual ~ForceCalc();

  /**
  * @brief Calculates the new position of the particle based on the Störmer-Verlet method
  * @param particles Particles to move
  * @param dt Velocity-Störmer-Verlet time step (delta t)
  */
  static void calculateX(ParticleContainer& particles, double dt);
  /**
  * @brief Calculates the new velocity of the particle based on the Störmer-Verlet method
  * @param particles Particles whose velocities are to be calculated
  * @param dt Velocity-Störmer-Verlet time step (delta t)
  */
  static void calculateV(ParticleContainer& particles, double dt);
  /**
  * @brief Calculates the new force that acts on the particles
  */
  virtual void calculateF() = 0;

  virtual void precomputeConstants();

 protected:
  /**
   * @brief Helper template to iterate over periodic boundary pairs and apply a force function
   * 
   * This eliminates code duplication between LennardJonesForce and SmoothedLJForce.
   * The ForceFunc callable should have signature: void(Particle* p1, Particle* p2)
   * 
   * @tparam ForceFunc Callable type for computing forces between particle pairs
   * @param lc LinkedCellParticleContainer to iterate over
   * @param calcForce Function to calculate and apply force between two particles
   */
  template <typename ForceFunc>
  static void applyPeriodicBoundariesImpl(class LinkedCellParticleContainer* lc, ForceFunc&& calcForce);
};

/**
 * @class GravityForce
 * @brief Models gravity forces between particles
 */
class GravityForce final : public ForceCalc {
 public:
  using ForceCalc::ForceCalc;

  /**
  * @brief Calculates the gravity forces acting on the particles
  */
  void calculateF() override;
};

/**
 * @class GlobalGravityForce
 * @brief Applies a constant gravitational acceleration to all particles along a specified axis
 */
class GlobalGravityForce final : public ForceCalc {
 private:
  double gravity;  ///< Gravitational acceleration magnitude
  int axis;        ///< Axis: 0=x, 1=y, 2=z

 public:
  /**
   * @param particles ParticleContainer that stores the particles
   * @param g Gravity acceleration value (e.g., -9.81 for downward gravity)
   * @param axis Axis direction: 0=x, 1=y (default), 2=z
   */
  GlobalGravityForce(ParticleContainer& particles, double g, int axis = 1);

  void calculateF() override;
};

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
   * @brief Strategy 1: C-coloring (domain decomposition)
   *
    * 2D: 3×2 coloring in X/Y (6 phases).
    * 3D: 2×2×2 coloring (8 phases).
   *
  * Pros: no atomics, cache-friendly, deterministic.
  * Cons: 6–8 barriers, can imbalance on inhomogeneous data.
   */
  void calculateFLinkedCellParallel1();

  /**
   * @brief Strategy 2: task-based with atomics
   *
   * One task per cell, work-stealing balances load.
   * Atomics protect inter-cell updates.
   *
   * Pros: fewer barriers, adapts to inhomogeneous data.
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

/**
 * @class SmoothedLJForce
 * @brief Smoothed Lennard-Jones potential with continuous force at cutoff
 *
 * Uses the standard Lennard-Jones force up to the smoothing radius, then
 * tapers to zero at the cutoff radius. This avoids discontinuities
 * in the force. Mixed particle types use Lorentz-Berthelot rules.
 */
class SmoothedLJForce final : public ForceCalc {
 private:
  /// Global Lennard-Jones epsilon parameter
  const double epsilon;
  /// Global Lennard-Jones sigma parameter
  const double sigma;
  /// Cutoff radius (r_c) beyond which interactions are zero
  const double cutoffRadius;
  /// Smoothing start radius (r_l) below which standard LJ is used
  const double smoothingRadius;

  // Precomputed constants
  double cutoffRadiusSq;     ///< r_c^2 for squared distance comparisons
  double smoothingRadiusSq;  ///< r_l^2 for squared distance comparisons
  double rcMinusRlCubed;     ///< (r_c - r_l)^3 for smoothing function

  /// Lookup tables for pair interactions (indexed by type_i * tableWidth + type_j)
  std::vector<double> pairEpsilon;  ///< Mixed epsilon values
  std::vector<double> pairSigma6;   ///< sigma_ij^6
  std::vector<double> pairSigma12;  ///< sigma_ij^12
  int tableWidth = 0;

 public:
  /**
   * @param particles ParticleContainer that stores the particles
   * @param epsilon Epsilon in the Lennard-Jones potential formula
   * @param sigma Sigma in the Lennard-Jones potential formula
   * @param cutoffRadius Cutoff radius r_c beyond which interactions are ignored
   * @param smoothingRadius Smoothing radius r_l below which standard LJ is used
   */
  SmoothedLJForce(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius,
                  double smoothingRadius);

  void calculateF() override;
  void precomputeConstants() override;

 private:
  void calculateFDirectSum() const;
  void calculateFLinkedCell();
  void calcFPeriodicPair(Particle* p1, Particle* p2) const;
  void applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const;

  /**
   * @brief Compute and apply smoothed LJ pair force (must be defined in header for proper inlining)
   */
  __attribute__((always_inline)) inline void applySmoothedPairForce(Particle& p_i, Particle& p_j) const {
    const auto& xi = p_i.getX();
    const auto& xj = p_j.getX();

    const double dx = xj[0] - xi[0];
    const double dy = xj[1] - xi[1];
    const double dz = xj[2] - xi[2];
    const double distSq = dx * dx + dy * dy + dz * dz;

    if (distSq >= cutoffRadiusSq) {
      return;
    }

    // Use precomputed lookup tables
    const int idx = p_i.getType() * tableWidth + p_j.getType();
    const double epsilon_ij = pairEpsilon[idx];
    const double sigma6 = pairSigma6[idx];
    const double sigma12 = pairSigma12[idx];

    const double inv_distSq = 1.0 / distSq;
    const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
    const double sigma6_d6 = sigma6 * inv_distSq3;
    const double sigma12_d12 = sigma12 * inv_distSq3 * inv_distSq3;

    const double U_LJ = 4.0 * epsilon_ij * (sigma12_d12 - sigma6_d6);
    const double F_LJ_scalar = 24.0 * epsilon_ij * inv_distSq * (sigma6_d6 - 2.0 * sigma12_d12);

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    if (distSq <= smoothingRadiusSq) {
      fx = F_LJ_scalar * dx;
      fy = F_LJ_scalar * dy;
      fz = F_LJ_scalar * dz;
    } else {
      const double d = std::sqrt(distSq);
      const double dMinusRl = d - smoothingRadius;
      const double dMinusRl2 = dMinusRl * dMinusRl;
      const double S = 1.0 - dMinusRl2 * (3.0 * cutoffRadius - smoothingRadius - 2.0 * d) / rcMinusRlCubed;
      const double dS_dd = -6.0 * dMinusRl * (cutoffRadius - d) / rcMinusRlCubed;

      const double F_total_scalar = S * F_LJ_scalar + U_LJ * dS_dd / d;
      fx = F_total_scalar * dx;
      fy = F_total_scalar * dy;
      fz = F_total_scalar * dz;
    }

    auto& fi = p_i.getF();
    auto& fj = p_j.getF();
    fi[0] += fx;
    fi[1] += fy;
    fi[2] += fz;
    fj[0] -= fx;
    fj[1] -= fy;
    fj[2] -= fz;
  }
};

/**
 * @class TruncatedLJForce
 * @brief Repulsive-only Lennard-Jones potential, truncated at 2^(1/6)·sigma
 * Used for membrane simulations to prevent self-penetration without attraction.
 * Uses Lorentz-Berthelot mixing rules for mixed particle types.
 */
class TruncatedLJForce final : public ForceCalc {
 public:
  explicit TruncatedLJForce(ParticleContainer& particles);

  void calculateF() override;
};

/**
 * @class HarmonicMembraneForce
 * @brief Harmonic potential for membrane bonds between neighboring particles
 * Models interactions between direct and diagonal neighbors in a 2D membrane.
 */
class HarmonicMembraneForce final : public ForceCalc {
 private:
  double stiffness;           ///< Stiffness constant k
  double avgBondLength;       ///< Average bond length r0 for direct neighbors
  double diagonalBondLength;  ///< Bond length for diagonal neighbors: sqrt(2) * r0

 public:
  /**
   * @param particles ParticleContainer with membrane particles (must have neighbor info set)
   * @param k Stiffness constant
   * @param r0 Average bond length for direct neighbors
   */
  HarmonicMembraneForce(ParticleContainer& particles, double k, double r0);

  void calculateF() override;
};

/**
 * @class ConstantForce
 * @brief Applies a constant force to specific particles (identified by membrane x/y indices)
 * The force is only applied until a specified end time ("pulling" membrane particles)
 */
class ConstantForce final : public ForceCalc {
 private:
  std::array<double, 3> force;
  double endTime;
  double& currentTime;                             ///< Reference to simulation's current time
  std::vector<std::pair<int, int>> targetIndices;  ///< x/y indices of target particles
  int membraneDimY;                                ///< Y-dimension of membrane grid for index calculation

 public:
  /**
   * @param particles ParticleContainer
   * @param fx Force in x direction
   * @param fy Force in y direction
   * @param fz Force in z direction
   * @param endTime Time after which force stops being applied
   * @param currentTime Reference to the simulation's current time variable
   * @param targetIndices Vector of (x, y) index pairs identifying which particles to pull
   * @param membraneDimY Y-dimension of the membrane for calculating particle indices
   */
  ConstantForce(ParticleContainer& particles, double fx, double fy, double fz, double endTime, double& currentTime,
                std::vector<std::pair<int, int>> targetIndices, int membraneDimY);

  void calculateF() override;
};

/**
 * @class LennardJonesForceParallel
 * @brief Models the Lennard-Jones potential with parallelization
 */
class LennardJonesForceParallel final : public ForceCalc {
 private:
  const double epsilon, sigma, cutoffRadius;

 public:
  LennardJonesForceParallel(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius);
  void calculateF() override;
};

// Template implementation (must be in header)
template <typename ForceFunc>
void ForceCalc::applyPeriodicBoundariesImpl(LinkedCellParticleContainer* lc, ForceFunc&& calcForce) {
  const auto domainDims = lc->domain_dims();
  const auto boundaryTypes = lc->boundary_types();
  const auto numCells = lc->num_cells();

  // Face interactions
  // X-periodic: left wall (x=1) <-> right wall (x=numCells[0]-2)
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC) {
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++)
      for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
        auto& cell1 = lc->cell_at(1, c1y, c1z);
        for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++)
          for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
            if (c2y < 1 || c2y > numCells[1] - 2 || c2z < 1 || c2z > numCells[2] - 2)
              continue;
            auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, c2z);
            for (auto& p1 : cell1)
              for (auto& p2 : cell2) {
                p2->setX(p2->getX()[0] - domainDims[0], 0);
                calcForce(p1, p2);
                p2->setX(p2->getX()[0] + domainDims[0], 0);
              }
          }
      }
  }

  // Y-periodic: bottom wall (y=1) <-> top wall (y=numCells[1]-2)
  if (boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++)
      for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
        auto& cell1 = lc->cell_at(c1x, 1, c1z);
        for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++)
          for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
            if (c2x < 1 || c2x > numCells[0] - 2 || c2z < 1 || c2z > numCells[2] - 2)
              continue;
            auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, c2z);
            for (auto& p1 : cell1)
              for (auto& p2 : cell2) {
                p2->setX(p2->getX()[1] - domainDims[1], 1);
                calcForce(p1, p2);
                p2->setX(p2->getX()[1] + domainDims[1], 1);
              }
          }
      }
  }

  // Z-periodic: front wall (z=1) <-> back wall (z=numCells[2]-2)
  if (boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++)
      for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
        auto& cell1 = lc->cell_at(c1x, c1y, 1);
        for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++)
          for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
            if (c2x < 1 || c2x > numCells[0] - 2 || c2y < 1 || c2y > numCells[1] - 2)
              continue;
            auto& cell2 = lc->cell_at(c2x, c2y, numCells[2] - 2);
            for (auto& p1 : cell1)
              for (auto& p2 : cell2) {
                p2->setX(p2->getX()[2] - domainDims[2], 2);
                calcForce(p1, p2);
                p2->setX(p2->getX()[2] + domainDims[2], 2);
              }
          }
      }
  }

  // Edge interactions (two periodic dimensions)
  // X+Y periodic edges
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC) {
    for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
      auto& cell1 = lc->cell_at(1, 1, c1z);
      for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
        if (c2z < 1 || c2z > numCells[2] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, c2z);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[1] - domainDims[1], 1);
            calcForce(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
          }
      }
    }
    for (int c1z = 1; c1z < numCells[2] - 1; c1z++) {
      auto& cell1 = lc->cell_at(1, numCells[1] - 2, c1z);
      for (int c2z = c1z - 1; c2z <= c1z + 1; c2z++) {
        if (c2z < 1 || c2z > numCells[2] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, 1, c2z);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
            calcForce(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[1] - domainDims[1], 1);
          }
      }
    }
  }

  // X+Z periodic edges
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      auto& cell1 = lc->cell_at(1, c1y, 1);
      for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
        if (c2y < 1 || c2y > numCells[1] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, numCells[2] - 2);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
            calcForce(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
          }
      }
    }
    for (int c1y = 1; c1y < numCells[1] - 1; c1y++) {
      auto& cell1 = lc->cell_at(1, c1y, numCells[2] - 2);
      for (int c2y = c1y - 1; c2y <= c1y + 1; c2y++) {
        if (c2y < 1 || c2y > numCells[1] - 2)
          continue;
        auto& cell2 = lc->cell_at(numCells[0] - 2, c2y, 1);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[0] - domainDims[0], 0);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
            calcForce(p1, p2);
            p2->setX(p2->getX()[0] + domainDims[0], 0);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
          }
      }
    }
  }

  // Y+Z periodic edges
  if (boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      auto& cell1 = lc->cell_at(c1x, 1, 1);
      for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {
        if (c2x < 1 || c2x > numCells[0] - 2)
          continue;
        auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, numCells[2] - 2);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[1] - domainDims[1], 1);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
            calcForce(p1, p2);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
          }
      }
    }
    for (int c1x = 1; c1x < numCells[0] - 1; c1x++) {
      auto& cell1 = lc->cell_at(c1x, 1, numCells[2] - 2);
      for (int c2x = c1x - 1; c2x <= c1x + 1; c2x++) {
        if (c2x < 1 || c2x > numCells[0] - 2)
          continue;
        auto& cell2 = lc->cell_at(c2x, numCells[1] - 2, 1);
        for (auto& p1 : cell1)
          for (auto& p2 : cell2) {
            p2->setX(p2->getX()[1] - domainDims[1], 1);
            p2->setX(p2->getX()[2] + domainDims[2], 2);
            calcForce(p1, p2);
            p2->setX(p2->getX()[1] + domainDims[1], 1);
            p2->setX(p2->getX()[2] - domainDims[2], 2);
          }
      }
    }
  }

  // Corner interactions (all 8 corners in fully periodic 3D domain)
  if (boundaryTypes[0] == BoundaryType::PERIODIC && boundaryTypes[1] == BoundaryType::PERIODIC &&
      boundaryTypes[2] == BoundaryType::PERIODIC && boundaryTypes[3] == BoundaryType::PERIODIC &&
      boundaryTypes[4] == BoundaryType::PERIODIC && boundaryTypes[5] == BoundaryType::PERIODIC) {

    auto& cell1 = lc->cell_at(1, 1, 1);
    auto& cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, numCells[2] - 2);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
        calcForce(p1, p2);
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
      }

    cell1 = lc->cell_at(numCells[0] - 2, 1, 1);
    cell2 = lc->cell_at(1, numCells[1] - 2, numCells[2] - 2);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
        calcForce(p1, p2);
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
      }

    cell1 = lc->cell_at(1, 1, numCells[2] - 2);
    cell2 = lc->cell_at(numCells[0] - 2, numCells[1] - 2, 1);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
        calcForce(p1, p2);
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
      }

    cell1 = lc->cell_at(numCells[0] - 2, 1, numCells[2] - 2);
    cell2 = lc->cell_at(1, numCells[1] - 2, 1);
    for (auto& p1 : cell1)
      for (auto& p2 : cell2) {
        p2->setX(p2->getX()[0] + domainDims[0], 0);
        p2->setX(p2->getX()[1] - domainDims[1], 1);
        p2->setX(p2->getX()[2] + domainDims[2], 2);
        calcForce(p1, p2);
        p2->setX(p2->getX()[0] - domainDims[0], 0);
        p2->setX(p2->getX()[1] + domainDims[1], 1);
        p2->setX(p2->getX()[2] - domainDims[2], 2);
      }
  }
}
