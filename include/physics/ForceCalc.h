#pragma once

#include "physics/ParticleContainer.h"
#include "simulation/SimulationConfig.h"

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
  * @brief Constructs a force calculation strategy
  * @param particles Particles affected by the force
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
  * @note This is a pure virtual method.
  */
  virtual void calculateF() = 0;

  virtual void precomputeConstants();
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
  void calculateFLinkedCellParallel1();

  void precomputeConstants() override;

 private:
  void applyReflectiveBoundaries(const class LinkedCellParticleContainer* lc) const;
  void calcFPeriodicBoundary(Particle* p1, Particle* p2) const;
  void applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const;
};

/**
 * @class TruncatedLJForce
 * @brief Repulsive-only Lennard-Jones potential, truncated at 2^(1/6)·sigma
 *
 * Used for membrane simulations to prevent self-penetration without attraction.
 */
class TruncatedLJForce final : public ForceCalc {
 public:
  TruncatedLJForce(ParticleContainer& particles);

  void calculateF() override;
};

/**
 * @class HarmonicMembraneForce
 * @brief Harmonic potential for membrane bonds between neighboring particles
 * Models springs between direct and diagonal neighbors in a 2D membrane.
 */
class HarmonicMembraneForce final : public ForceCalc {
 private:
  double stiffness;      ///< Spring constant k
  double avgBondLength;  ///< Average bond length r0

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
 *
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
