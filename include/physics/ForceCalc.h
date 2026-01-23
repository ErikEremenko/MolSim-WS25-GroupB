#pragma once

#include "physics/ParticleContainer.h"

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
  explicit ForceCalc(ParticleContainer& particles) : particles(particles) {}
  virtual ~ForceCalc();

  /**
  * @brief Calculates the new x-coordinate of the particle based on the Störmer-Verlet method
  * @param dt double representing the Velocity-Störmer-Verlet time step (delta t)
  */
  virtual void calculateX(double dt);
  /**
  * @brief Calculates the new velocity of the particle based on the Störmer-Verlet method
  * @param dt double representing the Velocity-Störmer-Verlet time step (delta t)
  */
  virtual void calculateV(double dt);
  /**
  * @brief Calculates the new force that acts on the particles
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
  /// Gravitational acceleration applied to all particles
  const double gravity;

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
   * @brief Lookup table for precomputed sqaured repulsion distances per type: (2^(1/6) * sigma)^2
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
   *
   * @param particles ParticleContainer that stores the particles used by the calculation method
   * @param epsilon Epsilon in the Lennard-Jones potential formula
   * @param sigma Sigma in the Lennard-Jones potential formula
   * @param cutoffRadius Distance beyond which interactions between the particles are not calculated (ignored)
   */
  LennardJonesForce(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius, double gravity);

  /**
  * @brief Calculates the Lennard-Jones forces acting on the particles
  */
  void calculateF() override;

  /**
   * @brief Calculates the Lennard-Jones forces using the direct sum O(n^2) algorithm
   */
  void calculateFDirectSum();

  /**
   * @brief Calculates the Lennard-Jones forces acting on the particles using the Linked Cell method
   */
  void calculateFLinkedCell();

  void precomputeConstants();

 private:
  /**
   * @brief Applies reflective boundary forces using ghost particles
   * @param lc Pointer to the LinkedCellParticleContainer
   */
  void applyReflectiveBoundaries(const class LinkedCellParticleContainer* lc) const;

  // TODO: Docstring these methods
  void calcFPeriodicBoundary(Particle* p1, Particle* p2) const;
  void applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const;
};

/**
 * @class LennardJonesForceParallel
 * @brief Models the Lennard-Jones potential with parallelization
 */
class LennardJonesForceParallel final : public ForceCalc {
 private:
  // TODO: Either docstring these or inherit from LennardJonesForce
  const double epsilon, sigma, cutoffRadius;

 public:
  /**
   *
   * @param particles ParticleContainer that stores the particles used by the calculation method
   * @param epsilon Epsilon in the Lennard-Jones potential formula
   * @param sigma Sigma in the Lennard-Jones potential formula
   * @param cutoffRadius Distance beyond which interactions between the particles are not calculated (ignored)
   */
  LennardJonesForceParallel(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius);

  /**
  * @brief Calculates the Lennard-Jones forces acting on the particles using OpenMP
  */
  void calculateF() override;
};
