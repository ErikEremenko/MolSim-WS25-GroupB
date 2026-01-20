#pragma once

#include "physics/ParticleContainer.h"

/**
 * @class ForceCalc
 * @brief Abstract base class used for implementing different force calculation strategies
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
 */
class LennardJonesForce final : public ForceCalc {
 private:
  // TODO: Docstring these members
  const double epsilon, sigma, cutoffRadius, repulsionDistance, gravity;

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

/**
 * @class ConstantAccelerationForce
 * @brief Models a force that provides constant acceleration for every particle.
 *
 * This force can be used to model earth's gravity.
 */
class ConstantAccelerationForce final : public ForceCalc {
  private:
    /** @brief Constant acceleration value in the x-direction. */
    double accX;
    /** @brief Constant acceleration value in the y-direction. */
    double accY;
    /** @brief Constant acceleration value in the z-direction. */
    double accZ;

  public:
    /**
     * @brief Initializes acceleration values and the particle container.
     * @param particles
     * @param accX
     * @param accY
     * @param accZ
     */
    ConstantAccelerationForce(ParticleContainer& particles, double accX, double accY, double accZ);

    /**
     * @brief Calculates the applied force according to the particle's acceleration.
     * This is uses a simple physics formula (called Newton's 2nd Law): F=ma.
     */
    void calculateF() override;
};
