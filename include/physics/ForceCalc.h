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
  // TODO: Consider renaming epsilon and sigma to globalEpsilon and globalSigma respectively
  const double epsilon, sigma, cutoffRadius, repulsionDistance;

 public:
  /**
   *
   * @param particles ParticleContainer that stores the particles used by the calculation method
   * @param epsilon Epsilon in the Lennard-Jones potential formula
   * @param sigma Sigma in the Lennard-Jones potential formula
   * @param cutoffRadius Distance beyond which interactions between the particles are not calculated (ignored)
   */
  LennardJonesForce(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius);

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
    const double accX;
    /** @brief Constant acceleration value in the y-direction. */
    const double accY;
    /** @brief Constant acceleration value in the z-direction. */
    const double accZ;

  public:
    /**
     * @brief Initializes acceleration values and the particle container.
     * @param particles Particle container
     * @param accX Acceleration in the x-direction
     * @param accY Acceleration in the x-direction
     * @param accZ Acceleration in the x-direction
     */
    ConstantAccelerationForce(ParticleContainer& particles, double accX, double accY, double accZ);

    /**
     * @brief Calculates the applied force according to the particle's acceleration.
     * This is uses a simple physics formula (called Newton's 2nd Law): F=ma.
     */
    void calculateF() override;
};

/**
 * @class MembraneBondForce
 * @brief Models a membrane and the bond forces that hold the molecules together.
 * @note This force should only be used when simulating a membrane.
 */
class MembraneBondForce final : public ForceCalc {
  private:
    /** @brief Stiffness constant of the bonds (k in the formula). */
    const double stiffnessConstant;

    /** @brief Average bond length of a molecule pair (r_0 in the formula). */
    const double bondLength;

  public:
    /**
     * @brief Initializes the stiffness constant and the particle container.
     * @param particles Particle container
     * @param stiffnessConstant Stiffness constant
     */
    MembraneBondForce(ParticleContainer& particles, double stiffnessConstant, double bondLength);

   /**
     * @brief Calculates the bond forces acting on every particle of the membrane.
     * @note When calculating the forces on each particle, the list of diagonal and direct neighbors are iterated
     * to calculate the total force on that single particle. Currently, this force calculation is not optimized
     * to use Newton's 3rd Law, Actio est reactio.
     */
    void calculateF() override;
};

/**
 * @class MembraneConstantForce
 * @brief Models a constant force value as specified in Task 1 of Assignment 5.
 * @note The class below could be made more generic (e.g. with parameters). For now, we see no use for it.
 */
class MembraneConstantForce final : public ForceCalc {
  private:
    static constexpr double forceEndTime = 150;  // force effective until t=150
    static constexpr int yDim = 50;  // number of particles in each row (y-dim)
    static constexpr double forceValue = 0.8;

    bool particlesInitialized;
    std::array<Particle*, 4> affectedParticles;
    double currentTime;
    double dt;

  public:
    MembraneConstantForce(ParticleContainer& particles, double dt);

   /**
     * @brief In the membrane simulation experiment, the constant force F_(Z-UP) pulls the particle
     * with x/y-indices (17/24), (17/25), (18/24) and (18/25) upwards along the z-axis. Its value is
     * 0.8 and the force is only effective until time 150.
     */
    void calculateF() override;
};
