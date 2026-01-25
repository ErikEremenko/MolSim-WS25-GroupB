#pragma once

#include "physics/ParticleContainer.h"
#include "simulation/SimulationConfig.h"

#include <utility>
#include <vector>

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
 */
class LennardJonesForce final : public ForceCalc {
 private:
  const double epsilon, sigma, cutoffRadius, repulsionDistance;

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
  double stiffness;     ///< Spring constant k
  double avgBondLength; ///< Average bond length r0

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
  double& currentTime;  ///< Reference to simulation's current time
  std::vector<std::pair<int, int>> targetIndices;  ///< x/y indices of target particles
  int membraneDimY;  ///< Y-dimension of membrane grid for index calculation

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
  ConstantForce(ParticleContainer& particles, double fx, double fy, double fz, 
                double endTime, double& currentTime, 
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
