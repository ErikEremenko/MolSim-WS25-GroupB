/**
 * @file CalcMethod.h
 *
 */

#pragma once

#include "ParticleContainer.h"
/**
 * @class ForceCalc
 * @brief Virtual class used as a base for different calculation methods for simulation
 *
 */
class ForceCalc {
 protected:
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
  const double epsilon, sigma, cutoffRadius;

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
* @brief Calculates the Lennard-Jones forces acting on the particles using the Linked Cell method
*/
  void calculateFLinkedCell(ParticleContainer& particles) const;
};

/**
 * @class LennardJonesForceParallel
 * @brief Models the Lennard-Jones potential with parallelization
 */
class LennardJonesForceParallel final : public ForceCalc {
private:
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

