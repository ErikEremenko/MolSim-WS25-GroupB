/**
 * @file CalcMethod.h
 *
 */

#pragma once

#include "ParticleContainer.h"
/**
 * @class CalcMethod
 * @brief Virtual class used as a base for different calculation methods for simulation
 *
 */
class CalcMethod {
 protected:
  ParticleContainer& particles;

 public:
  /**
 * @brief Constructor
 * @param particles ParticleContainer storing the particles used by the calculation method
 */
  explicit CalcMethod(ParticleContainer& particles) : particles(particles) {}
  virtual ~CalcMethod();

  /**
 * @brief Calculates the new x-coordinate of the particle
 * @param dt double representing the Velocity-Störmer-Verlet time step (delta t)
 */
  virtual void calculateX(double dt) = 0;
  /**
* @brief Calculates the new velocity of the particle
* @param dt double representing the Velocity-Störmer-Verlet time step (delta t)
*/
  virtual void calculateV(double dt) = 0;
  /**
* @brief Calculates the new force that acts on the particle
*/
  virtual void calculateGravityF() = 0;
};

/**
 * @brief Stormer-Verlet calculation method
 */
class StormerVerletMethod final : public CalcMethod {
 public:
  using CalcMethod::CalcMethod;

  void calculateX(double dt) override;
  void calculateV(double dt) override;
  void calculateGravityF() override;
  void calculateLennardJonesF(double epsilon, double sigma) const;
};
