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
  * @param particles ParticleContainer storing the particles used by the calculation method
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
 * @brief Stormer-Verlet calculation method
 */
class GravityForce final : public ForceCalc {
 public:
  using ForceCalc::ForceCalc;

  void calculateF() override;
};

class LennardJonesForce final : public ForceCalc {
private:
  const double epsilon, sigma;

public:
  LennardJonesForce(ParticleContainer& particles, const double epsilon,
                   const double sigma)
    : ForceCalc(particles), epsilon(epsilon), sigma(sigma) {}

  void calculateF() override;
};
