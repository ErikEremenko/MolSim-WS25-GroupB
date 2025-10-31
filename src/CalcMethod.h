/**
 * @file CalcMethod.h
 *
 */

#pragma once

#define EPSILON 1e-12  // TODO: We could move this to a 'defines.h' header

#include "ParticleContainer.h"
/**
 * @class CalcMethod
 * @brief Virtual class used as a base for different calculation methods for simulation
 *
 */
class CalcMethod {
 protected:
  ParticleContainer& particles;
  // TODO: Maybe move constructor and destructor to protected?
 public:
/**
 * @brief Constructor
 * @param particles ParticleContainer storing the particles used by the calculation method
 */
  explicit CalcMethod(ParticleContainer& particles) : particles(particles) {}
  virtual ~CalcMethod();

  virtual void calculateX(const double dt) = 0;
  virtual void calculateV(const double dt) = 0;
  virtual void calculateF(
      const double dt) = 0;  // TODO: for now we don't need dt, maybe remove?
};

/**
 * @brief Stormer-Verlet calculation method
 */
class StormerVerletMethod : public CalcMethod {
 public:
  using CalcMethod::CalcMethod;

  void calculateX(const double dt) override;
  void calculateV(const double dt) override;
  void calculateF(const double dt) override;
};
