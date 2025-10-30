#pragma once

#define EPSILON 1e-12  // TODO: We could move this to a 'defines.h' header

#include "ParticleContainer.h"

class CalcMethod {
 protected:
  ParticleContainer& particles;
  // TODO: Maybe move constructor and destructor to protected?
 public:
  explicit CalcMethod(ParticleContainer& particles) : particles(particles) {}
  virtual ~CalcMethod();

  virtual void calculateX(const double dt) = 0;
  virtual void calculateV(const double dt) = 0;
  virtual void calculateF(
      const double dt) = 0;  // TODO: for now we don't need dt, maybe remove?
};

class StormerVerletMethod : public CalcMethod {
 public:
  using CalcMethod::CalcMethod;

  void calculateX(const double dt) override;
  void calculateV(const double dt) override;
  void calculateF(const double dt) override;
};
