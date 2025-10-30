#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "Particle.h"

class ParticleContainer {
 private:
  std::vector<Particle> particles;

 public:
  ParticleContainer() = default;
  ~ParticleContainer() = default;

  [[nodiscard]] std::size_t size() const;
  void addParticle(std::array<double, 3> x, std::array<double, 3> v,
                   double m);  // function called in FileReader

  using iterator = std::vector<Particle>::iterator;
  using const_iterator = std::vector<Particle>::const_iterator;

  // Iteration over single particles
  iterator begin();
  iterator end();

  [[nodiscard]] const_iterator begin() const;
  [[nodiscard]] const_iterator end() const;

  Particle& operator[](std::size_t i);
  const Particle& operator[](std::size_t i) const;
};
