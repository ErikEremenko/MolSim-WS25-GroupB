/**
 * @file ParticleContainer.h
 *
 *
 */

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "Particle.h"
/**
 * @class ParticleContainer
 * @brief Iterable container class used for storing particles for a simulation.
 * It offers methods for adding, accessing particles and iterating through a set of particles
 * in an easy and efficient manner.
 *
 */
class ParticleContainer {
 private:
  /**
   * @brief A set of particles
   */
  std::vector<Particle> particles;

 public:
  ParticleContainer() = default;
  ~ParticleContainer() = default;
  /**
   * @brief Returns number of particles in the container
   */
  [[nodiscard]] std::size_t size() const;
  /**
   * @brief Adds a particle to the container
   * @param x position vector as a 3 element array
   * @param v velocity vector as a 3 element array
   * @param m mass
   */
  void addParticle(std::array<double, 3> x, std::array<double, 3> v,
                   double m);    // function called in FileReader
  void addParticle(Particle p);  // function called in FileReader

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
