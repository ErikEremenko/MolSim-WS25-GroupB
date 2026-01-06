/**
 * @file ParticleContainer.h
 *
 *
 */

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include <array>
#include "physics/Particle.h"
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
  virtual ~ParticleContainer() = default;
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
  virtual void addParticle(std::array<double, 3> x, std::array<double, 3> v,
                   double m);  // function called in FileReader
  virtual void addParticle(std::array<double, 3> x, std::array<double, 3> v,
                   double m, double sigma, double epsilon);
  virtual void addParticle(const Particle* p);

  virtual void removeParticle(size_t idx);

  /**
   * @brief Adds a particle from checkpoint with complete state to the container
   * @param x position vector
   * @param v velocity vector
   * @param m mass
   * @param f current force vector
   * @param oldF old force vector from previous iteration
   * @param type particle type
   * @param sigma Lennard-Jones sigma parameter
   * @param epsilon Lennard-Jones epsilon parameter
   */
  virtual void addParticle(std::array<double, 3> x, std::array<double, 3> v, double m,
                           std::array<double, 3> f, std::array<double, 3> oldF,
                           int type, double sigma, double epsilon);

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
