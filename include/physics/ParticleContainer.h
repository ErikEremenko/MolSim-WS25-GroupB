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
class ParticleContainer {  // TODO: Add missing docstrings
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
  virtual Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v,
                                double m);  // function called in FileReader
  /**
   * @brief Adds a particle with Lennard-Jones parameters to the container.
   * @param x position vector as a 3 element array
   * @param v velocity vector as a 3 element array
   * @param m mass
   * @param sigma Lennard-Jones sigma parameter
   * @param epsilon Lennard-Jones epsilon parameter
   */
  virtual Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, double sigma,
                                double epsilon);
  /**
   * @brief Adds a typed particle with Lennard-Jones parameters to the container.
   * @param x position vector as a 3 element array
   * @param v velocity vector as a 3 element array
   * @param m mass
   * @param type particle type
   * @param sigma Lennard-Jones sigma parameter
   * @param epsilon Lennard-Jones epsilon parameter
   */
  virtual Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, int type, double sigma,
                                double epsilon);
  /**
   * @brief Adds an existing particle to the container.
   * @param p pointer to particle to be added
   */
  virtual Particle* addParticle(const Particle* p);

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
  virtual Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, std::array<double, 3> f,
                                std::array<double, 3> oldF, int type, double sigma, double epsilon);

  virtual void removeParticle(size_t idx);

  using iterator = std::vector<Particle>::iterator;
  using const_iterator = std::vector<Particle>::const_iterator;

  /**
   * @brief Returns an iterator to the first particle.
   */
  iterator begin();

  /**
   * @brief Returns an iterator to one past the last particle.
   */
  iterator end();

  /**
   * @brief Returns a const iterator to the first particle.
   */
  [[nodiscard]] const_iterator begin() const;

  /**
   * @brief Returns a const iterator to one past the last particle.
   */
  [[nodiscard]] const_iterator end() const;

  /**
   * @brief Returns a mutable reference to the particle at index @p i.
   * @param i particle index
   * @return reference to the particle at index @p i
   */
  Particle& operator[](std::size_t i);

  /**
   * @brief Returns a const reference to the particle at index @p i.
   * @param i particle index
   * @return const reference to the particle at index @p i
   */
  const Particle& operator[](std::size_t i) const;
};
