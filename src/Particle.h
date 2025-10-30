/**
 * @file Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <ostream>
#include <string>

/**
 *
 *  @class Particle
 *  @brief A class representing a particle in 3 dimensional space with position, velocity,
 *  mass, type and force.
 *
 *  A user may use this class to perform a particle simulation in discrete time steps by
 *  updating the particles' position and speed based on their relative positions.
 */
class Particle {
 private:
  /** @name Data members */
  ///@{
  /**
   * @brief Position of the particle
   *
   *  A 3 component vector storing the particle's position
   */
  std::array<double, 3> x;

  /**
   * @brief Current velocity of the particle
   *
   *  A 3 component vector storing the particle's movement direction and speed
   */
  std::array<double, 3> v;

  /**
   * @brief Force effective on the particle
   *
   *  A 3 component vector storing the force vector of the force effective on the particle
   */
  std::array<double, 3> f;

  /**
   * @brief Force vector of the force that was effective on the particle in the previous iteration
   */
  std::array<double, 3> old_f;

  /**
   * @brief Mass of the particle
   *
   */
  double m;

  /**
   * @brief Type of the particle
   *
   */
  int type;
  ///@}
 public:


  /**@name Constructors */
  ///@{
  Particle() = default;
  /**
   * @brief Constructor which specifies the type of the particle
   * @param type_arg integer value that offers the ability to differentiate particles
   */
  explicit Particle(int type_arg);
  /**
   * @brief Constructor which specifies position, velocity, mass and type of the particle
   * @param x_arg position of particle
   * @param v_arg velocity vector of particle
   * @param m_arg mass of particle
   * @param type_arg integer value that offers the ability to differentiate particles
   */
  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second

      const std::array<double, 3>& x_arg, const std::array<double, 3>& v_arg
      , double m_arg, int type_arg = 0);
  ///@}
  /**
   * @brief Rule of Five holds
   */
  Particle(const Particle& other) = default;
  Particle &operator=(const Particle& other) = default;
  Particle(Particle&& other) noexcept = default;
  Particle &operator=(Particle&& other) noexcept = default;
  ~Particle() = default;

  /** @name Getter methods */
  ///@{
  /** @brief get position of particle */
  [[nodiscard]] const std::array<double, 3>& getX() const;
  /** @brief get velocity of particle */
  [[nodiscard]] const std::array<double, 3>& getV() const;
  /** @brief get force effective on particle */
  [[nodiscard]] const std::array<double, 3>& getF() const;
  /** @brief get old force effective on particle */
  [[nodiscard]] const std::array<double, 3>& getOldF() const;
  /** @brief get mass of particle */

  [[nodiscard]] double getM() const;
  /** @brief get type of particle */
  [[nodiscard]] int getType() const;
  ///@}

  /** @name Setter methods */
  ///@{
  /** @brief set particle position vector
   *  @param val velocity vector as 3 element array
   */
  void setX(const std::array<double, 3> &val) { this->x = val; }
  /** @brief set particle velocity vector
   *  @param val velocity vector as 3 element array
   */
  void setV(const std::array<double, 3> &val) { this->v = val; }
  /** @brief set force effective on particle
   *  @param val force vector as 3 element array
   */
  void setF(const std::array<double, 3> &val) { this->f = val; }
  /** @brief set old force effective on particle
   *  @param val force vector as 3 element array
   */
  void setOldF(const std::array<double, 3> &val) { this->old_f = val; }
  ///@}

  /** @name Operators + Utilities */
  bool operator==(const Particle &other) const;

  std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, const Particle& p);
