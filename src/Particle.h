/*
 * Particle.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <array>
#include <ostream>
#include <string>

class Particle {
 private:
  // Data Members
  /**
   * Position of the particle
   */
  std::array<double, 3> x;

  /**
   * Velocity of the particle
   */
  std::array<double, 3> v;

  /**
   * Force effective on this particle
   */
  std::array<double, 3> f;

  /**
   * Force which was effective on this particle
   */
  std::array<double, 3> old_f;

  /**
   * Mass of this particle
   */
  double m;

  /**
   * Type of the particle. Use it for whatever you want (e.g. to separate
   * molecules belonging to different bodies, matters, and so on)
   */
  int type;

 public:
  // Constructors
  Particle() = default;
  explicit Particle(int type_arg);
  Particle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      const std::array<double, 3>& x_arg, const std::array<double, 3>& v_arg,
      double m_arg, int type_arg = 0);

  // Rule of Five
  Particle(const Particle& other) = default;
  Particle& operator=(const Particle& other) = default;
  Particle(Particle&& other) noexcept = default;
  Particle& operator=(Particle&& other) noexcept = default;
  ~Particle() = default;

  // Getters
  [[nodiscard]] const std::array<double, 3>& getX() const;
  [[nodiscard]] const std::array<double, 3>& getV() const;
  [[nodiscard]] const std::array<double, 3>& getF() const;
  [[nodiscard]] const std::array<double, 3>& getOldF() const;
  [[nodiscard]] double getM() const;
  [[nodiscard]] int getType() const;

  // Setters
  void setF(const std::array<double, 3>& val) { this->f = val; }

  void setOldF(const std::array<double, 3>& val) { this->old_f = val; }

  void setX(const std::array<double, 3>& val) { this->x = val; }

  void setV(const std::array<double, 3>& val) { this->v = val; }

  // Operators + Utilities
  bool operator==(const Particle& other) const;
  std::string toString() const;
};

std::ostream& operator<<(std::ostream& stream, const Particle& p);
