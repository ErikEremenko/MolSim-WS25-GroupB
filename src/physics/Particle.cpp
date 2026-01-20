/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "physics/Particle.h"

#include <ostream>

#include "utils/ArrayUtils.h"

Particle::Particle(int type_arg)
    : x{0.0, 0.0, 0.0},
      v{0.0, 0.0, 0.0},
      f{0.0, 0.0, 0.0},
      old_f{0.0, 0.0, 0.0},
      m(0.0),
      type(type_arg),
      sigma(1.0),
      epsilon(5.0) {}

Particle::Particle(const std::array<double, 3>& x_arg, const std::array<double, 3>& v_arg, const double m_arg,
                   const int type_arg, const double sigma_arg, const double epsilon_arg)
    : x(x_arg),
      v(v_arg),
      f{0., 0., 0.},
      old_f{0., 0., 0.},
      m(m_arg),
      type(type_arg),
      sigma(sigma_arg),
      epsilon(epsilon_arg) {}

double Particle::getSigma() const {
  return sigma;
}

double Particle::getEpsilon() const {
  return epsilon;
}

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

bool Particle::operator==(const Particle& other) const {
  return (x == other.x) and (v == other.v) and (f == other.f) and (type == other.type) and (m == other.m) and
         (old_f == other.old_f);
}

std::ostream& operator<<(std::ostream& stream, const Particle& p) {
  stream << p.toString();
  return stream;
}
