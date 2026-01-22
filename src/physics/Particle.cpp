/*
 * Particle.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "physics/Particle.h"

#include <ostream>

#include "utils/ArrayUtils.h"

int Particle::idCounter = 0;

Particle::Particle(int type_arg)
    : x{0.0, 0.0, 0.0},
      v{0.0, 0.0, 0.0},
      f{0.0, 0.0, 0.0},
      old_f{0.0, 0.0, 0.0},
      m(0.0),
      type(type_arg),
      sigma(1.0),
      epsilon(5.0),
      id(idCounter++) {
  diagonalNeighbors.reserve(4);
  directNeighbors.reserve(4);
}

Particle::Particle(const std::array<double, 3>& x_arg, const std::array<double, 3>& v_arg, const double m_arg,
                   const int type_arg, const double sigma_arg, const double epsilon_arg)
    : x(x_arg),
      v(v_arg),
      f{0., 0., 0.},
      old_f{0., 0., 0.},
      m(m_arg),
      type(type_arg),
      sigma(sigma_arg),
      epsilon(epsilon_arg),
      id(idCounter++) {
  diagonalNeighbors.reserve(4);
  directNeighbors.reserve(4);
}

const std::array<double, 3>& Particle::getX() const {
  return x;
}

const std::array<double, 3>& Particle::getV() const {
  return v;
}

const std::array<double, 3>& Particle::getF() const {
  return f;
}
std::array<double, 3>& Particle::getF() {
  return f;
}

const std::array<double, 3>& Particle::getOldF() const {
  return old_f;
}

double Particle::getM() const {
  return m;
}

int Particle::getType() const {
  return type;
}

double Particle::getSigma() const {
  return sigma;
}

double Particle::getEpsilon() const {
  return epsilon;
}

int Particle::getID() const {
  return id;
}

std::string Particle::toString() const {
  std::stringstream stream;
  stream << "Particle: X:" << x << " v: " << v << " f: " << f << " old_f: " << old_f << " type: " << type;
  return stream.str();
}

void Particle::addF(const std::array<double, 3>& force) {
  f = f + force;
}

std::vector<int>& Particle::getDirectNeighbors() {
  return directNeighbors;
}

std::vector<int>& Particle::getDiagonalNeighbors() {
  return diagonalNeighbors;
}

bool Particle::operator==(const Particle& other) const {
  return (x == other.x) and (v == other.v) and (f == other.f) and (type == other.type) and (m == other.m) and
         (old_f == other.old_f);
}

std::ostream& operator<<(std::ostream& stream, const Particle& p) {
  stream << p.toString();
  return stream;
}
