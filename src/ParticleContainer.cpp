#include "ParticleContainer.h"

std::size_t ParticleContainer::size() const {
  return particles.size();
}

void ParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m) {
  particles.emplace_back(x, v, m);
}

void ParticleContainer::addParticle(const Particle* p) {
  particles.emplace_back(*p);
}

void ParticleContainer::removeParticle(const size_t idx) {
  if (idx < particles.size())
    particles.erase(particles.begin() + idx);
}

ParticleContainer::iterator ParticleContainer::begin() {
  return particles.begin();
}

ParticleContainer::iterator ParticleContainer::end() {
  return particles.end();
}

ParticleContainer::const_iterator ParticleContainer::begin() const {
  return particles.begin();
}

ParticleContainer::const_iterator ParticleContainer::end() const {
  return particles.end();
}

Particle& ParticleContainer::operator[](std::size_t i) {
  return particles[i];
}

const Particle& ParticleContainer::operator[](std::size_t i) const {
  return particles[i];
}