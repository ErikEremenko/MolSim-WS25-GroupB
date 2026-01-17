#include "physics/ParticleContainer.h"

std::size_t ParticleContainer::size() const {
  return particles.size();
}

Particle* ParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m) {
  return &particles.emplace_back(x, v, m);
}

Particle* ParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, double sigma,
                                    double epsilon) {
  return &particles.emplace_back(x, v, m, 0, sigma, epsilon);
}

Particle* ParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, int type, double sigma,
                                    double epsilon) {
  return &particles.emplace_back(x, v, m, type, sigma, epsilon);
}

Particle* ParticleContainer::addParticle(const Particle* p) {
  return &particles.emplace_back(*p);
}

Particle* ParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, std::array<double, 3> f,
                                    std::array<double, 3> oldF, int type, double sigma, double epsilon) {
  // TODO: This code can be made simpler by just adding a constructor that oldF to be f in an initializer list
  Particle p(x, v, m, type, sigma, epsilon);
  p.setF(f);
  p.setOldF(oldF);
  return &particles.emplace_back(std::move(p));
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