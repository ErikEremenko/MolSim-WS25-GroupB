
#include "../include/ParticleGenerator.h"
#include "../include/utils/MaxwellBoltzmannDistribution.h"

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#endif  // SPDLOG_ACTIVE_LEVEL
#include <spdlog/spdlog.h>

ParticleGenerator::ParticleGenerator(ParticleContainer& particles) : particles(particles) {}
ParticleGenerator::~ParticleGenerator() = default;

void ParticleGenerator::generateCuboid(std::array<double, 3> cx, std::array<double, 3> cv, std::array<int, 3> n,
                                       double h, double m, double t) const {
  for (int nx = 0; nx < n[0]; nx++) {
    for (int ny = 0; ny < n[1]; ny++) {
      for (int nz = 0; nz < n[2]; nz++) {
        std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(t, m, 2);
        std::array<double, 3> particleVelocity = {cv[0] + temperatureVel[0], cv[1] + temperatureVel[1], cv[2]};
        std::array<double, 3> particlePosition = {cx[0] + h * nx, cx[1] + h * ny, cx[2] + h * nz};
        particles.addParticle(particlePosition, particleVelocity, m);
      }
    }
  }
}

void ParticleGenerator::generateDisc(std::array<double, 3> cx, std::array<double, 3> cv, int rn, double h, double m,
                                     double t) const {
  double radius_sq = rn * rn * h * h;
  for (double px = -rn * h; px <= rn * h; px += h) {
    for (double py = -rn * h; py <= rn * h; py += h) {
      if (px * px + py * py <= radius_sq) {  // if inside the 'circle'
        std::array<double, 3> temperatureVelocity = maxwellBoltzmannDistributedVelocity(t, m, 2);
        std::array<double, 3> particleVelocity = {cv[0] + temperatureVelocity[0], cv[1] + temperatureVelocity[1], cv[2]};
        std::array<double, 3> particlePosition = {cx[0] + px, cx[1] + py, cx[2]};
        particles.addParticle(particlePosition, particleVelocity, m);
      }
    }
  }
}
