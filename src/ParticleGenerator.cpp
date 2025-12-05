
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

  for (int i = -rn; i <= rn; ++i) {
    for (int j = -rn; j <= rn; ++j) {
      double px = i * h;
      double py = j * h;
      if (px * px + py * py <= radius_sq) {
        std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(t, m, 2);
        std::array<double, 3> tempv = {cv[0] + temperatureVel[0], cv[1] + temperatureVel[1], cv[2] + temperatureVel[2]};
        std::array<double, 3> tempx = {cx[0] + px, cx[1] + py, cx[2]};
        particles.addParticle(tempx, tempv, m);
      }
    }
  }
}
