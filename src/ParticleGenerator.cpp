
#include "../include/ParticleGenerator.h"
#include "../include/utils/MaxwellBoltzmannDistribution.h"

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#endif  // SPDLOG_ACTIVE_LEVEL
#include <spdlog/spdlog.h>

ParticleGenerator::ParticleGenerator(ParticleContainer& particles) : particles(particles) {}
ParticleGenerator::~ParticleGenerator() = default;

void ParticleGenerator::generateCuboid(std::array<double, 3> cx, std::array<double, 3> cv, std::array<int, 3> n,
                                       double h, double m, double t) {

  for (int nx = 0; nx < n[0]; nx++) {
    for (int ny = 0; ny < n[1]; ny++) {
      for (int nz = 0; nz < n[2]; nz++) {

        std::array<double, 3> tempx = {cx[0] + h * nx, cx[1] + h * ny, cx[2] + h * nz};
        std::array<double, 3> tempv = {cv[0], cv[1], cv[2]};
        std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(t, 2);
        tempv[0] += temperatureVel[0];
        tempv[1] += temperatureVel[1];
        tempv[2] += temperatureVel[2];
        particles.addParticle(tempx, tempv, m);
      }
    }
  }
}

void ParticleGenerator::generateDisc(std::array<double, 3> cx, std::array<double, 3> cv, int rn, double h, double m) {

  double radius_sq = rn * rn * h * h;

  for (double px = -rn * h; px <= rn * h; px += h) {
    for (double py = -rn * h; py <= rn * h; py += h) {

      //std::array<double, 3> tempv = {cv[0], cv[1], cv[2]};
      double dist_sq = px * px + py * py;
      if (dist_sq <= radius_sq) {
        std::array<double, 3> tempx = {cx[0] + px, cx[1] + py, cx[2]};
        particles.addParticle(tempx, cv, m);
      }
    }
  }
}
