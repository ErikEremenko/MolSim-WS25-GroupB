#include "physics/ParticleGenerator.h"

#include <cmath>

#include "utils/MaxwellBoltzmannDistribution.h"

// TODO: Why do we set the logging level here again? Isn't it enough in main.cpp?
#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#endif  // SPDLOG_ACTIVE_LEVEL
#include <spdlog/spdlog.h>

ParticleGenerator::ParticleGenerator() = default;
ParticleGenerator::~ParticleGenerator() = default;

void ParticleGenerator::generateCuboid(ParticleContainer& container, const CuboidOrder& order) {
  for (int nx = 0; nx < order.dim[0]; nx++) {
    for (int ny = 0; ny < order.dim[1]; ny++) {
      for (int nz = 0; nz < order.dim[2]; nz++) {
        std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(std::sqrt(order.temp / order.m), 3);
        std::array<double, 3> particleVelocity = {order.vel[0] + temperatureVel[0], order.vel[1] + temperatureVel[1],
                                                  order.vel[2] + temperatureVel[2]};
        std::array<double, 3> particlePosition = {order.pos[0] + order.h * nx, order.pos[1] + order.h * ny,
                                                  order.pos[2] + order.h * nz};
        container.addParticle(particlePosition, particleVelocity, order.m, order.type, order.sigma, order.epsilon);
      }
    }
  }
}

void ParticleGenerator::generateMembrane(ParticleContainer& container, const CuboidOrder& order) {
  // TODO: Maybe write tests for this?
  // --- Particle generation ---
  const int xDim = order.dim[0];
  const int yDim = order.dim[1];

  // Use a separate list to keep track of the newly created particles
  std::vector<Particle*> pointers;  // 2D-array flattened into 1D
  pointers.reserve(xDim * yDim);
  for (int nx = 0; nx < xDim; nx++) {
    for (int ny = 0; ny < yDim; ny++) {
      std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(std::sqrt(order.temp / order.m), 3);
      std::array<double, 3> particleVelocity = {order.vel[0] + temperatureVel[0], order.vel[1] + temperatureVel[1],
                                                order.vel[2] + temperatureVel[2]};
      std::array<double, 3> particlePosition = {order.pos[0] + order.h * nx, order.pos[1] + order.h * ny, order.pos[2]};
      Particle* ptr = container.addParticle(particlePosition, particleVelocity, order.m, order.type, order.sigma, order.epsilon);
      pointers.push_back(ptr);
    }
  }

  // --- Neighbor initialization ---
  const auto& addNeighbor = [xDim, yDim, pointers](std::vector<Particle*>& neighbors, const int x, const int y) {
    if (x >= 0 && x < xDim && y >= 0 && y < yDim) {  // out of bounds check
      neighbors.push_back(pointers[x*yDim + y]);
    }
  };

  for (int x = 0; x < order.dim[0]; x++) {
    for (int y = 0; y < order.dim[1]; y++) {
      Particle& particle = *pointers[x*yDim + y];

      // Direct neighbors
      auto& directNeighbors = particle.getDirectNeighbors();
      addNeighbor(directNeighbors, x-1, y);  // LEFT
      addNeighbor(directNeighbors, x+1, y);  // RIGHT
      addNeighbor(directNeighbors, x, y-1);  // TOP
      addNeighbor(directNeighbors, x, y+1);  // BOTTOM

      // Diagonal neighbors
      auto& diagonalNeighbors = particle.getDiagonalNeighbors();
      addNeighbor(diagonalNeighbors, x-1, y-1);  // TOP LEFT
      addNeighbor(diagonalNeighbors, x+1, y-1);  // TOP RIGHT
      addNeighbor(diagonalNeighbors, x-1, y+1);  // BOTTOM LEFT
      addNeighbor(diagonalNeighbors, x+1, y+1);  // BOTTOM RIGHT
    }
  }
}

void ParticleGenerator::generateDisc(ParticleContainer& container, const DiscOrder& order) {

  double radius_sq = order.r_n * order.r_n * order.h * order.h;

  for (int i = -order.r_n; i <= order.r_n; ++i) {
    for (int j = -order.r_n; j <= order.r_n; ++j) {
      double px = i * order.h;
      double py = j * order.h;
      if (px * px + py * py <= radius_sq) {
        std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(order.temp, 3);
        std::array<double, 3> tempv = {order.vel[0] + temperatureVel[0], order.vel[1] + temperatureVel[1],
                                       order.vel[2] + temperatureVel[2]};
        std::array<double, 3> tempx = {order.center[0] + px, order.center[1] + py, order.center[2]};
        container.addParticle(tempx, tempv, order.m, order.type, order.sigma, order.epsilon);
      }
    }
  }
}

void ParticleGenerator::queueCuboid(std::array<double, 3> position, std::array<double, 3> velocity,
                                    std::array<int, 3> dimensions, double mesh_width, double mass, double temperature,
                                    int type, double sigma, double epsilon, bool isMembrane) {
  cuboidOrders.push_back({position, velocity, dimensions, mesh_width, mass, temperature, type, sigma, epsilon, isMembrane});
}

void ParticleGenerator::queueDisc(std::array<double, 3> cx, std::array<double, 3> cv, int rn, double h, double m,
                                  double t, int type, double sigma, double epsilon) {
  discOrders.push_back({cx, cv, rn, h, m, t, type, sigma, epsilon});
}

void ParticleGenerator::queueParticle(std::array<double, 3> position, std::array<double, 3> velocity, double mass,
                                      std::array<double, 3> force, std::array<double, 3> old_force, int type,
                                      double sigma, double epsilon) {
  explicitOrders.push_back({position, velocity, mass, force, old_force, type, sigma, epsilon});
}

void ParticleGenerator::generate(ParticleContainer& container) {
  SPDLOG_INFO("Executing particle generation orders...");

  // Generate cuboids
  for (const auto& order : cuboidOrders) {
    if (order.isMembrane) {
      generateMembrane(container, order);
    } else {
      generateCuboid(container, order);
    }
  }
  cuboidOrders.clear();
  cuboidOrders.shrink_to_fit();

  // Generate discs
  for (const auto& order : discOrders) {
    generateDisc(container, order);
  }
  discOrders.clear();
  discOrders.shrink_to_fit();

  // Generate explicit particles
  for (const auto& p : explicitOrders) {
    container.addParticle(p.x, p.v, p.m, p.f, p.oldF, p.type, p.sigma, p.epsilon);
  }
  explicitOrders.clear();
  explicitOrders.shrink_to_fit();

  SPDLOG_INFO("Particle generation complete.");
}
