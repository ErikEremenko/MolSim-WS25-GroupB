#pragma once

#include <array>
#include <vector>

#include "physics/ParticleContainer.h"
/**
 * @class ParticleGenerator
 * @brief Class for lazy particle generation
 *
 * Offers different methods to generate sets of particles in the domain area
 */
class ParticleGenerator {
 private:
  // TODO: Add docstrings to 'orders'
  struct CuboidOrder {
    std::array<double, 3> pos;
    std::array<double, 3> vel;
    std::array<int, 3> dim;
    double h;
    double m;
    double temp;
    int type;
    double sigma;
    double epsilon;
    bool isMembrane;
  };

  struct DiscOrder {
    std::array<double, 3> center;
    std::array<double, 3> vel;
    int r_n;
    double h;
    double m;
    double temp;
    int type;
    double sigma;
    double epsilon;
  };

  struct ExplicitOrder {
    std::array<double, 3> x;
    std::array<double, 3> v;
    double m;
    std::array<double, 3> f;
    std::array<double, 3> oldF;
    int type;
    double sigma;
    double epsilon;
  };

  std::vector<CuboidOrder> cuboidOrders;
  std::vector<DiscOrder> discOrders;
  std::vector<ExplicitOrder> explicitOrders;  // TODO: Use linked list for this

  // Internal helper functions to keep generate() clean
  /**
   * @brief Generates the specified cuboid and populates the container with its particles
   * @param container Target container to populate
   * @param order Cuboid to generate
   */
  static void generateCuboid(ParticleContainer& container, const CuboidOrder& order);

  /**
   * @brief Generates the specified 2D-membrane and populates the container with its particles
   * @param container Target container to populate
   * @param order 2D-Membrane to generate
   */
  static void generateMembrane(ParticleContainer& container, const CuboidOrder& order);

  /**
   * @brief Generates the specified disc and populates the container with its particles
   * @param container Target container to populate
   * @param order Disc to generate
   */
  static void generateDisc(ParticleContainer& container, const DiscOrder& order);

 public:
  /**
   * @brief Initializes the particle generator with an empty list of particle orders.
   */
  ParticleGenerator();
  ~ParticleGenerator();

  /**
   * @brief Queues the generation of a cuboid of particles at the specified position.
   *
   * @param cx position vector of the lower left front-side corner of the cuboid.
   * @param cv velocity vector of the particles in the cuboid.
   * @param n number of particles in each dimension of the cuboid.
   * @param h distance between the particles in the cuboid.
   * @param m mass of a single particle in the cuboid.
   * @param t temperature of the cuboid.
   * @param type type of the particles in the cuboid.
   * @param sigma Lennard-Jones sigma for particles in the cuboid (default: 1.0).
   * @param epsilon Lennard-Jones epsilon for particles in the cuboid (default: 5.0).
   */
  void queueCuboid(std::array<double, 3> cx, std::array<double, 3> cv, std::array<int, 3> n, double h, double m,
                   double t, int type = 0, double sigma = 1.0, double epsilon = 5.0, bool isMembrane = true);

  /**
   * @brief Queues the generation of a disc of particles at the specified position.
   *
   * @param cx position vector of the center of the disc.
   * @param cv velocity vector of the particles in the disc.
   * @param rn radius of the disc in terms of particles along the radius.
   * @param h distance between the particles in the disc.
   * @param m mass of a single particle in the disc.
   * @param t temperature of the disc.
   * @param type type of the particles in the disc.
   * @param sigma Lennard-Jones sigma for particles in the disc (default: 1.0).
   * @param epsilon Lennard-Jones epsilon for particles in the disc (default: 5.0).
   */
  void queueDisc(std::array<double, 3> cx, std::array<double, 3> cv, int rn, double h, double m, double t, int type = 0,
                 double sigma = 1.0, double epsilon = 5.0);

  /**
   * @brief Queues the generation of a single particle.
   * TODO: Maybe use another ParticleContainer and don't lazy generate explicit particles?
   *
   * @param position position of the particle
   * @param velocity velocity of the particle
   * @param mass mass of the particle
   * @param force force acting on the particle in this time step
   * @param old_force force acting on the particle in the previous time step
   * @param type type of the particle
   * @param sigma sigma of the particle's substance
   * @param epsilon epsilon of the particle's substance
   */
  void queueParticle(std::array<double, 3> position, std::array<double, 3> velocity, double mass,
                     std::array<double, 3> force, std::array<double, 3> old_force, int type = 0, double sigma = 1.0,
                     double epsilon = 5.0);

  /**
   * @brief Executes all queued orders into the provided container.
   * @note Exhausts the order lists (vectors) and shrinks them to ensure minimal memory usage.
   * @param container Target container to populate
   */
  void generate(ParticleContainer& container);
};
