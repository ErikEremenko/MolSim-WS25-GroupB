#pragma once

#include <array>

#include "../include/ParticleContainer.h"
/**
 * @class ParticleGenerator
 * @brief Class for particle generation
 *
 * Offers different methods to generate sets of particles in the domain area
 */
class ParticleGenerator {
private:
    ParticleContainer& particles;
public:
    /**
     * @brief Constructor for \ref ParticleGenerator.
     * @param particles ParticleContainer into which the particles would be generated.
     */
    explicit ParticleGenerator(ParticleContainer& particles);
    ~ParticleGenerator();

    /**
     * @brief Generates a disc of particles at the specified position.
     *
     * @param cx position vector of the center of the disc.
     * @param cv velocity vector of the particles in the disc.
     * @param rn radius of the disc in terms of particles along the radius.
     * @param h distance between the particles in the disc.
     * @param m mass of each particle in the disc.
     * @param t the mean value of the velocity of the Brownian Motion.
     */
    void generateDisc(std::array<double, 3> cx, std::array<double, 3> cv,
                        int rn, double h, double m, double t) const;

    /**
     * @brief Generates a cuboid of particles at the specified position.
     *
     * @param cx position vector of the lower left front-side corner of the cuboid.
     * @param cv velocity vector of the particles in the cuboid.
     * @param n number of particles in each dimension of the cuboid.
     * @param h distance between the particles in the cuboid.
     * @param m mass of each particle in the disc.
     * @param t the mean value of the velocity of the Brownian Motion.
     */
    void generateCuboid(std::array<double, 3> cx, std::array<double, 3> cv,
                        std::array<int, 3> n, double h, double m,
                        double t) const;

};
