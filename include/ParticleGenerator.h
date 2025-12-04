#pragma once

#include <array>

#include "../include/ParticleContainer.h"

class ParticleGenerator {
private:
    ParticleContainer& particles;
public:

    explicit ParticleGenerator(ParticleContainer& particles);
    ~ParticleGenerator();

    void generateDisc(std::array<double, 3> cx, std::array<double, 3> cv,
                        int rn, double h, double m);

    void generateCuboid(std::array<double, 3> cx, std::array<double, 3> cv,
                        std::array<int, 3> n, double h, double m,
                        double t);

};
