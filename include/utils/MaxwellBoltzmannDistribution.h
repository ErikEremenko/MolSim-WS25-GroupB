/*
 * MaxwellBoltzmannDistribution.h
 *
 * @Date: 13.12.2019
 * @Author: F. Gratl
 */

#pragma once

#include <array>
#include <random>
#include <cmath>

#define BOLTZMANN_CONSTANT 1.3806503e-23

/**
 * Generate a random velocity vector according to the Maxwell-Boltzmann distribution, with given temperature and mass.
 *
 * @param temp The temperature of the whole object in Kelvin
 * @param mass The mass of a single particle
 * @param dimensions Number of dimensions for which the velocity vector shall be generated. Set this to 2 or 3.
 * @return Array containing the generated velocity vector.
 */
inline std::array<double, 3> maxwellBoltzmannDistributedVelocity(
    double temp, double mass, size_t dimensions) {
  // we use a constant seed for repeatability.
  // random engine needs static lifetime otherwise it would be recreated for every call.
  static std::default_random_engine randomEngine(42069);

  // when adding independent normally distributed values to all velocity components
  // the velocity change is maxwell boltzmann distributed
  const double stddev = std::sqrt(BOLTZMANN_CONSTANT * temp / mass);
  std::normal_distribution<double> normalDistribution{0, stddev};
  std::array<double, 3> randomVelocity{};
  for (size_t i = 0; i < dimensions; ++i) {
    randomVelocity[i] = normalDistribution(randomEngine);
  }
  return randomVelocity;
}
