#include <gtest/gtest.h>

#include <cmath>

#include "ForceCalc.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

class ForceCalcTest : public ::testing::Test {
 protected:
  ParticleContainer pc;
};

// Tests if overflow errors is thrown when the calculations are run on particles with the same coords.
TEST_F(ForceCalcTest, ExpectNormError) {
  pc.addParticle(std::array<double, 3>{0.}, std::array<double, 3>{0.}, 0.);
  pc.addParticle(std::array<double, 3>{0.}, std::array<double, 3>{1.}, 0.);
  EXPECT_THROW(GravityForce(pc).calculateF(), std::overflow_error);
  EXPECT_THROW(LennardJonesForce(pc, 1., 1., INFINITY).calculateF(), std::overflow_error);
}

// Tests the gravitational force between two particles if one particle has zero mass
TEST_F(ForceCalcTest, GravityF_ZeroMass) {
  auto p1 = Particle(0);
  auto p2 = Particle(std::array<double, 3>{1., 1., 1.}, std::array<double, 3>{0.}, 1.);
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  GravityForce(pc).calculateF();
  EXPECT_EQ(pc[0].getF(), (std::array<double, 3>{0.}));
}

// Tests the gravitational force between two particles with valid mass
TEST_F(ForceCalcTest, GravityF_TwoBody) {
  auto p1 = Particle(std::array<double, 3>{0., 0., 0.}, std::array<double, 3>{0.}, 1.);
  auto p2 = Particle(std::array<double, 3>{1., 1., 1.}, std::array<double, 3>{0.}, 1.);
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  double fr = 1. / pow(sqrt(3.), 3);
  std::array<double, 3> result = {fr, fr, fr};

  GravityForce(pc).calculateF();
  for (int i = 0; i < pc.size(); i++) {
    EXPECT_DOUBLE_EQ(pc[0].getF()[i], result[i]);
    EXPECT_DOUBLE_EQ(pc[1].getF()[i], -1. * result[i]);
  }
}

// Tests the gravitational force between two particles with a valid mass
TEST_F(ForceCalcTest, GravityF_TwoBody2) {
  auto p1 = Particle(std::array<double, 3>{10., 20., 30.}, std::array<double, 3>{1., 2., 3.}, 1000.);
  auto p2 = Particle(
      std::array<double, 3>{
          103.,
          202.,
          301,
      },
      std::array<double, 3>{500., 500., 500.}, 10000.);

  pc.addParticle(&p1);
  pc.addParticle(&p2);

  double norm_inv = 1. / ArrayUtils::L2Norm(p2.getX() - p1.getX());
  std::array<double, 3> F = 1000 * 10000 * norm_inv * norm_inv * norm_inv * (p2.getX() - p1.getX());
  GravityForce(pc).calculateF();
  for (int i = 0; i < pc.size(); i++) {
    EXPECT_DOUBLE_EQ(pc[0].getF()[i], F[i]);
    EXPECT_DOUBLE_EQ(pc[1].getF()[i], -1. * F[i]);
  }
}

// Tests the Lennard-Jones-Force calculation between two particles with valid arguments up to an error of 10e-6 simulation units
TEST_F(ForceCalcTest, LJ_F_TwoBody) {
  // LJ-Potential factor precomputed using WolframAlpha
  constexpr double factor = -41145. / 13176688.;

  auto p1 = Particle(std::array<double, 3>{0., 0., 0}, std::array<double, 3>{0., 0., 0.}, 1.);
  auto p2 = Particle(
      std::array<double, 3>{
          1.,
          2.,
          3,
      },
      std::array<double, 3>{100., 200., 300.}, 100.);

  pc.addParticle(&p1);
  pc.addParticle(&p2);

  const std::array<double, 3> F = factor * (p1.getX() - p2.getX());
  LennardJonesForce(pc, 5, 1, INFINITY).calculateF();
  for (int i = 0; i < pc.size(); i++) {
    EXPECT_NEAR(pc[0].getF()[i], F[i], 10e-6);
    EXPECT_NEAR(pc[1].getF()[i], -1. * F[i], 10e-6);
  }
}
