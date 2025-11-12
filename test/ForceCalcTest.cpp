#include <gtest/gtest.h>

#include <cmath>

#include "ForceCalc.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

// Tests if overflow errors is thrown when the calculations are run on particles with the same coords.
TEST(ForceCalcTest, ExpectNormError) {
  auto pc = ParticleContainer();
  pc.addParticle(std::array<double, 3>{0.}, std::array<double, 3>{0.}, 0.);
  pc.addParticle(std::array<double, 3>{0.}, std::array<double, 3>{1.}, 0.);

  GravityForce gf(pc);
  EXPECT_THROW(gf.calculateF(), std::overflow_error);

  LennardJonesForce ljf(pc, 1.0, 1.0);
  EXPECT_THROW(ljf.calculateF(), std::overflow_error);
}

// Tests the gravitational force between two particles if one particle has zero mass
TEST(ForceCalcTest, GravityF_ZeroMass) {
  auto p1 = Particle(0);
  auto p2 = Particle(std::array<double, 3>{1., 1., 1.},
                     std::array<double, 3>{0.}, 1.);
  auto pc = ParticleContainer();
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  GravityForce gf(pc);
  gf.calculateF();
  EXPECT_EQ(pc[0].getF(), (std::array<double, 3>{0.}));
}

// Tests the gravitational force between two particles with a valid mass
TEST(ForceCalcTest, GravityF_TwoBody) {
  auto p1 = Particle(std::array<double, 3>{0., 0., 0.},
                     std::array<double, 3>{0.}, 1.);
  auto p2 = Particle(std::array<double, 3>{1., 1., 1.},
                     std::array<double, 3>{0.}, 1.);

  auto pc = ParticleContainer();
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  double fr = 1. / pow(sqrt(3.), 3);
  std::array<double, 3> result = {fr, fr, fr};

  GravityForce gf(pc);
  gf.calculateF();
  for (int i = 0; i < pc.size(); i++) {
    EXPECT_DOUBLE_EQ(pc[0].getF()[i], result[i]);
    EXPECT_DOUBLE_EQ(pc[1].getF()[i], -1. * result[i]);
  }
}

// Tests the gravitational force between two particles with a valid mass
TEST(ForceCalcTest, GravityF_TwoBody2) {
  auto p1 = Particle(std::array<double, 3>{10., 20., 30.},
                     std::array<double, 3>{1., 2., 3.}, 1000.);
  auto p2 = Particle(
      std::array<double, 3>{
          103.,
          202.,
          301,
      },
      std::array<double, 3>{500., 500., 500.}, 10000.);

  auto pc = ParticleContainer();
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  double norm_inv = 1. / ArrayUtils::L2Norm(p2.getX() - p1.getX());
  std::array<double, 3> F =
      1000 * 10000 * norm_inv * norm_inv * norm_inv * (p2.getX() - p1.getX());
  GravityForce gf(pc);
  gf.calculateF();
  for (int i = 0; i < pc.size(); i++) {
    EXPECT_DOUBLE_EQ(pc[0].getF()[i], F[i]);
    EXPECT_DOUBLE_EQ(pc[1].getF()[i], -1. * F[i]);
  }
}

TEST(ForceCalcTest, LJ_F_TwoBody) {
  // precomputed factor
  double f = -41145./13176688.;

  auto p1 = Particle(std::array<double, 3>{0., 0., 0},
                     std::array<double, 3>{0., 0., 0.}, 1.);
  auto p2 = Particle(
      std::array<double, 3>{
          1.,
          2.,
          3,
      },
      std::array<double, 3>{100., 200., 300.}, 100.);

  auto pc = ParticleContainer();
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  std::array<double, 3> F = f * (p1.getX() - p2.getX());

  LennardJonesForce ljf(pc, 5.0, 1.0);
  ljf.calculateF();
  for (int i = 0; i < pc.size(); i++) {
    EXPECT_NEAR(pc[0].getF()[i], F[i], 10e-6);
    EXPECT_NEAR(pc[1].getF()[i], -1. * F[i], 10e-6);
  }
}
