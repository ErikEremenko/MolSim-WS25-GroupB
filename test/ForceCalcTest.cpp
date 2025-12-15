#include <gtest/gtest.h>
#include <cmath>

#include "ForceCalc.h"
#include "LinkedCellParticleContainer.h"
#include "ParticleContainer.h"
#include "utils/ArrayUtils.h"

class ForceCalcTest : public ::testing::Test {
 protected:
  ParticleContainer pc;
};

// Test if overflow error is thrown when the calculations are run on particles with the same coords.
TEST_F(ForceCalcTest, ExpectNormError) {
  pc.addParticle(std::array<double, 3>{0.}, std::array<double, 3>{0.}, 0.);
  pc.addParticle(std::array<double, 3>{0.}, std::array<double, 3>{1.}, 0.);
  EXPECT_THROW(GravityForce(pc).calculateF(), std::overflow_error);
  EXPECT_THROW(LennardJonesForce(pc, 1., 1., INFINITY).calculateF(), std::overflow_error);
}

// Test the gravitational force between two particles if one particle has zero mass
TEST_F(ForceCalcTest, GravityF_ZeroMass) {
  auto p1 = Particle(0);
  auto p2 = Particle(std::array<double, 3>{1., 1., 1.}, std::array<double, 3>{0.}, 1.);
  pc.addParticle(&p1);
  pc.addParticle(&p2);

  GravityForce(pc).calculateF();
  EXPECT_EQ(pc[0].getF(), (std::array<double, 3>{0.}));
}

// Test the gravitational force between two particles with valid mass
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

// Test the gravitational force between two particles with a valid mass
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

// Test the Lennard-Jones-Force calculation between two particles with valid arguments up to an error of 10e-6 simulation units
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
// Boundary Condition Tests using LinkedCellParticleContainer

class BoundaryConditionTest : public ::testing::Test {
 protected:
  std::array<double, 3> domainDims = {10.0, 10.0, 10.0};
  double cutoffRadius = 3.0;
  double epsilon = 5.0;
  double sigma = 1.0;
};

// Test that outflow boundary removes particles that move outside the domain
TEST_F(BoundaryConditionTest, OutflowRemovesParticles) {
  std::array<LinkedCellParticleContainer::BoundaryType, 6> outflowBoundaries = {
      LinkedCellParticleContainer::BoundaryType::OUTFLOW, LinkedCellParticleContainer::BoundaryType::OUTFLOW,
      LinkedCellParticleContainer::BoundaryType::OUTFLOW, LinkedCellParticleContainer::BoundaryType::OUTFLOW,
      LinkedCellParticleContainer::BoundaryType::OUTFLOW, LinkedCellParticleContainer::BoundaryType::OUTFLOW};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  // Particle inside domain
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle outside domain
  lpc.addParticle({-0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle at upper boundary (outside)
  lpc.addParticle({5.0, 10.5, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 3);

  lpc.applyBoundaryConditions();

  // Only the particle inside domain should remain
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[0], 5.0);
}

// Test that reflective boundaries apply repulsive force near walls
TEST_F(BoundaryConditionTest, ReflectiveAppliesForce) {
  std::array<LinkedCellParticleContainer::BoundaryType, 6> reflectiveBoundaries = {
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, reflectiveBoundaries);

  // Add particle close to the left wall (x = 0), repulsionDistance = 2^(1/6) * sigma ~ 1.1225
  lpc.addParticle({0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius);
  forceCalc.calculateF();

  // Particle close to wall should experience repulsive force pushing it away from wall -> positive force in x-direction (away from wall)
  EXPECT_GT(lpc[0].getF()[0], 0.0);
}

// Test mixed boundary conditions
TEST_F(BoundaryConditionTest, MixedBoundaries) {
  std::array<LinkedCellParticleContainer::BoundaryType, 6> mixedBoundaries = {
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::OUTFLOW,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::OUTFLOW,
      LinkedCellParticleContainer::BoundaryType::OUTFLOW,    LinkedCellParticleContainer::BoundaryType::OUTFLOW};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, mixedBoundaries);

  // Particle outside at x_max (outflow) should be removed
  lpc.addParticle({10.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle inside should stay
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 2);

  lpc.applyBoundaryConditions();

  // Particle outside outflow boundary should be removed
  EXPECT_EQ(lpc.size(), 1);
}

// Test that particles far from reflective walls don't experience extra forces
TEST_F(BoundaryConditionTest, ReflectiveNoForceWhenFar) {
  std::array<LinkedCellParticleContainer::BoundaryType, 6> reflectiveBoundaries = {
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, reflectiveBoundaries);

  // Add particle in the center (far from walls)
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius);
  forceCalc.calculateF();

  // Single particle in center should have zero force
  EXPECT_DOUBLE_EQ(lpc[0].getF()[0], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[1], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[2], 0.0);
}
