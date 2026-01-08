#include <gtest/gtest.h>
#include <cmath>

#include "physics/ForceCalc.h"
#include "physics/LinkedCellParticleContainer.h"
#include "physics/ParticleContainer.h"
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
  EXPECT_THROW(LennardJonesForce(pc, 1., 1., INFINITY, 0).calculateF(), std::overflow_error);
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
  LennardJonesForce(pc, 5, 1, INFINITY, 0).calculateF();
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
  std::array<BoundaryType, 6> outflowBoundaries = {
      BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
      BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
      BoundaryType::OUTFLOW, BoundaryType::OUTFLOW};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  // Particle inside domain
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle outside domain
  lpc.addParticle({-0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle at upper boundary (outside)
  lpc.addParticle({5.0, 10.5, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 3);

  lpc.handleOutflowBoundaries();

  // Only the particle inside domain should remain
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[0], 5.0);
}

// Test that reflective boundaries apply repulsive force near walls
TEST_F(BoundaryConditionTest, ReflectiveAppliesForce) {
  std::array<BoundaryType, 6> reflectiveBoundaries = {
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE,
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE,
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, reflectiveBoundaries);

  // Add particle close to the left wall (x = 0), repulsionDistance = 2^(1/6) * sigma ~ 1.1225
  lpc.addParticle({0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius, 0);
  forceCalc.calculateF();

  // Particle close to wall should experience repulsive force pushing it away from wall -> positive force in x-direction (away from wall)
  EXPECT_GT(lpc[0].getF()[0], 0.0);
}

// Test mixed boundary conditions
TEST_F(BoundaryConditionTest, MixedBoundaries) {
  std::array<BoundaryType, 6> mixedBoundaries = {
      BoundaryType::REFLECTIVE, BoundaryType::OUTFLOW,
      BoundaryType::REFLECTIVE, BoundaryType::OUTFLOW,
      BoundaryType::OUTFLOW,    BoundaryType::OUTFLOW};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, mixedBoundaries);

  // Particle outside at x_max (outflow) should be removed
  lpc.addParticle({10.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle inside should stay
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 2);

  lpc.handleOutflowBoundaries();

  // Particle outside outflow boundary should be removed
  EXPECT_EQ(lpc.size(), 1);
}

// Test that particles far from reflective walls don't experience extra forces
TEST_F(BoundaryConditionTest, ReflectiveNoForceWhenFar) {
  std::array<BoundaryType, 6> reflectiveBoundaries = {
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE,
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE,
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, reflectiveBoundaries);

  // Add particle in the center (far from walls)
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius, 0);
  forceCalc.calculateF();

  // Single particle in center should have zero force
  EXPECT_DOUBLE_EQ(lpc[0].getF()[0], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[1], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[2], 0.0);
}

class PeriodicBoundaryTest : public ::testing::Test {
 protected:
  std::array<double, 3> domainDims = {10.0, 10.0, 10.0};
  double cutoffRadius = 3.0;
  double epsilon = 5.0;
  double sigma = 1.0;
  std::array<BoundaryType, 6> periodicBoundaries = {
      BoundaryType::PERIODIC, BoundaryType::PERIODIC,
      BoundaryType::PERIODIC, BoundaryType::PERIODIC,
      BoundaryType::PERIODIC, BoundaryType::PERIODIC};
};

// Test that particles outside left boundary (x < 0) wrap to right side
TEST_F(PeriodicBoundaryTest, ParticleWrapsFromLeftToRight) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Particle slightly outside left boundary
  lpc.addParticle({-0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 1);

  lpc.handleOutflowBoundaries();

  // Particle should be wrapped to righ side (-0.5 + 10.0 = 9.5)
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_NEAR(lpc[0].getX()[0], 9.5, 1e-10);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[1], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[2], 5.0);
}

// Test that particles outside right boundary (x > domain) wrap to the left side
TEST_F(PeriodicBoundaryTest, ParticleWrapsFromRightToLeft) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Particle slightly outside the right boundary
  lpc.addParticle({10.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 1);

  lpc.handleOutflowBoundaries();

  // Particle should be wrapped to left side (10.5 - 10.0 = 0.5)
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_NEAR(lpc[0].getX()[0], 0.5, 1e-10);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[1], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[2], 5.0);
}

// Test that particles wrap correctly in y-dimension
TEST_F(PeriodicBoundaryTest, ParticleWrapsInYDimension) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Particle outside bottom boundary
  lpc.addParticle({5.0, -1.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  lpc.handleOutflowBoundaries();

  // Should wrap to top: -1.0 + 10.0 = 9.0
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[0], 5.0);
  EXPECT_NEAR(lpc[0].getX()[1], 9.0, 1e-10);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[2], 5.0);
}

// Test that particles wrap correctly in z-dimension
TEST_F(PeriodicBoundaryTest, ParticleWrapsInZDimension) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Particle outside front boundary
  lpc.addParticle({5.0, 5.0, 11.0}, {0.0, 0.0, 0.0}, 1.0);

  lpc.handleOutflowBoundaries();

  // Should wrap to back (11.0 - 10.0 = 1.0)
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[0], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[1], 5.0);
  EXPECT_NEAR(lpc[0].getX()[2], 1.0, 1e-10);
}

// Test that particles at corners wrap correctly in multiple dimensions
TEST_F(PeriodicBoundaryTest, ParticleWrapsInMultipleDimensions) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Particle outside in both X and Y
  lpc.addParticle({-0.5, 10.5, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  lpc.handleOutflowBoundaries();

  // Should wrap in both dims
  EXPECT_EQ(lpc.size(), 1);
  EXPECT_NEAR(lpc[0].getX()[0], 9.5, 1e-10);
  EXPECT_NEAR(lpc[0].getX()[1], 0.5, 1e-10);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[2], 5.0);
}

// Test that two particles near opposite periodic boundaries interact correctly
TEST_F(PeriodicBoundaryTest, CrossBoundaryForceInteraction) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Two particles near opposite x boundaries (within cutoff distance across boundary)
  // Distance across periodic boundary: (10-9.5) + 0.5 = 1.0 < 3.0
  lpc.addParticle({0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  lpc.addParticle({9.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius, 0);
  forceCalc.calculateF();

  // Both particles should experience non-zero forces (periodic interaction)
  // Particle at 0.5 should be pulled towards negative x (toward the wrapped particle at 9.5)
  // Particle at 9.5 should be pulled towards positive x (toward the wrapped particle at 0.5)
  double f0_x = lpc[0].getF()[0];
  double f1_x = lpc[1].getF()[0];

  // Forces should be opposite (N3L)
  EXPECT_NEAR(f0_x + f1_x, 0.0, 1e-10);

  // At distance 1.0 with sigma=1.0, inside repulsion distance (d < 2^(1/6)*sigma ~ 1.12) -> forces should push particles apart
  EXPECT_GT(f0_x, 0.0);  // Particle 0 pushed toward +x (away from wrapped 9.5)
  EXPECT_LT(f1_x, 0.0);  // Particle 1 pushed toward -x (away from wrapped 0.5)
}

// Test that particles inside domain are not affected by periodic wrapping
TEST_F(PeriodicBoundaryTest, ParticlesInsideDomainUnchanged) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Particle well inside the domain
  lpc.addParticle({5.0, 5.0, 5.0}, {1.0, 2.0, 3.0}, 1.5);

  lpc.handleOutflowBoundaries();

  // Position should remain unchanged
  EXPECT_DOUBLE_EQ(lpc[0].getX()[0], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[1], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[2], 5.0);
}

// Test mixed boundary types: x:periodic, y:reflective, z:outflow
TEST_F(PeriodicBoundaryTest, MixedBoundaryWithPeriodic) {
  std::array<BoundaryType, 6> mixedBoundaries = {
      BoundaryType::PERIODIC,   BoundaryType::PERIODIC,
      BoundaryType::REFLECTIVE, BoundaryType::REFLECTIVE,
      BoundaryType::OUTFLOW,    BoundaryType::OUTFLOW};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, mixedBoundaries);

  // Particle outside X boundary (periodic) -> should wrap
  lpc.addParticle({-0.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle inside
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Particle outside Z boundary (outflow) -> should be removed
  lpc.addParticle({5.0, 5.0, 11.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 3);

  lpc.handleOutflowBoundaries();

  // First particle should wrap, second unchanged, third removed
  EXPECT_EQ(lpc.size(), 2);
}

// Test single particle with periodic boundaries has no force
TEST_F(PeriodicBoundaryTest, SingleParticleNoForce) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius, 0);
  forceCalc.calculateF();

  // Particle should have zero force
  EXPECT_DOUBLE_EQ(lpc[0].getF()[0], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[1], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[2], 0.0);
}

// Test particles too far across periodic boundary don't interact
TEST_F(PeriodicBoundaryTest, ParticlesBeyondCutoffNoInteraction) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, periodicBoundaries);

  // Placed s.t. minimum periodic distance > cutoff (cutoff: 3.0, domain: 10.0)
  // Particles at x = 2.0 and x = 6.0: direct distance = 4.0, periodic distance = 6.0
  // Both > cutoff ->no interaction
  lpc.addParticle({2.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  lpc.addParticle({6.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  LennardJonesForce forceCalc(lpc, epsilon, sigma, cutoffRadius, 0);
  forceCalc.calculateF();

  // Particles should experience zero force (beyond cutoff in all directions)
  EXPECT_DOUBLE_EQ(lpc[0].getF()[0], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[1], 0.0);
  EXPECT_DOUBLE_EQ(lpc[0].getF()[2], 0.0);
  EXPECT_DOUBLE_EQ(lpc[1].getF()[0], 0.0);
  EXPECT_DOUBLE_EQ(lpc[1].getF()[1], 0.0);
  EXPECT_DOUBLE_EQ(lpc[1].getF()[2], 0.0);
}
