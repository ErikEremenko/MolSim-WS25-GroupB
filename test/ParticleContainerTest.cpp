#include <gtest/gtest.h>
#include "physics/ParticleContainer.h"

// Check if ParticleContainer saves new particles
class ParticleContainerTest : public ::testing::Test {
 protected:
  ParticleContainer pc;
};

TEST_F(ParticleContainerTest, IsInitiallyEmpty) {
  EXPECT_EQ(pc.size(), 0);
}

TEST_F(ParticleContainerTest, SizeAfterAdd) {
  Particle p(0);
  pc.addParticle(&p);
  EXPECT_EQ(pc.size(), 1);
}

TEST_F(ParticleContainerTest, AddMultiple) {
  Particle p1(0);
  Particle p2(1);
  pc.addParticle(&p1);
  pc.addParticle(&p1);
  pc.addParticle(&p2);
  EXPECT_EQ(pc.size(), 3);
}

TEST_F(ParticleContainerTest, AccessParticles) {
  Particle p1(0);
  p1.setX({1., 2., 3.});
  pc.addParticle(&p1);

  EXPECT_EQ(pc[0].getX()[0], 1.);
  EXPECT_EQ(pc[0].getX()[1], 2.);
  EXPECT_EQ(pc[0].getX()[2], 3.);
}

// --- Checkpoint addParticle Tests ---

// Test adding a particle with full checkpoint state
TEST_F(ParticleContainerTest, AddParticleWithFullState) {
  std::array<double, 3> x = {10.0, 20.0, 30.0};
  std::array<double, 3> v = {1.0, 2.0, 3.0};
  std::array<double, 3> f = {0.1, 0.2, 0.3};
  std::array<double, 3> oldF = {0.05, 0.1, 0.15};
  double m = 2.5;
  int type = 1;
  double sigma = 1.5;
  double epsilon = 4.0;

  pc.addParticle(x, v, m, f, oldF, type, sigma, epsilon);

  EXPECT_EQ(pc.size(), 1);
}

// Test that checkpoint particle has correct position and velocity
TEST_F(ParticleContainerTest, CheckpointParticlePositionVelocity) {
  std::array<double, 3> x = {10.0, 20.0, 30.0};
  std::array<double, 3> v = {1.0, 2.0, 3.0};
  std::array<double, 3> f = {0.1, 0.2, 0.3};
  std::array<double, 3> oldF = {0.05, 0.1, 0.15};

  pc.addParticle(x, v, 1.0, f, oldF, 0, 1.0, 5.0);

  EXPECT_DOUBLE_EQ(pc[0].getX()[0], 10.0);
  EXPECT_DOUBLE_EQ(pc[0].getX()[1], 20.0);
  EXPECT_DOUBLE_EQ(pc[0].getX()[2], 30.0);
  EXPECT_DOUBLE_EQ(pc[0].getV()[0], 1.0);
  EXPECT_DOUBLE_EQ(pc[0].getV()[1], 2.0);
  EXPECT_DOUBLE_EQ(pc[0].getV()[2], 3.0);
}

// Test that checkpoint particle has correct force vectors
TEST_F(ParticleContainerTest, CheckpointParticleForces) {
  std::array<double, 3> x = {0.0, 0.0, 0.0};
  std::array<double, 3> v = {0.0, 0.0, 0.0};
  std::array<double, 3> f = {0.1, 0.2, 0.3};
  std::array<double, 3> oldF = {0.05, 0.1, 0.15};

  pc.addParticle(x, v, 1.0, f, oldF, 0, 1.0, 5.0);

  // Check current force
  EXPECT_DOUBLE_EQ(pc[0].getF()[0], 0.1);
  EXPECT_DOUBLE_EQ(pc[0].getF()[1], 0.2);
  EXPECT_DOUBLE_EQ(pc[0].getF()[2], 0.3);

  // Check old force
  EXPECT_DOUBLE_EQ(pc[0].getOldF()[0], 0.05);
  EXPECT_DOUBLE_EQ(pc[0].getOldF()[1], 0.1);
  EXPECT_DOUBLE_EQ(pc[0].getOldF()[2], 0.15);
}

// Test that checkpoint particle has correct mass, type, sigma, epsilon
TEST_F(ParticleContainerTest, CheckpointParticleProperties) {
  std::array<double, 3> x = {0.0, 0.0, 0.0};
  std::array<double, 3> v = {0.0, 0.0, 0.0};
  std::array<double, 3> f = {0.0, 0.0, 0.0};
  std::array<double, 3> oldF = {0.0, 0.0, 0.0};
  double m = 2.5;
  int type = 3;
  double sigma = 1.5;
  double epsilon = 4.0;

  pc.addParticle(x, v, m, f, oldF, type, sigma, epsilon);

  EXPECT_DOUBLE_EQ(pc[0].getM(), 2.5);
  EXPECT_EQ(pc[0].getType(), 3);
  EXPECT_DOUBLE_EQ(pc[0].getSigma(), 1.5);
  EXPECT_DOUBLE_EQ(pc[0].getEpsilon(), 4.0);
}

// Test adding multiple checkpoint particles
TEST_F(ParticleContainerTest, AddMultipleCheckpointParticles) {
  std::array<double, 3> x1 = {1.0, 0.0, 0.0};
  std::array<double, 3> x2 = {2.0, 0.0, 0.0};
  std::array<double, 3> v = {0.0, 0.0, 0.0};
  std::array<double, 3> f = {0.0, 0.0, 0.0};
  std::array<double, 3> oldF = {0.0, 0.0, 0.0};

  pc.addParticle(x1, v, 1.0, f, oldF, 0, 1.0, 5.0);
  pc.addParticle(x2, v, 2.0, f, oldF, 1, 1.2, 4.0);

  EXPECT_EQ(pc.size(), 2);
  EXPECT_DOUBLE_EQ(pc[0].getX()[0], 1.0);
  EXPECT_DOUBLE_EQ(pc[1].getX()[0], 2.0);
  EXPECT_DOUBLE_EQ(pc[0].getM(), 1.0);
  EXPECT_DOUBLE_EQ(pc[1].getM(), 2.0);
}