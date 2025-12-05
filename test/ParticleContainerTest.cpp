#include <gtest/gtest.h>
#include "ParticleContainer.h"

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