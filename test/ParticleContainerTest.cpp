#include <gtest/gtest.h>

#include "ParticleContainer.h"

// Test if ParticleContainer saves new particles
TEST(ParticleContainerTest, SizeAfterAdd) {
  ParticleContainer container;
  container.addParticle(Particle(0));

  EXPECT_EQ(container.size(), 1);
}