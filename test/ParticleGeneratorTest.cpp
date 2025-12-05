#include <gtest/gtest.h>
#include <array>

#include "ParticleGenerator.h"

class ParticleGeneratorTest : public ::testing::Test {
 protected:
  ParticleContainer pc;
};

TEST_F(ParticleGeneratorTest, CuboidCorrectNrElements) {
  ParticleGenerator pg(pc);

  std::array<double, 3> cx = {1, 2, 3};
  std::array<double, 3> vx = {1, 1, 1};
  std::array<int, 3> n = {3, 4, 2};
  pg.generateCuboid(cx, vx, n, 1.225, 1, 5);
  EXPECT_EQ(pc.size(), n[0] * n[1] * n[2]);
}

TEST_F(ParticleGeneratorTest, AllParticlesWithinDisc) {
  ParticleGenerator pg(pc);
  double epsilon = 1e-9;

  std::array<double, 3> cx = {1, 2, 3};
  std::array<double, 3> vx = {1, 1, 1};
  std::array<int, 3> n = {3, 4, 2};
  pg.generateDisc(cx, vx, 5, 1.225, 1, 5);
  bool invalid = false;
  for (Particle p : pc) {
    std::array<double, 3> distv = {cx[0] - p.getX()[0], cx[1] - p.getX()[1], cx[2] - p.getX()[2]};
    double dist_sq = distv[0] * distv[0] + distv[1] * distv[1] + distv[2] * distv[2];
    if (5 * 5 * 1.225 * 1.225 < dist_sq * epsilon)
      invalid = true;
  }
  EXPECT_FALSE(invalid);
}
