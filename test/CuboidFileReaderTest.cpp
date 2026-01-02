#include <gtest/gtest.h>

#include <cmath>
#include "ParticleContainer.h"
#include "Simulation.h"
#include "io/FileReader.h"
#include "utils/ArrayUtils.h"
class CuboidFileReaderTest : public ::testing::Test {
 protected:
  std::string project_dir = PROJ_SRC_DIR;
  std::string inputFilename = project_dir + "/input/cuboid-test.txt";
  ParticleContainer pc;
};

// Tests if CuboidFileReader reads the input file as expected
TEST_F(CuboidFileReaderTest, BasicRead) {
  CuboidFileReader reader(inputFilename);
  reader.readFile(pc);

  ASSERT_EQ(pc.size(), 384);

  EXPECT_EQ(pc[0].getX()[0], 0.0);
  EXPECT_EQ(pc[0].getX()[1], 0.0);
  EXPECT_EQ(pc[0].getX()[2], 0.0);
  EXPECT_EQ(pc[0].getType(), 0);

  EXPECT_EQ(pc[1].getX()[0], 0.0);
  EXPECT_EQ(pc[1].getX()[1], 1.1225);
  EXPECT_EQ(pc[1].getX()[2], 0.0);
  EXPECT_EQ(pc[1].getType(), 0);

  // velocities are randomized, checking for sane values
  for (int i = 0; i < 3; ++i) {
    EXPECT_TRUE(std::isfinite(pc[0].getV()[i]));
    EXPECT_TRUE(std::isfinite(pc[1].getV()[i]));
  }
}