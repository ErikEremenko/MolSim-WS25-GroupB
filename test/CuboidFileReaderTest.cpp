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
  EXPECT_EQ(pc[0].toString(), "Particle: X:[0, 0, 0] v: [-0.171411, 0.0178057, 0] f: [0, 0, 0] old_f: [0, 0, 0] type: 0");
  EXPECT_EQ(pc[1].toString(), "Particle: X:[0, 1.1225, 0] v: [0.00571789, -0.14098, 0] f: [0, 0, 0] old_f: [0, 0, 0] type: 0");
}