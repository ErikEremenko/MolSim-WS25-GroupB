#include <gtest/gtest.h>

#include "io/YAMLFileReader.h"
#include "ParticleContainer.h"

class YAMLFileReaderTest : public ::testing::Test {
 protected:
  std::string project_dir = PROJ_SRC_DIR;
  // Use an existing YAML file for testing
  std::string inputFilename = project_dir + "/input/collision2_reflective.yaml";
};

// Test that the reader successfully loads a valid YAML file
TEST_F(YAMLFileReaderTest, LoadValidFile) {
  EXPECT_NO_THROW(YAMLFileReader reader(inputFilename));
}

// Test simulation parameters are read correctly
TEST_F(YAMLFileReaderTest, ReadSimulationParameters) {
  YAMLFileReader reader(inputFilename);

  EXPECT_DOUBLE_EQ(reader.getTend(), 20.0);
  EXPECT_DOUBLE_EQ(reader.getDeltaT(), 0.0005);
  EXPECT_DOUBLE_EQ(reader.getEpsilon(), 5.0);
  EXPECT_DOUBLE_EQ(reader.getSigma(), 1.0);
  EXPECT_DOUBLE_EQ(reader.getCutoff(), 3.0);
}

// Test output parameters are read correctly
TEST_F(YAMLFileReaderTest, ReadOutputParameters) {
  YAMLFileReader reader(inputFilename);

  EXPECT_EQ(reader.getOutputBaseName(), "collision2_reflective");
  EXPECT_EQ(reader.getWriteFrequency(), 100);
}

// Test domain size is read correctly
TEST_F(YAMLFileReaderTest, ReadDomainSize) {
  YAMLFileReader reader(inputFilename);

  auto domainSize = reader.getDomainSize();
  EXPECT_DOUBLE_EQ(domainSize[0], 180.0);
  EXPECT_DOUBLE_EQ(domainSize[1], 90.0);
  EXPECT_DOUBLE_EQ(domainSize[2], 1.0);
}

// Test boundary types are read correctly
TEST_F(YAMLFileReaderTest, ReadBoundaryTypes) {
  YAMLFileReader reader(inputFilename);

  auto boundaries = reader.getBoundaryTypesRaw();
  EXPECT_EQ(boundaries[0], "REFLECTIVE");
  EXPECT_EQ(boundaries[1], "REFLECTIVE");
  EXPECT_EQ(boundaries[2], "REFLECTIVE");
  EXPECT_EQ(boundaries[3], "REFLECTIVE");
  EXPECT_EQ(boundaries[4], "REFLECTIVE");
  EXPECT_EQ(boundaries[5], "REFLECTIVE");
}

// Test reading particles from cuboids
TEST_F(YAMLFileReaderTest, ReadCuboidParticles) {
  YAMLFileReader reader(inputFilename);
  ParticleContainer pc;

  reader.readFile(pc);

  // collision2_reflective.yaml has 20*20 + 100*20 = 2400 particles in total
  EXPECT_EQ(pc.size(), 2400);
}

// Test outflow boundary file
TEST_F(YAMLFileReaderTest, ReadOutflowBoundaries) {
  std::string outflowFile = project_dir + "/input/collision2_outflow.yaml";
  YAMLFileReader reader(outflowFile);

  auto boundaries = reader.getBoundaryTypesRaw();
  EXPECT_EQ(boundaries[0], "OUTFLOW");
  EXPECT_EQ(boundaries[1], "OUTFLOW");
  EXPECT_EQ(boundaries[2], "OUTFLOW");
  EXPECT_EQ(boundaries[3], "OUTFLOW");
  EXPECT_EQ(boundaries[4], "OUTFLOW");
  EXPECT_EQ(boundaries[5], "OUTFLOW");
}
