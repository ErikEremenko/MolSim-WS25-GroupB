#include <gtest/gtest.h>

#include "physics/ParticleContainer.h"
#include "io/YAMLFileReader.h"

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

  EXPECT_DOUBLE_EQ(reader.getTEnd(), 20.0);
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

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

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

// --- Checkpoint Tests ---

// Test that a non-checkpoint file is correctly identified as such and vice versa
TEST_F(YAMLFileReaderTest, IsCheckpointReturnsFalseForRegularFile) {
  YAMLFileReader reader(inputFilename);
  EXPECT_FALSE(reader.isCheckpoint());
}
TEST_F(YAMLFileReaderTest, IsCheckpointReturnsTrueForCheckpointFile) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  EXPECT_TRUE(reader.isCheckpoint());
}

// Test reading checkpoint metadata (iteration and time)
TEST_F(YAMLFileReaderTest, ReadCheckpointMetadata) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);

  EXPECT_EQ(reader.getCheckpointIteration(), 500);
  EXPECT_DOUBLE_EQ(reader.getCheckpointTime(), 0.25);
}

// Test that regular files return default checkpoint values
TEST_F(YAMLFileReaderTest, RegularFileReturnsDefaultCheckpointValues) {
  YAMLFileReader reader(inputFilename);

  // Regular files should return 0 for checkpoint iteration and time
  EXPECT_EQ(reader.getCheckpointIteration(), 0);
  EXPECT_DOUBLE_EQ(reader.getCheckpointTime(), 0.0);
  // Regular file doesn't have checkpoint_frequency, should return 0
  EXPECT_EQ(reader.getCheckpointFrequency(), 0);
}

// Test reading checkpoint frequency from checkpoint file
TEST_F(YAMLFileReaderTest, ReadCheckpointFrequency) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);

  EXPECT_EQ(reader.getCheckpointFrequency(), 200);
}

// Test reading individual particles from checkpoint file
TEST_F(YAMLFileReaderTest, ReadCheckpointParticles) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  ParticleContainer pc;

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

  // The test checkpoint has 3 particles
  EXPECT_EQ(pc.size(), 3);
}

// Test that checkpoint particles have correct position
TEST_F(YAMLFileReaderTest, CheckpointParticlePosition) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  ParticleContainer pc;

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

  // Check first particle position
  EXPECT_DOUBLE_EQ(pc[0].getX()[0], 10.0);
  EXPECT_DOUBLE_EQ(pc[0].getX()[1], 10.0);
  EXPECT_DOUBLE_EQ(pc[0].getX()[2], 0.0);

  // Check second particle position
  EXPECT_DOUBLE_EQ(pc[1].getX()[0], 20.0);
  EXPECT_DOUBLE_EQ(pc[1].getX()[1], 15.0);
  EXPECT_DOUBLE_EQ(pc[1].getX()[2], 0.0);
}

// Test that checkpoint particles have correct velocity
TEST_F(YAMLFileReaderTest, CheckpointParticleVelocity) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  ParticleContainer pc;

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

  // Check first particle velocity
  EXPECT_DOUBLE_EQ(pc[0].getV()[0], 1.0);
  EXPECT_DOUBLE_EQ(pc[0].getV()[1], 0.5);
  EXPECT_DOUBLE_EQ(pc[0].getV()[2], 0.0);
}

// Test that checkpoint particles have correct force vectors
TEST_F(YAMLFileReaderTest, CheckpointParticleForces) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  ParticleContainer pc;

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

  // Check first particle current force
  EXPECT_DOUBLE_EQ(pc[0].getF()[0], 0.1);
  EXPECT_DOUBLE_EQ(pc[0].getF()[1], 0.2);
  EXPECT_DOUBLE_EQ(pc[0].getF()[2], 0.0);

  // Check first particle old force
  EXPECT_DOUBLE_EQ(pc[0].getOldF()[0], 0.05);
  EXPECT_DOUBLE_EQ(pc[0].getOldF()[1], 0.1);
  EXPECT_DOUBLE_EQ(pc[0].getOldF()[2], 0.0);
}

// Test that checkpoint particles have correct mass and type
TEST_F(YAMLFileReaderTest, CheckpointParticleMassAndType) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  ParticleContainer pc;

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

  // First particle: mass=1.0, type=0
  EXPECT_DOUBLE_EQ(pc[0].getM(), 1.0);
  EXPECT_EQ(pc[0].getType(), 0);

  // Second particle: mass=2.0, type=1
  EXPECT_DOUBLE_EQ(pc[1].getM(), 2.0);
  EXPECT_EQ(pc[1].getType(), 1);
}

// Test that checkpoint particles have correct per-particle sigma/epsilon
TEST_F(YAMLFileReaderTest, CheckpointParticleSigmaEpsilon) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);
  ParticleContainer pc;

  auto config = reader.getConfig();
  config.particleGenerator->generate(pc);

  // First particle: sigma=1.0, epsilon=5.0
  EXPECT_DOUBLE_EQ(pc[0].getSigma(), 1.0);
  EXPECT_DOUBLE_EQ(pc[0].getEpsilon(), 5.0);

  // Second particle: sigma=1.2, epsilon=4.0
  EXPECT_DOUBLE_EQ(pc[1].getSigma(), 1.2);
  EXPECT_DOUBLE_EQ(pc[1].getEpsilon(), 4.0);
}

// Test reading checkpoint simulation parameters
TEST_F(YAMLFileReaderTest, CheckpointSimulationParameters) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);

  EXPECT_DOUBLE_EQ(reader.getTEnd(), 1.0);
  EXPECT_DOUBLE_EQ(reader.getDeltaT(), 0.0005);
  EXPECT_EQ(reader.getOutputBaseName(), "test_checkpoint");
  EXPECT_EQ(reader.getWriteFrequency(), 50);
}

// Test checkpoint with mixed boundary types
TEST_F(YAMLFileReaderTest, CheckpointMixedBoundaries) {
  std::string checkpointFile = project_dir + "/input/test_checkpoint.yaml";
  YAMLFileReader reader(checkpointFile);

  auto boundaries = reader.getBoundaryTypesRaw();
  EXPECT_EQ(boundaries[0], "OUTFLOW");
  EXPECT_EQ(boundaries[1], "OUTFLOW");
  EXPECT_EQ(boundaries[2], "REFLECTIVE");
  EXPECT_EQ(boundaries[3], "REFLECTIVE");
  EXPECT_EQ(boundaries[4], "OUTFLOW");
  EXPECT_EQ(boundaries[5], "OUTFLOW");
}
