#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>

#include "io/CheckpointWriter.h"
#include "io/YAMLFileReader.h"
#include "physics/ParticleContainer.h"

class CheckpointWriterTest : public ::testing::Test {
 protected:
  std::string testOutputDir = "test/output";
  std::string testFilename = "test_writer_checkpoint.yaml";

  void SetUp() override {
    // Ensure output directory exists
    std::filesystem::create_directories(testOutputDir);
  }

  void TearDown() override {
    // Clean up test file
    std::string fullPath = testOutputDir + "/" + testFilename;
    if (std::filesystem::exists(fullPath)) {
      std::filesystem::remove(fullPath);
    }
    if (std::filesystem::exists(testOutputDir)) {
      std::filesystem::remove_all(testOutputDir);
    }
  }

  ParticleContainer createTestParticles() {
    ParticleContainer pc;
    std::array<double, 3> x1 = {1.0, 2.0, 3.0};
    std::array<double, 3> v1 = {0.1, 0.2, 0.3};
    std::array<double, 3> f1 = {0.01, 0.02, 0.03};
    std::array<double, 3> oldF1 = {0.005, 0.01, 0.015};
    pc.addParticle(x1, v1, 1.0, f1, oldF1, 0, 1.0, 5.0);

    std::array<double, 3> x2 = {4.0, 5.0, 6.0};
    std::array<double, 3> v2 = {0.4, 0.5, 0.6};
    std::array<double, 3> f2 = {0.04, 0.05, 0.06};
    std::array<double, 3> oldF2 = {0.02, 0.025, 0.03};
    pc.addParticle(x2, v2, 2.0, f2, oldF2, 1, 1.2, 4.0);

    return pc;
  }
};

// Test that writeCheckpoint creates a file
TEST_F(CheckpointWriterTest, CreatesFile) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  EXPECT_TRUE(std::filesystem::exists(fullPath));
}

// Test that checkpoint file is valid YAML and can be parsed
TEST_F(CheckpointWriterTest, CreatesValidYAML) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;

  // Should not throw when parsing
  EXPECT_NO_THROW({ YAML::Node config = YAML::LoadFile(fullPath); });
}

// Test that checkpoint metadata is correctly written
TEST_F(CheckpointWriterTest, WritesCheckpointMetadata) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 1234, 0.617, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  EXPECT_EQ(config["checkpoint"]["iteration"].as<int>(), 1234);
  EXPECT_DOUBLE_EQ(config["checkpoint"]["time"].as<double>(), 0.617);
}

// Test that output parameters are correctly written
TEST_F(CheckpointWriterTest, WritesOutputParameters) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "my_simulation", 50, 1000, 2.0, 0.001,
                                                  5.0, 1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  EXPECT_EQ(config["output"]["base_name"].as<std::string>(), "my_simulation");
  EXPECT_EQ(config["output"]["write_frequency"].as<int>(), 50);
  EXPECT_EQ(config["output"]["checkpoint_frequency"].as<int>(), 1000);
}

// Test that simulation parameters are correctly written
TEST_F(CheckpointWriterTest, WritesSimulationParameters) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 2.0, 0.001, 5.0,
                                                  1.0, 3.0, -9.81, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  EXPECT_DOUBLE_EQ(config["simulation"]["t_end"].as<double>(), 2.0);
  EXPECT_DOUBLE_EQ(config["simulation"]["delta_t"].as<double>(), 0.001);
  EXPECT_DOUBLE_EQ(config["simulation"]["epsilon"].as<double>(), 5.0);
  EXPECT_DOUBLE_EQ(config["simulation"]["sigma"].as<double>(), 1.0);
  EXPECT_DOUBLE_EQ(config["simulation"]["cutoff_radius"].as<double>(), 3.0);
  EXPECT_DOUBLE_EQ(config["simulation"]["gravity"].as<double>(), -9.81);
}

// Test that domain size is correctly written
TEST_F(CheckpointWriterTest, WritesDomainSize) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {150.0, 75.0, 2.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  auto size = config["domain"]["size"].as<std::array<double, 3>>();
  EXPECT_DOUBLE_EQ(size[0], 150.0);
  EXPECT_DOUBLE_EQ(size[1], 75.0);
  EXPECT_DOUBLE_EQ(size[2], 2.0);
}

// Test that boundary types are correctly written
TEST_F(CheckpointWriterTest, WritesBoundaryTypes) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "REFLECTIVE", "PERIODIC", "OUTFLOW", "REFLECTIVE", "PERIODIC"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  EXPECT_EQ(config["boundaries"]["x_min"].as<std::string>(), "OUTFLOW");
  EXPECT_EQ(config["boundaries"]["x_max"].as<std::string>(), "REFLECTIVE");
  EXPECT_EQ(config["boundaries"]["y_min"].as<std::string>(), "PERIODIC");
  EXPECT_EQ(config["boundaries"]["y_max"].as<std::string>(), "OUTFLOW");
  EXPECT_EQ(config["boundaries"]["z_min"].as<std::string>(), "REFLECTIVE");
  EXPECT_EQ(config["boundaries"]["z_max"].as<std::string>(), "PERIODIC");
}

// Test that correct number of particles are written
TEST_F(CheckpointWriterTest, WritesCorrectParticleCount) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  EXPECT_EQ(config["particles"].size(), 2);
}

// Test that particle data is correctly written
TEST_F(CheckpointWriterTest, WritesParticleData) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  // Check first particle
  auto p1 = config["particles"][0];
  auto x1 = p1["x"].as<std::array<double, 3>>();
  EXPECT_DOUBLE_EQ(x1[0], 1.0);
  EXPECT_DOUBLE_EQ(x1[1], 2.0);
  EXPECT_DOUBLE_EQ(x1[2], 3.0);

  auto v1 = p1["v"].as<std::array<double, 3>>();
  EXPECT_DOUBLE_EQ(v1[0], 0.1);
  EXPECT_DOUBLE_EQ(v1[1], 0.2);
  EXPECT_DOUBLE_EQ(v1[2], 0.3);

  EXPECT_DOUBLE_EQ(p1["m"].as<double>(), 1.0);
  EXPECT_EQ(p1["type"].as<int>(), 0);
  EXPECT_DOUBLE_EQ(p1["sigma"].as<double>(), 1.0);
  EXPECT_DOUBLE_EQ(p1["epsilon"].as<double>(), 5.0);
}

// Test that force vectors are correctly written
TEST_F(CheckpointWriterTest, WritesForceVectors) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 100, 0.5, "test_base", 10, 500, 1.0, 0.0005, 5.0,
                                                  1.0, 3.0, 0.0, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;
  YAML::Node config = YAML::LoadFile(fullPath);

  // Check first particle force vectors
  auto p1 = config["particles"][0];

  auto f1 = p1["f"].as<std::array<double, 3>>();
  EXPECT_DOUBLE_EQ(f1[0], 0.01);
  EXPECT_DOUBLE_EQ(f1[1], 0.02);
  EXPECT_DOUBLE_EQ(f1[2], 0.03);

  auto oldF1 = p1["oldF"].as<std::array<double, 3>>();
  EXPECT_DOUBLE_EQ(oldF1[0], 0.005);
  EXPECT_DOUBLE_EQ(oldF1[1], 0.01);
  EXPECT_DOUBLE_EQ(oldF1[2], 0.015);
}

// Test round-trip: write then read checkpoint
TEST_F(CheckpointWriterTest, RoundTripPreservesData) {
  ParticleContainer pc = createTestParticles();
  std::array<double, 3> domainSize = {100.0, 100.0, 1.0};
  std::array<std::string, 6> boundaryTypes = {"OUTFLOW", "REFLECTIVE", "OUTFLOW", "REFLECTIVE", "OUTFLOW", "OUTFLOW"};

  outputWriter::CheckpointWriter::writeCheckpoint(pc, testFilename, 500, 0.25, "roundtrip_test", 50, 200, 1.0, 0.0005,
                                                  5.0, 1.0, 3.0, -9.81, domainSize, boundaryTypes, testOutputDir);

  std::string fullPath = testOutputDir + "/" + testFilename;

  // Read back using YAMLFileReader
  YAMLFileReader reader(fullPath);
  auto config = reader.getConfig();

  // Verify checkpoint metadata via config
  EXPECT_EQ(config.startIteration, 500);
  EXPECT_DOUBLE_EQ(config.startTime, 0.25);

  // Verify parameters via config
  EXPECT_EQ(config.outputBasename, "roundtrip_test");
  EXPECT_EQ(config.writeFrequency, 50);
  EXPECT_EQ(config.checkpointFrequency, 200);
  // Note: gravity is no longer stored in SimulationConfig directly; it's part of force configs

  // Verify particles can be loaded
  ParticleContainer loadedPc;
  config.particleGenerator->generate(loadedPc);
  EXPECT_EQ(loadedPc.size(), 2);

  // Verify first particle data matches
  EXPECT_DOUBLE_EQ(loadedPc[0].getX()[0], 1.0);
  EXPECT_DOUBLE_EQ(loadedPc[0].getV()[0], 0.1);
  EXPECT_DOUBLE_EQ(loadedPc[0].getF()[0], 0.01);
  EXPECT_DOUBLE_EQ(loadedPc[0].getOldF()[0], 0.005);
  EXPECT_EQ(loadedPc[0].getType(), 0);
  EXPECT_DOUBLE_EQ(loadedPc[0].getSigma(), 1.0);
  EXPECT_DOUBLE_EQ(loadedPc[0].getEpsilon(), 5.0);

  // Verify second particle data matches (different type, sigma, epsilon)
  EXPECT_DOUBLE_EQ(loadedPc[1].getX()[0], 4.0);
  EXPECT_DOUBLE_EQ(loadedPc[1].getV()[0], 0.4);
  EXPECT_EQ(loadedPc[1].getType(), 1);
  EXPECT_DOUBLE_EQ(loadedPc[1].getSigma(), 1.2);
  EXPECT_DOUBLE_EQ(loadedPc[1].getEpsilon(), 4.0);
}
