#include "io/CheckpointWriter.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace outputWriter {

void CheckpointWriter::writeCheckpoint(const ParticleContainer& particles, const std::string& filename, int iteration,
                                       double currentTime, const std::string& baseName, int writeFrequency,
                                       int checkpointFrequency, double tEnd, double deltaT, double epsilon,
                                       double sigma, double cutoffRadius, double gravity,
                                       const std::array<double, 3>& domainSize,
                                       const std::array<std::string, 6>& boundaryTypes,
                                       const std::string& outputDirectory) {

  // Use the passed output directory
  const std::string& output_directory = outputDirectory;

  try {
    std::filesystem::create_directories(output_directory);
  } catch (const std::filesystem::filesystem_error& err) {
    SPDLOG_ERROR("Error creating checkpoint directory {}: {}", output_directory, err.what());
    return;
  }

  YAML::Emitter out;
  out << YAML::BeginMap;

  // Checkpoint metadata
  out << YAML::Key << "checkpoint" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "iteration" << YAML::Value << iteration;
  out << YAML::Key << "time" << YAML::Value << currentTime;
  out << YAML::EndMap;

  // Output configuration
  out << YAML::Key << "output" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "base_name" << YAML::Value << baseName;
  out << YAML::Key << "write_frequency" << YAML::Value << writeFrequency;
  out << YAML::Key << "checkpoint_frequency" << YAML::Value << checkpointFrequency;
  out << YAML::EndMap;

  // Simulation parameters
  out << YAML::Key << "simulation" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "t_end" << YAML::Value << tEnd;
  out << YAML::Key << "delta_t" << YAML::Value << deltaT;
  out << YAML::Key << "epsilon" << YAML::Value << epsilon;
  out << YAML::Key << "sigma" << YAML::Value << sigma;
  out << YAML::Key << "cutoff_radius" << YAML::Value << cutoffRadius;
  out << YAML::Key << "gravity" << YAML::Value << gravity;
  out << YAML::EndMap;

  // Domain configuration
  out << YAML::Key << "domain" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "size" << YAML::Value << YAML::Flow << YAML::BeginSeq << domainSize[0] << domainSize[1]
      << domainSize[2] << YAML::EndSeq;
  out << YAML::EndMap;

  // Boundary conditions
  out << YAML::Key << "boundaries" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "x_min" << YAML::Value << boundaryTypes[0];
  out << YAML::Key << "x_max" << YAML::Value << boundaryTypes[1];
  out << YAML::Key << "y_min" << YAML::Value << boundaryTypes[2];
  out << YAML::Key << "y_max" << YAML::Value << boundaryTypes[3];
  out << YAML::Key << "z_min" << YAML::Value << boundaryTypes[4];
  out << YAML::Key << "z_max" << YAML::Value << boundaryTypes[5];
  out << YAML::EndMap;

  // Empty cuboids and spheres (not needed for checkpoint, but required by format)
  out << YAML::Key << "cuboids" << YAML::Value << YAML::BeginSeq << YAML::EndSeq;
  out << YAML::Key << "spheres" << YAML::Value << YAML::BeginSeq << YAML::EndSeq;

  // Particles: complete phase space for each particle
  out << YAML::Key << "particles" << YAML::Value << YAML::BeginSeq;
  for (const auto& p : particles) {
    out << YAML::BeginMap;

    // Position
    out << YAML::Key << "x" << YAML::Value << YAML::Flow << YAML::BeginSeq << p.getX()[0] << p.getX()[1] << p.getX()[2]
        << YAML::EndSeq;

    // Velocity
    out << YAML::Key << "v" << YAML::Value << YAML::Flow << YAML::BeginSeq << p.getV()[0] << p.getV()[1] << p.getV()[2]
        << YAML::EndSeq;

    // Mass
    out << YAML::Key << "m" << YAML::Value << p.getM();

    // Current force
    out << YAML::Key << "f" << YAML::Value << YAML::Flow << YAML::BeginSeq << p.getF()[0] << p.getF()[1] << p.getF()[2]
        << YAML::EndSeq;

    // Old force (from previous iteration, needed for Störmer-Verlet)
    out << YAML::Key << "oldF" << YAML::Value << YAML::Flow << YAML::BeginSeq << p.getOldF()[0] << p.getOldF()[1]
        << p.getOldF()[2] << YAML::EndSeq;

    // Type
    out << YAML::Key << "type" << YAML::Value << p.getType();

    // Lennard-Jones parameters
    out << YAML::Key << "sigma" << YAML::Value << p.getSigma();
    out << YAML::Key << "epsilon" << YAML::Value << p.getEpsilon();

    out << YAML::EndMap;
  }
  out << YAML::EndSeq;

  out << YAML::EndMap;

  // Write to file
  std::string fullPath = output_directory + "/" + filename;
  std::ofstream fout(fullPath);
  if (!fout.is_open()) {
    SPDLOG_ERROR("Failed to open checkpoint file for writing: {}", fullPath);
    return;
  }

  fout << out.c_str();
  fout.close();

  SPDLOG_INFO("Checkpoint written: {} ({} particles, iteration {}, time {})",
              std::filesystem::absolute(fullPath).string(), particles.size(), iteration, currentTime);
}

}  // namespace outputWriter
