#include "io/CheckpointWriter.h"

#include <spdlog/spdlog.h>
#include <yaml-cpp/yaml.h>
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace outputWriter {

void CheckpointWriter::writeCheckpoint(const ParticleContainer& particles, const std::string& filename, int iteration,
                                       double currentTime, const std::string& baseName, int writeFrequency,
                                       int checkpointFrequency, double tEnd, double deltaT,
                                       const std::vector<ForceConfig>& forceConfigs,
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

  // Forces configuration - write all force types
  out << YAML::Key << "forces" << YAML::Value << YAML::BeginSeq;

  // Infer cutoff from first LJ-type force for container config
  double cutoffRadius = 3.0;  // default

  for (const auto& fc : forceConfigs) {
    out << YAML::BeginMap;
    switch (fc.forceType) {
      case ForceType::LENNARD_JONES:
        out << YAML::Key << "force_type" << YAML::Value << "lennard_jones";
        out << YAML::Key << "epsilon" << YAML::Value << fc.epsilon.value_or(1.0);
        out << YAML::Key << "sigma" << YAML::Value << fc.sigma.value_or(1.0);
        out << YAML::Key << "cutoff_radius" << YAML::Value << fc.cutoff.value_or(3.0);
        cutoffRadius = fc.cutoff.value_or(3.0);
        break;

      case ForceType::SMOOTHED_LJ:
        out << YAML::Key << "force_type" << YAML::Value << "smoothed_lj";
        out << YAML::Key << "epsilon" << YAML::Value << fc.epsilon.value_or(1.0);
        out << YAML::Key << "sigma" << YAML::Value << fc.sigma.value_or(1.0);
        out << YAML::Key << "cutoff_radius" << YAML::Value << fc.cutoff.value_or(3.0);
        out << YAML::Key << "smoothing_radius" << YAML::Value << fc.rl.value_or(1.8);
        cutoffRadius = fc.cutoff.value_or(3.0);
        break;

      case ForceType::TRUNCATED_LJ:
        out << YAML::Key << "force_type" << YAML::Value << "truncated_lj";
        break;

      case ForceType::GLOBAL_GRAVITY:
        out << YAML::Key << "force_type" << YAML::Value << "global_gravity";
        out << YAML::Key << "g" << YAML::Value << fc.gravity.value_or(0.0);
        if (fc.gravityAxis.has_value() && fc.gravityAxis.value() != 1) {
          out << YAML::Key << "axis" << YAML::Value << fc.gravityAxis.value();
        }
        break;

      case ForceType::HARMONIC_MEMBRANE:
        out << YAML::Key << "force_type" << YAML::Value << "harmonic_membrane";
        out << YAML::Key << "stiffness" << YAML::Value << fc.stiffness.value_or(300.0);
        out << YAML::Key << "avg_bond_length" << YAML::Value << fc.avgBondLength.value_or(2.2);
        break;

      case ForceType::CONSTANT_FORCE:
        out << YAML::Key << "force_type" << YAML::Value << "constant_force";
        out << YAML::Key << "fx" << YAML::Value << fc.forceX.value_or(0.0);
        out << YAML::Key << "fy" << YAML::Value << fc.forceY.value_or(0.0);
        out << YAML::Key << "fz" << YAML::Value << fc.forceZ.value_or(0.0);
        out << YAML::Key << "end_time" << YAML::Value << fc.endTime.value_or(0.0);
        if (!fc.targetIndices.empty()) {
          out << YAML::Key << "target_indices" << YAML::Value << YAML::BeginSeq;
          for (const auto& [x, y] : fc.targetIndices) {
            out << YAML::Flow << YAML::BeginSeq << x << y << YAML::EndSeq;
          }
          out << YAML::EndSeq;
        }
        break;
    }
    out << YAML::EndMap;
  }
  out << YAML::EndSeq;

  // Container configuration (required by YAMLFileReader for linked cell)
  out << YAML::Key << "container" << YAML::Value << YAML::BeginMap;
  out << YAML::Key << "container_type" << YAML::Value << "linked";
  out << YAML::Key << "cutoff" << YAML::Value << cutoffRadius;
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
