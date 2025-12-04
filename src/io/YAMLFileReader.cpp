#include "../../include/io/YAMLFileReader.h"

#include <spdlog/spdlog.h>

YAMLFileReader::YAMLFileReader(std::string filename) : BaseFileReader(filename) {
  try {
    config = YAML::LoadFile(filename);
    checkRequiredKeys();
    SPDLOG_INFO("Successfully loaded YAML file: {}", filename);
  } catch (const YAML::Exception& e) {
    SPDLOG_ERROR("Failed to parse YAML file: {}", e.what());
    exit(-1);
  }
}

void YAMLFileReader::checkRequiredKeys() const {
  if (!config["output"] || !config["simulation"] || !config["cuboids"]) {
    SPDLOG_ERROR("YAML file missing required top level keys (output, simulation, or cuboids)!");
    exit(-1);
  }
}

std::string YAMLFileReader::getOutputBaseName() const {
  return config["output"]["base_name"].as<std::string>();
}

int YAMLFileReader::getWriteFrequency() const {
  return config["output"]["write_frequency"].as<int>();
}

double YAMLFileReader::getTend() const {
  return config["simulation"]["t_end"].as<double>();
}

double YAMLFileReader::getDeltaT() const {
  return config["simulation"]["delta_t"].as<double>();
}

double YAMLFileReader::getEpsilon() const {
  return config["simulation"]["epsilon"].as<double>();
}

double YAMLFileReader::getSigma() const {
  return config["simulation"]["sigma"].as<double>();
}

double YAMLFileReader::getCutoff() const {
  return config["simulation"]["cutoff_radius"].as<double>();
}

void YAMLFileReader::readFile(ParticleContainer& particles) {
  const auto& cuboids = config["cuboids"];

  for (std::size_t i = 0; i < cuboids.size(); ++i) {
    const auto& cuboid = cuboids[i];

    // read Cuboid Parameters
    auto pos = cuboid["position"].as<std::vector<double>>();
    auto vel = cuboid["velocity"].as<std::vector<double>>();
    auto dim = cuboid["dimensions"].as<std::vector<int>>();
    const auto h = cuboid["mesh_width"].as<double>();
    const auto m = cuboid["mass"].as<double>();
    const auto meanV = cuboid["mean_velocity"].as<double>();

    // generate particles for cuboid
    for (int nx = 0; nx < dim[0]; nx++) {
      for (int ny = 0; ny < dim[1]; ny++) {
        for (int nz = 0; nz < dim[2]; nz++) {
          const std::array<double, 3> p_pos = {pos[0] + nx * h, pos[1] + ny * h, pos[2] + nz * h};

          std::array<double, 3> p_vel = {vel[0], vel[1], vel[2]};

          // Brownian Motion
          std::array<double, 3> brownian_vel = maxwellBoltzmannDistributedVelocity(meanV, 3);
          p_vel[0] += brownian_vel[0];
          p_vel[1] += brownian_vel[1];
          p_vel[2] += brownian_vel[2];

          particles.addParticle(p_pos, p_vel, m);
          SPDLOG_DEBUG("Generated particle with x={{{:.2f}, {:.2f}, {:.2f}}}, v={{{:.2f}, {:.2f}, {:.2f}}}, m={:.2f}",
                       p_pos[0], p_pos[1], p_pos[2], p_vel[0], p_vel[1], p_vel[2], m);
        }
      }
    }
    SPDLOG_DEBUG("Loaded cuboid {} with {} particles.", i, dim[0] * dim[1] * dim[2]);
  }
}
