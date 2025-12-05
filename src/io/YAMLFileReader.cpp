#include "../../include/io/YAMLFileReader.h"
#include "../include/ParticleGenerator.h"

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
  if (!config["domain"] || !config["domain"]["size"]) {
    SPDLOG_ERROR("YAML file missing required key (domain.size)!");
    exit(-1);
  }
  if (!config["boundaries"]) {
    SPDLOG_ERROR("YAML file missing required key (boundaries)!");
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

std::array<double, 3> YAMLFileReader::getDomainSize() const {
  return config["domain"]["size"].as<std::array<double, 3>>();
}
std::array<std::string, 6> YAMLFileReader::getBoundaryTypesRaw() const {
  std::array<std::string, 6> b{};

  const auto node = config["boundaries"];
  b[0] = node["x_min"].as<std::string>();
  b[1] = node["x_max"].as<std::string>();
  b[2] = node["y_min"].as<std::string>();
  b[3] = node["y_max"].as<std::string>();
  b[4] = node["z_min"].as<std::string>();
  b[5] = node["z_max"].as<std::string>();
  return b;
}
double YAMLFileReader::getLJRepulsionDistance() const {
  if (const auto sim = config["simulation"]; sim["lj_repulsion_distance"]) {
    return sim["lj_repulsion_distance"].as<double>();
  }
  // defaults to 2^(1/6)*sigma if not specified
  return std::pow(2.0, 1.0 / 6.0) * getSigma();
}

void YAMLFileReader::readFile(ParticleContainer& particles) {

  ParticleGenerator particleGenerator(particles);

  const auto& cuboids = config["cuboids"];

  for (std::size_t i = 0; i < cuboids.size(); ++i) {
    const auto& cuboid = cuboids[i];

    // read Cuboid Parameters
    auto pos = cuboid["position"].as<std::array<double, 3>>();
    auto vel = cuboid["velocity"].as<std::array<double, 3>>();
    auto dim = cuboid["dimensions"].as<std::array<int, 3>>();
    const auto h = cuboid["mesh_width"].as<double>();
    const auto m = cuboid["mass"].as<double>();
    const auto meanV = cuboid["mean_velocity"].as<double>();

    particleGenerator.generateCuboid(pos, vel, dim, h, m, meanV);

    SPDLOG_DEBUG("Loaded cuboid {} with {} particles.", i, dim[0] * dim[1] * dim[2]);
  }

  const auto& spheres = config["spheres"];

  for (std::size_t i = 0; i < spheres.size(); ++i) {
    const auto& sphere = spheres[i];

    // read Sphere Parameters
    auto pos = sphere["position"].as<std::array<double, 3>>();
    auto vel = sphere["velocity"].as<std::array<double, 3>>();
    auto rn = sphere["radius_particles"].as<int>();
    const auto h = sphere["mesh_width"].as<double>();
    const auto m = sphere["mass"].as<double>();
    const auto meanV = sphere["mean_velocity"].as<double>();

    particleGenerator.generateDisc(pos, vel, rn, h, m, meanV);

    SPDLOG_DEBUG("Loaded sphere.");
  }
}
