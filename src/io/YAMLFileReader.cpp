#include "io/YAMLFileReader.h"
#include "ParticleGenerator.h"

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
    const auto pos = sphere["position"].as<std::array<double, 3>>();
    const auto vel = sphere["velocity"].as<std::array<double, 3>>();
    const auto rn = sphere["radius_particles"].as<int>();
    const auto h = sphere["mesh_width"].as<double>();
    const auto m = sphere["mass"].as<double>();
    const auto meanV = sphere["mean_velocity"].as<double>();

    particleGenerator.generateDisc(pos, vel, rn, h, m, meanV);

    SPDLOG_DEBUG("Loaded sphere.");
  }

  // TODO: The issue here is that we already initialize the particles with temperatures - so T_init is irrelevant!
  // 1. Parse Required Thermostat Parameter: n_thermostat
  // We expect this to exist if the thermostat block exists.
  if (const auto& thermostatConfig = config["thermostat"]) {
    if (!thermostatConfig["n_thermostat"]) {
      SPDLOG_ERROR("Thermostat config found, but 'n_thermostat' is missing!");
      throw std::invalid_argument("Thermostat config not found.");
    }
    int n_thermostat = thermostatConfig["n_thermostat"].as<int>();

    double t_target = 0.0;
    if (thermostatConfig["T_target"]) {
      t_target = thermostatConfig["T_target"].as<double>();
    } else {
      SPDLOG_WARN("T_target not specified in thermostat. Defaulting to 0.");
    }

    // TODO: Add the rest of the thermostat parameters here
  } else {
    SPDLOG_INFO("No thermostat configuration found.");
  }
  // TODO: Add thermostat initialization now that we have thermostat config parsing
}
