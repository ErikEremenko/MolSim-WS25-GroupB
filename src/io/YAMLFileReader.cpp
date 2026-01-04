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
  if (!config["output"] || !config["simulation"]) {
    SPDLOG_ERROR("YAML file missing required top level keys (output or simulation)!");
    exit(-1);
  }
  // Either cuboids or particles (checkpoint) must be present
  const bool hasCuboids = config["cuboids"] && config["cuboids"].size() > 0;
  const bool hasSpheres = config["spheres"] && config["spheres"].size() > 0;
  const bool hasParticles = config["particles"] && config["particles"].size() > 0;
  if (!hasCuboids && !hasSpheres && !hasParticles) {
    SPDLOG_ERROR("YAML file must contain at least one of: cuboids, spheres, or particles!");
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

int YAMLFileReader::getCheckpointFrequency() const {
  if (config["output"]["checkpoint_frequency"]) {
    return config["output"]["checkpoint_frequency"].as<int>();
  }
  return 0;  // Default to 0 (disabled) if not specified
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

bool YAMLFileReader::isCheckpoint() const {
  return config["particles"] && config["particles"].size() > 0;
}

int YAMLFileReader::getCheckpointIteration() const {
  if (config["checkpoint"] && config["checkpoint"]["iteration"]) {
    return config["checkpoint"]["iteration"].as<int>();
  }
  return 0;
}

double YAMLFileReader::getCheckpointTime() const {
  if (config["checkpoint"] && config["checkpoint"]["time"]) {
    return config["checkpoint"]["time"].as<double>();
  }
  return 0.0;
}

void YAMLFileReader::readFile(ParticleContainer& particles) {
  ParticleGenerator particleGenerator(particles);

  // Get global sigma/epsilon as defaults
  const double globalSigma = getSigma();
  const double globalEpsilon = getEpsilon();

  const auto& cuboids = config["cuboids"];

  for (std::size_t i = 0; i < cuboids.size(); ++i) {
    const auto& cuboid = cuboids[i];

    // Read Cuboid Parameters
    auto pos = cuboid["position"].as<std::array<double, 3>>();
    auto vel = cuboid["velocity"].as<std::array<double, 3>>();
    auto dim = cuboid["dimensions"].as<std::array<int, 3>>();
    const auto h = cuboid["mesh_width"].as<double>();
    const auto m = cuboid["mass"].as<double>();
    const auto meanV = cuboid["mean_velocity"].as<double>();

    // Per-object sigma/epsilon with fallback to the global values
    const double sigma = cuboid["sigma"] ? cuboid["sigma"].as<double>() : globalSigma;
    const double epsilon = cuboid["epsilon"] ? cuboid["epsilon"].as<double>() : globalEpsilon;

    particleGenerator.generateCuboid(pos, vel, dim, h, m, meanV, sigma, epsilon);

    SPDLOG_DEBUG("Loaded cuboid {} with {} particles (sigma={}, epsilon={}).", i, dim[0] * dim[1] * dim[2], sigma,
                 epsilon);
  }

  const auto& spheres = config["spheres"];

  for (std::size_t i = 0; i < spheres.size(); ++i) {
    const auto& sphere = spheres[i];

    // Read Sphere Parameters
    const auto pos = sphere["position"].as<std::array<double, 3>>();
    const auto vel = sphere["velocity"].as<std::array<double, 3>>();
    const auto rn = sphere["radius_particles"].as<int>();
    const auto h = sphere["mesh_width"].as<double>();
    const auto m = sphere["mass"].as<double>();
    const auto meanV = sphere["mean_velocity"].as<double>();

    // Per-object sigma/epsilon with fallback to the global values
    const double sigma = sphere["sigma"] ? sphere["sigma"].as<double>() : globalSigma;
    const double epsilon = sphere["epsilon"] ? sphere["epsilon"].as<double>() : globalEpsilon;

    particleGenerator.generateDisc(pos, vel, rn, h, m, meanV, sigma, epsilon);

    SPDLOG_DEBUG("Loaded sphere (sigma={}, epsilon={}).", sigma, epsilon);
  }

  // Checkpoint loading: if "particles" section exists, load individual particles
  // precedence over cuboid/sphere generation
  if (config["particles"] && config["particles"].size() > 0) {
    SPDLOG_INFO("Loading {} particles from checkpoint...", config["particles"].size());

    for (const auto& p : config["particles"]) {
      auto x = p["x"].as<std::array<double, 3>>();
      auto v = p["v"].as<std::array<double, 3>>();
      double m = p["m"].as<double>();

      // Force vectors (required for proper restart)
      std::array<double, 3> f = {0.0, 0.0, 0.0};
      std::array<double, 3> oldF = {0.0, 0.0, 0.0};
      if (p["f"]) {
        f = p["f"].as<std::array<double, 3>>();
      }
      if (p["oldF"]) {
        oldF = p["oldF"].as<std::array<double, 3>>();
      }

      // Type defaults to 0
      int type = 0;
      if (p["type"]) {
        type = p["type"].as<int>();
      }

      // Per-particle sigma/epsilon with fallback to global values
      double sigma = p["sigma"] ? p["sigma"].as<double>() : globalSigma;
      double epsilon = p["epsilon"] ? p["epsilon"].as<double>() : globalEpsilon;

      particles.addParticle(x, v, m, f, oldF, type, sigma, epsilon);
    }

    SPDLOG_DEBUG("Loaded {} particles from checkpoint.", config["particles"].size());
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
