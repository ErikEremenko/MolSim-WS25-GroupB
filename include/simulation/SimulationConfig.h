#pragma once

#include <array>
#include <memory>
#include <optional>
#include <string>

#include "physics/LinkedCellParticleContainer.h"
#include "physics/ParticleGenerator.h"

/**
 * @enum SimulationMode
 * @brief Defines how the simulation is executed.
 *
 * Selects whether the simulation measures execution performance
 * or produces simulation results as output data.
 */
enum class SimulationMode {
  /**
   * @var SimulationMode::BENCHMARK
   * Runs the simulation without producing output and measures execution time.
   */
  BENCHMARK,

  /**
   * @var SimulationMode::FILE_OUTPUT
   * Runs the simulation and writes simulation results to an output file.
   */
  FILE_OUTPUT
};

/**
 * @enum ContainerType
 * @brief Specifies the container implementation used for particle storage.
 *
 * Determines the underlying data structure used to store and access particles.
 */
enum class ContainerType {
  /**
   * @var ContainerType::DIRECT
   * Uses a direct-access container (e.g. array-based storage).
   */
  DIRECT,

  /**
   * @var ContainerType::LINKED
   * Uses a linked-cell container for spatial partitioning.
   */
  LINKED
};

/**
 * @enum Parallelization
 * @brief Controls whether the simulation is executed using multithreading.
 */
enum class Parallelization {
  /**
   * @var Parallelization::OFF
   * Executes the simulation sequentially.
   */
  OFF,

  /**
   * @var Parallelization::ON
   * Executes the simulation in parallel using OpenMP.
   */
  ON
};

enum class ForceType {
  LENNARD_JONES,
  ACCELERATION,
  MEMBRANE_BONDS
  // TODO: Add the other forces here
};

struct ForceConfig {
  // Which force strategy was picked
  ForceType forceType = ForceType::LENNARD_JONES; // TODO: This could also be done with a std::type_index

  // Lennard-Jones potential
  std::optional<double> epsilon = std::nullopt;
  std::optional<double> sigma = std::nullopt;
  std::optional<double> cutoff = std::nullopt;  // cutoff for the force calculations

  // Acceleration
  std::optional<double> accX = std::nullopt;
  std::optional<double> accY = std::nullopt;
  std::optional<double> accZ = std::nullopt;

  // Membrane bonds
  std::optional<double> stiffnessConstant = std::nullopt;
  std::optional<double> bondLength = std::nullopt;

  // TODO: Constant force

  // TODO: Lennard-Jones truncated (or just use the truncated cutoff)

  // TODO: Or not to do? We could add inter-particular gravity but it's not used in any of the new simulations
};

struct ThermostatConfig {
  std::optional<double> tempInit = std::nullopt;
  int nThermostat;
  std::optional<double> tempTarget = std::nullopt;
  std::optional<double> tempDelta = std::nullopt;
};

struct SimulationConfig {
  // Basic simulation parameters
  double tEnd = 0.0;
  double deltaT = 0.0;
  int startIteration = 0;
  double startTime = 0.0;
  SimulationMode simulationMode = SimulationMode::FILE_OUTPUT;
  std::optional<int> dimensions = 3;  // Default to 3D

  // Particle initialization
  std::unique_ptr<ParticleGenerator> particleGenerator = std::make_unique<ParticleGenerator>();

  // File output
  int writeFrequency = 10;  // not used when benchmarking
  int checkpointFrequency = 0;
  std::string outputBasename = "MD_vtk";

  // Forces
  std::vector<ForceConfig> forceConfigs{};

  bool useParallelization = false;  // TODO: Replace this with compiler flag

  // Container and Linked Cell parameters
  ContainerType containerType = ContainerType::DIRECT;
  std::optional<double> linkedCellCutoff = std::nullopt;  // cutoff for linked cell
  std::optional<std::array<double, 3>> domainSize = std::nullopt;
  std::optional<std::array<BoundaryType, 6>> boundaryTypes = std::nullopt;

  // Thermostat
  std::optional<ThermostatConfig> thermostatConfig = std::nullopt;
};
