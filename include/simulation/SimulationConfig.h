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
enum class ContainerType {  // TODO: Consider moving this to a separate header file
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

  // Particle initialization
  std::unique_ptr<ParticleGenerator> particleGenerator = std::make_unique<ParticleGenerator>();
  ;

  // File output
  int writeFrequency = 10;  // not used when benchmarking
  int checkpointFrequency = 0;
  std::string outputBasename = "MD_vtk";

  // Force parameters
  // TODO: Are these really all 'optional' parameters?
  std::optional<double> epsilon = std::nullopt;
  std::optional<double> sigma = std::nullopt;
  std::optional<double> cutoff = std::nullopt;
  std::optional<double> gravity = std::nullopt;

  bool useParallelization = false;

  // Container and Linked Cell parameters
  ContainerType containerType = ContainerType::DIRECT;
  std::optional<std::array<double, 3>> domainSize = std::nullopt;
  std::optional<std::array<BoundaryType, 6>> boundaryTypes = std::nullopt;

  // Thermostat
  std::optional<ThermostatConfig> thermostatConfig = std::nullopt;
};
