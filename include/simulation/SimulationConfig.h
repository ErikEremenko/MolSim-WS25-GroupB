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
  LENNARD_JONES,      ///< Full Lennard-Jones potential
  SMOOTHED_LJ,        ///< Smoothed Lennard-Jones potential with continuous force at cutoff
  TRUNCATED_LJ,       ///< Repulsive-only LJ (truncated at 2^(1/6)*sigma)
  GLOBAL_GRAVITY,     ///< Constant gravitational acceleration on all particles
  HARMONIC_MEMBRANE,  ///< Harmonic bonds between membrane neighbors
  CONSTANT_FORCE      ///< Constant force on specific particles (e.g., pulling)
};

struct ForceConfig {
  ForceType forceType = ForceType::LENNARD_JONES;

  // Lennard-Jones potential (full, smoothed, or truncated)
  std::optional<double> epsilon = std::nullopt;
  std::optional<double> sigma = std::nullopt;
  std::optional<double> cutoff = std::nullopt;  // cutoff radius (r_c) for force calculations
  std::optional<double> rl = std::nullopt;      // smoothing start radius (r_l) for smoothed LJ

  // Global gravity (constant acceleration along specified axis)
  std::optional<double> gravity = std::nullopt;   ///< Gravity acceleration value
  std::optional<int> gravityAxis = std::nullopt;  ///< Axis: 0=x, 1=y (default), 2=z

  // Harmonic membrane bonds
  std::optional<double> stiffness = std::nullopt;      ///< Stiffness constant k
  std::optional<double> avgBondLength = std::nullopt;  ///< Average bond length r0

  // Constant force on specific particles
  std::optional<double> forceX = std::nullopt;
  std::optional<double> forceY = std::nullopt;
  std::optional<double> forceZ = std::nullopt;
  std::optional<double> endTime = std::nullopt;    ///< Time after which force stops
  std::vector<std::pair<int, int>> targetIndices;  ///< x/y indices of particles to apply force to
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
  ContainerType containerType = ContainerType::LINKED;
  std::optional<double> linkedCellCutoff = std::nullopt;  // cutoff for linked cell
  std::optional<std::array<double, 3>> domainSize = std::nullopt;
  std::optional<std::array<BoundaryType, 6>> boundaryTypes = std::nullopt;

  // Thermostat
  std::optional<ThermostatConfig> thermostatConfig = std::nullopt;

  // Thermodynamics statistics
  bool calculateThermodynamics = false;

  // Membrane parameters (for ConstantForce target index calculation)
  std::optional<int> membraneDimY = std::nullopt;  ///< Y-dimension of membrane grid
};
