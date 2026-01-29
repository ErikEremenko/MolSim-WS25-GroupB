#pragma once

#include "io/YAMLFileReader.h"
#include "physics/ForceCalc.h"
#include "physics/ParticleGenerator.h"
#include "physics/Thermostat.h"
#include "simulation/SimulationConfig.h"

#include <memory>
#include <string>

/**
 * @class Simulation
 * Base class providing the core structure for building particle simulations.
 *
 * Holds global simulation parameters such as time step length and end time,
 * a particle container and a force acting on the particles, and provides the main run-loop
 * and output functionality.
 *
 * Other classes can inherit from this class to build more specific simulations, but this is strongly discouraged.
 */
class Simulation {
 protected:
  /**
   * @brief Total simulated time.
   * @note The simulation ends when the internal time counter passes this value.
   */
  double endTime;

  /**
   * @brief Time step size used for integration.
   */
  double dt;

  /**
   * @brief Current simulation time (updated each iteration).
   * Used by time-dependent forces (e.g., ConstantForce that stops after a certain time).
   */
  double currentTime;

  /**
   * @brief Start value of the internal time counter, this is usually 0 when starting a new simulation.
   * @note A checkpoint file will have a starting time greater than 0, e.g. 10.
   */
  double startTime;
  /**
   * @brief Starting iteration number, this is usually 0 when starting a new simulation.
   * @note A checkpoint file will have a starting iteration greater than 0, e.g. 15000.
   */
  int startIteration;

  // Simulation parameters needed for checkpointing
  // TODO: Move these to a CheckpointInfo struct, as they are not relevant to this class
  // TODO: A different idea is to have Simulation own a CheckpointWriter object
  double epsilon;
  double sigma;
  double cutoff;
  double gravity;
  int dimensions;
  std::array<double, 3> domainSize;
  std::array<std::string, 6> boundaryTypeStrings;

  /**
   * @brief Membrane Y dimension for force target index calculation.
   */
  int membraneDimY;

  /**
   * @brief Chosen simulation execution mode (benchmark/file output).
   */
  SimulationMode simulationMode;

  /**
   * @brief Frequency defines after how many simulation steps the output is written to file.
   */
  int writeFrequency;

  /**
   * @brief Frequency defines after how many simulation steps a checkpoint is written to file.
   */
  int checkpointFrequency;

  std::unique_ptr<ParticleGenerator> particleGenerator;

  /**
   * @brief Defines the base name of the simulation output files.
   */
  std::string outputBasename;

  /**
   * @brief Strategy defining how particles are stored and accessed.
   */
  std::unique_ptr<ParticleContainer> particles;

  /**
   * @brief Forces acting on or between the particles.
   */
  std::vector<std::unique_ptr<ForceCalc>> forces;

  /**
   * @brief Original force configurations (needed for checkpoint writing).
   */
  std::vector<ForceConfig> forceConfigs;

  // Thermostat-related members, TODO: Add docstrings for them
  std::unique_ptr<Thermostat> thermostat;
  bool needToAutoSetTargetTemperature = false;
  std::optional<double> initialTemperature = std::nullopt;

  /**
   * @brief Outputs the state of the particles for visualization in ParaView.
   * @param iteration Current simulation step in ticks.
   */
  void plotParticles(int iteration) const;

  /**
   * @brief Writes a checkpoint file.
   * @param iteration Current simulation step.
   * @param time Current simulation time.
   */
  void writeCheckpoint(int iteration, double time) const;

  /**
   * @brief Creates/loads the particles in the simulation.
   *
   * Must be implemented by subclasses to define a specific simulation scenario.
   */
  virtual void setupSimulation();

  /**
   * @name Simulation run methods
   * @{
   * @brief Runs the simulation in benchmark mode (no file output).
   */
  virtual void runFileOutput();

  /**
   * @brief Runs the simulation in benchmark mode (no file output).
   * @}
   */
  virtual void runBenchmark();

 public:
  /**
   * @brief Constructs a simulation object from the given configuration.
   * @param config Simulation configuration
   */
  explicit Simulation(SimulationConfig& config);

  /**
   * @brief Destructor.
   */
  virtual ~Simulation();

  /**
   * @brief The entry-point of the simulation.
   */
  void run();
};
