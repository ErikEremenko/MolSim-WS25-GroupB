#pragma once

#include "io/CheckpointWriter.h"  // TODO: Re-implement checkpointing
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
class Simulation {  // TODO: For now, we keep the virtual methods for thermostat testing - make class final!
 protected:
  /**
   * @brief Total simulated time.
   */
  double endTime;

  /**
   * @brief Time step size used for integration.
   */
  double dt;

  /**
   * @brief Gravitational acceleration applied to particles.
   */
  // double gravity; // Moved below

  double startTime;
  int startIteration;

  // Simulation parameters needed for checkpointing
  double epsilon;
  double sigma;
  double cutoff;
  double gravity;
  int dimensions;
  std::array<double, 3> domainSize;
  std::array<std::string, 6> boundaryTypeStrings;

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

  // Thermostat-related members, TODO: Add docstrings for them
  std::unique_ptr<Thermostat> thermostat;
  bool needToAutoSetTargetTemperature = false;
  std::optional<double> initialTemperature = std::nullopt;

  /**
   * @brief Outputs the state of the particles for visualization in ParaView.
   * @param iteration Current simulation step in ticks.
   * @param outputBaseName Name for the file output
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
