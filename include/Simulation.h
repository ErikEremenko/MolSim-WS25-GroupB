#pragma once

#include "ForceCalc.h"
#include "io/YAMLFileReader.h"

#include <memory>
#include <string>

/**
 * @enum SimulationMode
 * @brief This enum allows the user to choose if they want to benchmark the execution time or if
 * they want to get the output of the simulation in a file.
 */
enum class SimulationMode {
  BENCHMARK,
  FILE_OUTPUT
};

/**
 * @class BaseSimulation
 * Abstract base class providing the core structure for building particle simulations.
 *
 * Holds global simulation parameters such as time step length and end time,
 * a particle container and a force acting on the particles, and provides the main run-loop
 * and output functionality. Subclasses implement \ref setupSimulation "setupSimulation()"
 * to define a simulation scenario.
 */
class BaseSimulation {
protected:
  /**
   * @brief Total simulated time.
   */
  double end_time;

  /**
   * @brief Time step size used for integration.
   */
  double dt;

  /**
   * @brief Frequency defines after how many simulation steps the output is written to file.
   */
  int write_frequency;

  /**
   * @brief Defines the base name of the simulation output files.
   */
  std::string base_name;

  /**
   * @brief Strategy defining how particles are stored and accessed.
   */
  std::unique_ptr<ParticleContainer> particles;

  /**
   * @brief Strategy defining how forces are computed between particles.
   */
  std::unique_ptr<ForceCalc> forceCalc;

  /**
   * @brief Chosen simulation execution mode (benchmark/file output).
   */
  SimulationMode simulationMode;

  /**
   * @brief Outputs the state of the particles for visualization in ParaView.
   * @param iteration Current simulation step in ticks.
   * @param outputBaseName Name for the file output
   */
  void plotParticles(int iteration, const std::string& outputBaseName) const;

  /**
   * @brief Creates/loads the particles in the simulation.
   *
   * Must be implemented by subclasses to define a specific simulation scenario.
   */
  virtual void setupSimulation() = 0;

  /**
   * @name Simulation run methods
   * @{
   * @brief Runs the simulation in benchmark mode (no file output).
   */
  void runFileOutput(int frequency, const std::string& outputBaseName) const;

  /**
   * @brief Runs the simulation in benchmark mode (no file output).
   * @}
   */
  void runBenchmark() const;
public:
  /**
   * @brief Constructor for \ref Simulation.
   * @param end_time Total simulation time.
   * @param dt Time step size.
   * @param write_frequency
   * @param base_name
   * @param simulationMode Selected simulation mode.
   */
  BaseSimulation(double end_time, double dt, int write_frequency, const std::string& base_name,  SimulationMode simulationMode);

  /**
 * @brief Constructor for \ref Simulation used in \ref YAMLSimulation. End_time, Write_frequency, dt, base_name, are given defaults, but initialized in \ref YAMLSimulation::YAMLSimulation
 * @param simulationMode Selected simulation mode.
 */
  explicit BaseSimulation(SimulationMode simulationMode);
 virtual ~BaseSimulation();

  /**
   * @brief The entry-point of the simulation.
   */
  void run();
};

/**
 * @class CollisionSimulation
 * @brief Simulation setup for a collision scenario using cuboids from an input file.
 *
 * Implements the initialization of particles for the collision experiments in Assignment 2.
 */
class CollisionSimulation : public BaseSimulation {
private:
  /**
   * @brief Path to the input file containing the cuboids.
   */
  std::string inputFilename;
public:
  /**
   * @brief Constructor for \ref CollisionSimulation.
   * @param inputFilename Path to input cuboid file.
   * @param end_time Total simulation time.
   * @param dt Time step size.
   * @param simulationMode Selected simulation mode.
   */
  CollisionSimulation(std::string inputFilename, double end_time, double dt, SimulationMode simulationMode);
protected:
  /**
   * @brief Loads the cuboids from the input file and populates the particle container.
   */
  void setupSimulation() override;
};

class CollisionSimulationParallel final : public BaseSimulation {
private:
  std::string inputFilename;
public:
  CollisionSimulationParallel(std::string inputFilename, double end_time, double dt, SimulationMode simulationMode);
protected:
  void setupSimulation() override;
};

class YAMLSimulation final : public BaseSimulation {
private:
  std::string inputFilename;
  YAMLFileReader reader;
public:
  enum class ContainerKind {DIRECT, LINKED};
  enum class Parallelization {OFF, ON};
  YAMLSimulation(std::string inputFilename,
               SimulationMode simulationMode,
               ContainerKind kind = ContainerKind::LINKED,
               Parallelization parallelization = Parallelization::OFF
               );
protected:
  void setupSimulation() override;
};