#pragma once

#include "ForceCalc.h"

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
  const double end_time;

  /**
   * @brief Time step size used for integration.
   */
  const double dt;

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
   */
  void plotParticles(int iteration) const;

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
  void runFileOutput() const;

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
   * @param simulationMode Selected simulation mode.
   */
  BaseSimulation(double end_time, double dt, SimulationMode simulationMode);
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
