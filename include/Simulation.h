#pragma once

#include "ForceCalc.h"

#include <memory>
#include <string>

enum class SimulationMode {
  BENCHMARK,
  FILE_OUTPUT
};

class BaseSimulation {
protected:
  const double end_time;
  const double dt;

  std::unique_ptr<ParticleContainer> particles;
  std::unique_ptr<ForceCalc> forceCalc;
  SimulationMode simulationMode;

  void plotParticles(int iteration) const;

  virtual void setupSimulation() = 0;  // Implement this function to create new simulations

  // Simulation run methods
  void runFileOutput() const;
  void runBenchmark() const;
public:
  BaseSimulation(double end_time, double dt, SimulationMode simulationMode);
  virtual ~BaseSimulation();

  void run();
};

class CollisionSimulation : public BaseSimulation {
private:
  std::string inputFilename;
public:
  CollisionSimulation(std::string inputFilename, double end_time, double dt, SimulationMode simulationMode);
protected:
  void setupSimulation() override;
};

class CollisionSimulationParallel : public BaseSimulation {
private:
  std::string inputFilename;
public:
  CollisionSimulationParallel(std::string inputFilename, double end_time, double dt, SimulationMode simulationMode);
protected:
  void setupSimulation() override;
};