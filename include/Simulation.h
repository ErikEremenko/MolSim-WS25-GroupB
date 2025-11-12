#pragma once

#include "ForceCalc.h"

enum class SimulationMode {
  BENCHMARK,
  FILE_OUTPUT
};

class Simulation {
private:
  char* filename;
  const double end_time;
  const double dt;

  ParticleContainer& particles;
  ForceCalc& calcMethod;
  SimulationMode simulationMode;

  void plotParticles(int iteration) const;

  // Simulation run methods
  void runFileOutput() const;
  void runBenchmark() const;
public:
  Simulation(
    char* filename,
    double end_time,
    double dt,
    ParticleContainer& particles,
    ForceCalc& calcMethod,
    SimulationMode simulationMode
    );
  ~Simulation();

  void loadParticles() const;
  void run() const;
};
