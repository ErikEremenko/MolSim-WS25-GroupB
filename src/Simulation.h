#pragma once

#include "FileCuboidReader.h"
#include "FileReader.h"
#include "ForceCalc.h"
#include "outputWriter/VTKWriter.h"

class Simulation {
private:
  char* filename;
  const double end_time;
  const double dt;

  ParticleContainer& particles;
  ForceCalc& calcMethod;

  void plotParticles(const int iteration) const;

public:
  Simulation(char* filename, const double end_time, const double dt,
             ParticleContainer& particles, ForceCalc& calcMethod);
  ~Simulation();

  void loadParticles() const;
  void run() const;
};
