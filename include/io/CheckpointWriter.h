#pragma once

#include <string>
#include "physics/ParticleContainer.h"

namespace outputWriter {

/**
 * @class CheckpointWriter
 * @brief Writes simulation state to a YAML checkpoint file for later restart.
 *
 * The checkpoint file contains the complete phase space (position, velocity,
 * force, old force, mass, type, sigma, epsilon) for each particle, plus
 * simulation metadata needed to resume the simulation from the checkpoint.
 */
class CheckpointWriter {
 public:
  /**
   * @brief Write a checkpoint file with the current simulation state.
   * @param particles The particle container with all particles
   * @param filename Output filename (should end in .yaml or .yml)
   * @param iteration Current simulation iteration
   * @param currentTime Current simulation time
   * @param baseName Output base name for resumed simulation
   * @param writeFrequency Write frequency for resumed simulation
   * @param tEnd End time for resumed simulation
   * @param deltaT Time step for resumed simulation
   * @param epsilon Global epsilon parameter
   * @param sigma Global sigma parameter
   * @param cutoffRadius Cutoff radius for interactions
   * @param domainSize Simulation domain size
   * @param boundaryTypes Boundary type strings (x_min, x_max, y_min, y_max, z_min, z_max)
   */
  static void writeCheckpoint(
      const ParticleContainer& particles,
      const std::string& filename,
      int iteration,
      double currentTime,
      const std::string& baseName,
      int writeFrequency,
      int checkpointFrequency,
      double tEnd,
      double deltaT,
      double epsilon,
      double sigma,
      double cutoffRadius,
      const std::array<double, 3>& domainSize,
      const std::array<std::string, 6>& boundaryTypes);
};

}  // namespace outputWriter
