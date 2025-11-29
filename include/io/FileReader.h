#pragma once

#include "ParticleContainer.h"

/**
 * @class BaseFileReader
 * @brief Abstract base class for reading particles from input files.
 *
 * Provides a common interface for file readers that load particles into
 * a ParticleContainer. Subclasses implement specific file formats or structures.
 */
class BaseFileReader {
protected:
  /**
   * @brief Path to the input file.
   */
  std::string filename;
public:
  /**
   * @brief Constructor for \ref BaseFileReader.
   * @param filename Path to the file that will be read.
   */
  explicit BaseFileReader(std::string filename);
  virtual ~BaseFileReader();

  /**
   * @brief Reads particle data from the file and inserts them into the container.
   *
   * This function should be overwritten by the inheriting classes.
   *
   * Default implementation reads a list of particles with:
   * XYZ position, XYZ velocity, and mass.
   *
   * @param particles ParticleContainer which will hold the particles read from the file.
   */
  virtual void readFile(ParticleContainer& particles);
};

/**
 * @class CuboidFileReader
 * @brief File reader that loads cuboids of particles.
 */
class CuboidFileReader : BaseFileReader {
public:
  using BaseFileReader::BaseFileReader;

  /**
   * @brief Reads the input file and accordingly populates ParticleContainer with particles forming cuboids.
   *
   * Each cuboid is expanded into particles on a 3D grid and assigned
   * velocities according to the Maxwell–Boltzmann distribution.
   *
   * @param particles ParticleContainer where particles are inserted.
   */
  void readFile(ParticleContainer& particles) override;
};
