#pragma once

#include <array>
#include <fstream>
#include <string>

#include "physics/thermodynamics/RDFCalculator.h"

// TODO: Write docstrings for this class
class ThermodynamicsWriter {
 private:
  inline static const std::string outputDirectory = "output/thermodynamics/";

  std::string filename;

  bool firstWrite = true;  // check if we are appending data or running a new experiment

 public:
  explicit ThermodynamicsWriter(const std::string& baseName);

  /**
     * @brief Writes the calculated thermodynamics data: Diffusion + RDF
     * @param simulationIterations Current iteration number of the simulation
     * @param temperature Current temperature of the system
     * @param diffusion The calculated diffusion coefficient (MSD)
     * @param rdfDensities The amount of particles in intervals - the radial distribution function (RDF)
     */
  void write(int simulationIterations, double temperature, double diffusion,
             const std::array<double, RDFCalculator::intervalCount>& rdfDensities);
};
