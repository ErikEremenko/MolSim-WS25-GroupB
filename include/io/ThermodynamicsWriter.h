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
  [[nodiscard]] bool fileExists() const {
    const std::ifstream f(filename.c_str());
    return f.good();
  }

 public:
  explicit ThermodynamicsWriter(const std::string& baseName);

  /**
     * @brief Writes the calculated thermodynamics data: Diffusion + RDF
     * @param diffusion The calculated diffusion coefficient (MSD)
     * @param rdfBins The amount of particles in intervals - the radial distribution function (RDF)
     */
  void write(double diffusion, const std::array<int, RDFCalculator::intervalCount>& rdfBins) const;
};
