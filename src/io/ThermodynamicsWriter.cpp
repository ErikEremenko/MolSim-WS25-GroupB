#include "io/ThermodynamicsWriter.h"

#include <stdexcept>

ThermodynamicsWriter::ThermodynamicsWriter(const std::string& baseName)
  : filename(outputDirectory + baseName + "_thermodynamics.csv") {}

void ThermodynamicsWriter::write(const double diffusion, const std::array<int, RDFCalculator::intervalCount>& rdfBins) const {
  const bool exists = fileExists();

  // Open file in append mode
  std::ofstream file(filename, std::ios::app);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  // Writer headers if writing for the first time
  if (!exists) {
    file << "Diffusion";
    for (size_t i = 1; i <= RDFCalculator::intervalCount; ++i) {
      file << ",Bin" << i;  // "bin1, bin2, bin3..."
    }
    file << "\n";
  }

  // Write diffusion and RDF values
  file << diffusion;
  for (const double val : rdfBins) {
    file << "," << val;
  }
  file << "\n";

  // File closes automatically
}
