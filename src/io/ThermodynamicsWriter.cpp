#include "io/ThermodynamicsWriter.h"

#include <filesystem>
#include <stdexcept>

ThermodynamicsWriter::ThermodynamicsWriter(const std::string& baseName)
    : filename(outputDirectory + baseName + "_thermodynamics.csv") {
  // Create the directory if it doesn't exist
  try {
    std::filesystem::create_directories(outputDirectory);
  } catch (const std::filesystem::filesystem_error& err) {
    throw std::runtime_error("Could not create output directory " + outputDirectory + ": " + err.what());
  }
}

void ThermodynamicsWriter::write(const int simulationIterations, const double temperature, const double diffusion,
                                 const std::array<double, RDFCalculator::intervalCount>& rdfDensities) {
  auto writeMode = std::ios::app;  // append
  if (firstWrite) {
    writeMode = std::ios::trunc;  // overwrite
  }

  // Open/create file
  std::ofstream file(filename, writeMode);
  if (!file.is_open()) {
    throw std::runtime_error("Could not open file: " + filename);
  }

  // Write headers if writing for the first time
  if (firstWrite) {
    file << "Iterations,Temperature,Diffusion";
    for (size_t i = 1; i <= RDFCalculator::intervalCount; ++i) {
      file << ",Bin" << i;  // "bin1, bin2, bin3..."
    }
    file << "\n";

    firstWrite = false;
  }

  // Write data
  file << simulationIterations << "," << temperature << "," << diffusion;
  for (const double val : rdfDensities) {
    file << "," << val;
  }
  file << "\n";

  // File closes automatically
}
