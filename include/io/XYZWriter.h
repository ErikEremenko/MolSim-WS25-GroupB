/**
 * @file XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include <fstream>
#include <vector>

#include "Particle.h"

namespace outputWriter {

/**
 * @brief A singleton class that can be used to describe the state of the simulation in a given moment.
 *
 */
class XYZWriter final {
 public:
  XYZWriter();

  virtual ~XYZWriter();
  /**
   * @brief creates a file, which describes the state of all particles in the system,
   * for a particular iteration.
   * @param particles A vector holding all particles in the system.
   * @param filename A reference to the string that represents the name of the file being generated.
   * @param iteration Integer number representing the number of the iteration being plotted.
   */
  static void plotParticles(const std::vector<Particle>& particles,
                            const std::string& filename, int iteration);
};

}  // namespace outputWriter
