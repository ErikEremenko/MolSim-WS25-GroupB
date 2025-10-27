/*
 * XYZWriter.h
 *
 *  Created on: 01.03.2010
 *      Author: eckhardw
 */

#pragma once

#include <fstream>
#include <vector>

#include "Particle.h"

namespace outputWriter {

class XYZWriter {
 public:
  XYZWriter();

  virtual ~XYZWriter();

  static void plotParticles(const std::vector<Particle>& particles, const std::string &filename, int iteration);
};

}  // namespace outputWriter
