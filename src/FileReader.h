/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#pragma once

#include <vector>

#include "Particle.h"

class FileReader {
 public:
  FileReader();
  virtual ~FileReader();

  static void readFile(std::vector<Particle>& particles, char* filename);
};
