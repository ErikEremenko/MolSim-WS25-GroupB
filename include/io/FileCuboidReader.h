//
// Created by User on 08/11/2025.
//

#pragma once

#include "ParticleContainer.h"

class FileCuboidReader {
 public:
  FileCuboidReader();
  virtual ~FileCuboidReader();

  static void readFile(ParticleContainer& particles, char* filename);
};
