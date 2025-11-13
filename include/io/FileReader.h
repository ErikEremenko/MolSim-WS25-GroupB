#pragma once

#include "ParticleContainer.h"

class BaseFileReader {  // TODO: Implement YAML reader next week
protected:
  std::string filename;
public:
  explicit BaseFileReader(std::string filename);
  virtual ~BaseFileReader();

  virtual void readFile(ParticleContainer& particles);
};

class FileCuboidReader : BaseFileReader {
public:
  using BaseFileReader::BaseFileReader;

  void readFile(ParticleContainer& particles) override;
};
