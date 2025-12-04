#include "io/FileReader.h"
#include "utils/MaxwellBoltzmannDistribution.h"  // include for testing (? TODO)
#include "../include/ParticleGenerator.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#endif  // SPDLOG_ACTIVE_LEVEL
#include <spdlog/spdlog.h>

BaseFileReader::BaseFileReader(std::string filename) : filename(std::move(filename)) {}

BaseFileReader::~BaseFileReader() = default;

void BaseFileReader::readFile(ParticleContainer& particles) {
  std::array<double, 3> x{};
  std::array<double, 3> v{};
  double m;
  int num_particles = 0;

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {
    getline(input_file, tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
    }
    SPDLOG_DEBUG("Read line: {}", tmp_string);

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    SPDLOG_DEBUG("Reading {}.", tmp_string);
    getline(input_file, tmp_string);
    SPDLOG_DEBUG("Read line: {}", tmp_string);

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto& xj : x) {
        datastream >> xj;
      }
      for (auto& vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        SPDLOG_ERROR("Error reading file: eof reached unexpectedly reading from line {}", i);
        exit(-1);
      }
      datastream >> m;
      particles.addParticle(x, v, m);
      //      SPDLOG_DEBUG("Generated particle with x={{{0:.2f}, {0:.2f}, {0:.2f}}}, v={{{0:.2f}, {0:.2f}, {0:.2f}}}, m={0:.2f}",
      //                   x[0], x[1], x[2], v[0], v[1], v[2], m);

      getline(input_file, tmp_string);
      SPDLOG_DEBUG("Read line: {}", tmp_string);
    }
  } else {
    SPDLOG_ERROR("Error: could not open file {}", filename);
    exit(-1);
  }
}

// CuboidFileReader class definition
void CuboidFileReader::readFile(ParticleContainer& particles) {
  std::array<double, 3> cx{};
  std::array<double, 3> cv{};
  std::array<int, 3> n{};
  double h;
  double m;
  double t;

  int num_cuboids = 0;

  std::ifstream input_file(filename);
  std::string tmp_string;

  if (input_file.is_open()) {
    getline(input_file, tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
    }
    SPDLOG_DEBUG("Read line: {}", tmp_string);

    std::istringstream numstream(tmp_string);
    numstream >> num_cuboids;

    SPDLOG_DEBUG("Reading {} cuboids.", num_cuboids);
    getline(input_file, tmp_string);
    SPDLOG_DEBUG("Read line: {}", tmp_string);

    ParticleGenerator particleGenerator(particles);

    for (int i = 0; i < num_cuboids; i++) {
      std::istringstream datastream(tmp_string);

      for (auto& cxj : cx) {
        datastream >> cxj;
      }
      for (auto& cvj : cv) {
        datastream >> cvj;
      }
      for (auto& nj : n) {
        datastream >> nj;
      }
      if (datastream.eof()) {
        SPDLOG_ERROR("Error reading file: eof reached unexpectedly reading from line {}", i);
        exit(-1);
      }
      datastream >> h;
      datastream >> m;
      datastream >> t;

      particleGenerator.generateCuboid(cx, cv, n, h, m, t);

      getline(input_file, tmp_string);
      SPDLOG_DEBUG("Read line: {}", tmp_string);
      SPDLOG_DEBUG(
          "Generated cuboid of {} particles with position {{{}, {}, {}}}, velocity {{{}, {}, {}}}, particle separation "
          "{} and average brownian motion velocity {}",
          n[0] * n[1] * n[2], cx[0], cx[1], cx[2], cv[0], cv[1], cv[2], h, t);
    }
  } else {
    SPDLOG_ERROR("Error: could not open file {}", filename);
    exit(-1);
  }
}
