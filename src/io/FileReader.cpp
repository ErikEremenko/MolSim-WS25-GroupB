#include "io/FileReader.h"
#include "utils/MaxwellBoltzmannDistribution.h"  // include for testing (? TODO)

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

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
    SPDLOG_DEBUG("Read line: {}", tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      SPDLOG_DEBUG("Read line: {}", tmp_string);
    }

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
    SPDLOG_DEBUG("Read line: {}", tmp_string);

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      SPDLOG_DEBUG("Read line: {}", tmp_string);
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_cuboids;
    SPDLOG_DEBUG("Reading {} cuboids.", num_cuboids);
    SPDLOG_DEBUG("Reading {} cuboids.", num_cuboids);
    getline(input_file, tmp_string);
    SPDLOG_DEBUG("Read line: {}", tmp_string);

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

      //code for generating particles:
      for (int nx = 0; nx < n[0]; nx++) {
        for (int ny = 0; ny < n[1]; ny++) {
          for (int nz = 0; nz < n[2]; nz++) {

            std::array<double, 3> tempx = {cx[0] + h * nx, cx[1] + h * ny, cx[2] + h * nz};
            std::array<double, 3> tempv = {cv[0], cv[1], cv[2]};
            std::array<double, 3> temperatureVel = maxwellBoltzmannDistributedVelocity(t, 2);
            tempv[0] += temperatureVel[0];
            tempv[1] += temperatureVel[1];
            tempv[2] += temperatureVel[2];
            particles.addParticle(tempx, tempv, m);
          }
        }
      }

      getline(input_file, tmp_string);
      SPDLOG_DEBUG("Read line: {}", tmp_string);
    }
  } else {
    SPDLOG_ERROR("Error: could not open file {}", filename);
    exit(-1);
  }
}
