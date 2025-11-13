#include "io/FileReader.h"
#include "utils/MaxwellBoltzmannDistribution.h"  // include for testing (? TODO)

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

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
    std::cout << "Read line: " << tmp_string << std::endl;

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_particles;
    std::cout << "Reading " << num_particles << "." << std::endl;
    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << std::endl;

    for (int i = 0; i < num_particles; i++) {
      std::istringstream datastream(tmp_string);

      for (auto& xj : x) {
        datastream >> xj;
      }
      for (auto& vj : v) {
        datastream >> vj;
      }
      if (datastream.eof()) {
        std::cout << "Error reading file: eof reached unexpectedly reading from line " << i << std::endl;
        exit(-1);
      }
      datastream >> m;
      particles.addParticle(x, v, m);

      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }
  } else {
    std::cout << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }
}

// FileCuboidReader class definition
void FileCuboidReader::readFile(ParticleContainer& particles) {
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
    std::cout << "Read line: " << tmp_string << std::endl;

    while (tmp_string.empty() or tmp_string[0] == '#') {
      getline(input_file, tmp_string);
      std::cout << "Read line: " << tmp_string << std::endl;
    }

    std::istringstream numstream(tmp_string);
    numstream >> num_cuboids;
    std::cout << "Reading " << num_cuboids << " cuboids" << "." << std::endl;
    getline(input_file, tmp_string);
    std::cout << "Read line: " << tmp_string << std::endl;

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
        std::cout << "Error reading file: eof reached unexpectedly reading from line " << i << std::endl;
        exit(-1);
      }
      datastream >> h;
      datastream >> m;
      datastream >> t;

      //particles.addParticle(x, v, m);
      //generateParticleCuboid(cx, cv, n, h, m);

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
      std::cout << "Read line: " << tmp_string << std::endl;
    }
  } else {
    std::cout << "Error: could not open file " << filename << std::endl;
    exit(-1);
  }
}
