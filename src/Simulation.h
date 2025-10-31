#pragma once

#include "CalcMethod.h"
#include "FileReader.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XYZWriter.h"

class Simulation {
 private:
  char* filename;
  const double end_time;
  const double dt;

  ParticleContainer& particles;

  CalcMethod& calcMethod;

  void plotParticles(const int iteration) const {
    const std::string out_name("MD_vtk");
    /*outputWriter::XYZWriter writer;
            outputWriter::XYZWriter::plotParticles(particles, out_name, iteration);*/
    outputWriter::VTKWriter writer;
    writer.plotParticles(particles, out_name, iteration);
  }

 public:
  Simulation(char* filename, const double end_time, const double dt,
             ParticleContainer& particles, CalcMethod& calcMethod)
      : filename(filename),
        end_time(end_time),
        dt(dt),
        particles(particles),
        calcMethod(calcMethod) {}
  ~Simulation() = default;

  void loadParticles() const {
    FileReader fileReader;
    fileReader.readFile(particles, filename);
  }

  void run() const {
    constexpr double start_time = 0;

    // Load particles from file

    double current_time = start_time;
    int iteration = 0;

    // for this loop, we assume: current x, current f and current v are known
    while (current_time < end_time) {
      // calculate new x
      calcMethod.calculateX(dt);
      for (auto& p : particles) {
        p.setOldF(p.getF());  // store f(t_n) for v update
      }
      // calculate new f
      calcMethod.calculateF(dt);
      // calculate new v
      calcMethod.calculateV(dt);

      iteration++;
      if (iteration % 10 == 0) {
        plotParticles(iteration);
      }
      // std::cout << "Iteration " << iteration << " finished." << std::endl;

      current_time += dt;
    }
  }
};
