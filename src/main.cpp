#include "Simulation.h"

#include <iostream>

#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO  // TODO: Make this a global define using CMake
#include "spdlog/spdlog.h"

int main(const int argc, char* argsv[]) {
  SPDLOG_INFO("Hello from MolSim for PSE!");
  // Read arguments from the command line
  if (argc != 6) {
    SPDLOG_ERROR("Erroneous programme call!");
    SPDLOG_ERROR("./MolSim filename t_end delta_t [file | benchmark]");
    return 1;
  }

  SimulationMode simulation_mode;
  std::string mode_arg = argsv[4];
  if (mode_arg == "file") {
    simulation_mode = SimulationMode::FILE_OUTPUT;
  } else if (mode_arg == "benchmark") {
    simulation_mode = SimulationMode::BENCHMARK;
  } else {
    SPDLOG_ERROR("Invalid simulation mode provided: {}", mode_arg);
    SPDLOG_ERROR("Please use 'file' or 'benchmark'");
    return 1;
  }

  std::string log_level = argsv[5];

  if (log_level == "info"){
    spdlog::set_level(spdlog::level::info);
  }else if(log_level == "off"){
    spdlog::set_level(spdlog::level::off);
  }else if(log_level == "trace"){
    spdlog::set_level(spdlog::level::trace);
  }else if(log_level == "error"){
    spdlog::set_level(spdlog::level::err);
  }else if(log_level == "debug"){
    spdlog::set_level(spdlog::level::debug);
  }else{
    return 1;
  }

  CollisionSimulation simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), simulation_mode);
  simulation.run();

  return 0;
}
