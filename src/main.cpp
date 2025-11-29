#include "Simulation.h"

#include <iostream>

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG  // TODO: Make this a global define using CMake
#endif                                          // SPDLOG_ACTIVE_LEVEL
#include "spdlog/spdlog.h"

int main(const int argc, char* argsv[]) {
  SPDLOG_INFO("Hello from MolSim for PSE!");
  // Read arguments from the command line
  if (argc != 7) {
    SPDLOG_ERROR("Erroneous programme call!");
    SPDLOG_ERROR(
        "./MolSim filename t_end delta_t [file | benchmark] [off | error | debug | trace | info] [P:OFF | "
        "P:ON]");
    return 1;
  }

  SimulationMode simulation_mode;
  if (std::string mode_arg = argsv[4]; mode_arg == "file") {
    simulation_mode = SimulationMode::FILE_OUTPUT;
  } else if (mode_arg == "benchmark") {
    simulation_mode = SimulationMode::BENCHMARK;
  } else {
    SPDLOG_ERROR("Invalid simulation mode provided: {}", mode_arg);
    SPDLOG_ERROR("Please use 'file' or 'benchmark'");
    return 1;
  }

  if (std::string log_level = argsv[5]; log_level == "info") {
    spdlog::set_level(spdlog::level::info);
  } else if (log_level == "off") {
    spdlog::set_level(spdlog::level::off);
  } else if (log_level == "error") {
    spdlog::set_level(spdlog::level::err);
  } else if (log_level == "warn") {
    spdlog::set_level(spdlog::level::warn);
  } else if (log_level == "debug") {
    spdlog::set_level(spdlog::level::debug);
  } else if (log_level == "trace") {
    spdlog::set_level(spdlog::level::trace);
  } else {
    SPDLOG_ERROR("Invalid log level. Valid options are off, error, warn, info, debug and trace.");
    return 1;
  }

  if (std::string parallelization = argsv[6]; parallelization == "P:OFF") {
    CollisionSimulation simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), simulation_mode);
    simulation.run();
  } else if (parallelization == "P:ON") {
    CollisionSimulationParallel simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), simulation_mode);
    simulation.run();
  } else {
    SPDLOG_ERROR("Invalid parallelization option. Valid options are P:OFF and P:ON.");
  }

  return 0;
}
