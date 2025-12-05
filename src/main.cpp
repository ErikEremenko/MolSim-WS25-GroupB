#include "Simulation.h"

#include <iostream>

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG  // TODO: Make this a global define using CMake
#endif                                          // SPDLOG_ACTIVE_LEVEL
#include "spdlog/spdlog.h"

int main(const int argc, char* argsv[]) {
  SPDLOG_INFO("Hello from MolSim for PSE!");
  // Read arguments from the command line
  if (argc < 2) {
    SPDLOG_ERROR("Erroneous programme call!");
    SPDLOG_ERROR(
        "YAML mode: ./MolSim filename [file | benchmark] [off | error | debug | trace | info] [linked | direct]");
    SPDLOG_ERROR(
        "Legacy mode: ./MolSim filename t_end delta_t [file | benchmark] [off | error | debug | trace | info] [P:OFF | "
        "P:ON]");
    return 1;
  }

  std::string filename = argsv[1];

  SimulationMode simulation_mode = SimulationMode::FILE_OUTPUT;
  int mode_arg_index = 2;  // for .yaml

  bool is_yaml = filename.find(".yaml") != std::string::npos || filename.find(".yml") != std::string::npos;

  if (!is_yaml && argc != 7) {
    SPDLOG_ERROR("Erroneous programme call!");
    SPDLOG_ERROR("Legacy .txt mode requires 7 arguments!");
    SPDLOG_ERROR(
        "./MolSim filename t_end delta_t [file | benchmark] [off | error | debug | trace | info] [P:OFF | "
        "P:ON]");
    return 1;
  }

  if (!is_yaml)
    mode_arg_index = 4;  // for legacy .txt support

  if (argc > mode_arg_index) {
    if (std::string mode_arg = argsv[mode_arg_index]; mode_arg == "file") {
      simulation_mode = SimulationMode::FILE_OUTPUT;
    } else if (mode_arg == "benchmark") {
      simulation_mode = SimulationMode::BENCHMARK;
    } else {
      SPDLOG_ERROR("Invalid simulation mode provided: {}", mode_arg);
      SPDLOG_ERROR("Please use 'file' or 'benchmark'");
      return 1;
    }
  }
  if (int log_arg_index = is_yaml ? 3 : 5; argc > log_arg_index) {
    std::string log_level = argsv[log_arg_index];
    if (log_level == "info") {
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
  }
  try {
    if (is_yaml) {
      YAMLSimulation::ContainerKind kind = YAMLSimulation::ContainerKind::LINKED;
      if (argc > 4) {
        std::string k = argsv[4];
        if (k == "direct")
          kind = YAMLSimulation::ContainerKind::DIRECT;
        else if (k == "linked")
          kind = YAMLSimulation::ContainerKind::LINKED;
        else {
          SPDLOG_ERROR("Invalid container kind: {}. Use 'direct' or 'linked'.", k);
          return 1;
        }
      }

      YAMLSimulation::Parallelization parallel = YAMLSimulation::Parallelization::OFF;
      if (argc > 5) {
        std::string p = argsv[5];
        if (p == "P:ON")
          parallel = YAMLSimulation::Parallelization::ON;
        else if (p == "P:OFF")
          parallel = YAMLSimulation::Parallelization::OFF;
        else {
          SPDLOG_ERROR("Invalid parallelization option: {}. Use 'P:ON' or 'P:OFF'.", p);
          return 1;
        }
      }

      YAMLSimulation simulation(argsv[1], simulation_mode, kind, parallel);
      simulation.run();
    } else {
      // legacy simulation
      if (std::string parallelization = argsv[6]; parallelization == "P:OFF") {
        CollisionSimulation simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), simulation_mode);
        simulation.run();
      } else if (parallelization == "P:ON") {
        CollisionSimulationParallel simulation(argsv[1], std::stod(argsv[2]), std::stod(argsv[3]), simulation_mode);
        simulation.run();
      } else {
        SPDLOG_ERROR(
            "Invalid parallelization option of file name. Valid options are P:OFF and P:ON. Input filename must not be "
            "empty");
      }
    }

  } catch (const std::invalid_argument& e) {
    SPDLOG_ERROR("Simulation failed: {}", e.what());
    return 1;
  }

  return 0;
}
