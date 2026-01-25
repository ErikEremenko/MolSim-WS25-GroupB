#include <iostream>

#include "io/ArgParser.h"
#include "io/YAMLFileReader.h"
#include "simulation/Simulation.h"
#include "simulation/SimulationConfig.h"

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL \
  SPDLOG_LEVEL_DEBUG  // TODO: Make this a global define using CMake or even remove completely?
#endif                // SPDLOG_ACTIVE_LEVEL
#include <spdlog/spdlog.h>

int main(const int argc, char* argv[]) {
  // Greeting
  spdlog::set_level(spdlog::level::info);
  SPDLOG_INFO("Hello from MolSim for PSE!");

  // --- Argument parsing and logger configuration ---
  ArgParser argParser(argc, argv);
  std::optional<CLIConfig> cliConfig = argParser.parse();
  if (!cliConfig) {
    return 1;  // CLI parsing failed
  }

  const CLIConfig& cli = cliConfig.value();
  if (cli.logLevel) {
    ArgParser::setLogLevel(cli.logLevel.value());
  }

  SimulationConfig simConfig;
  if (cli.isYaml) {  // TODO: If we want to add a cuboid file reader, we can introduce a file type enum
    // YAML mode
    SPDLOG_INFO("YAML mode, Reading file: {}", cli.filename);
    YAMLFileReader fileReader(cli.filename);
    simConfig = fileReader.getConfig();

    // Check possible CLI overrides (CLI > YAML)
    if (cli.simulationMode) {
      simConfig.simulationMode = *cli.simulationMode;
      SPDLOG_WARN("CLI Override: Simulation Mode");
    }
    if (cli.containerType) {
      simConfig.containerType = *cli.containerType;
      SPDLOG_WARN("CLI Override: Container Type");
    }
    if (cli.useParallelization) {
      simConfig.useParallelization = *cli.useParallelization;
      SPDLOG_WARN("CLI Override: Parallelization");
    }
  } else {
    // Legacy mode
    SPDLOG_INFO("Legacy mode, configuring from CLI arguments.");

    // Mandatory legacy arguments
    simConfig.tEnd = cli.tEnd.value();
    simConfig.deltaT = cli.deltaT.value();
    simConfig.simulationMode = cli.simulationMode.value();

    // Parallelization
    simConfig.useParallelization = cli.useParallelization.value();

    // Container Type, legacy mode defaults to LINKED
    simConfig.containerType = ContainerType::LINKED;

    // Load particles (cuboids from file) - TODO: This is not implemented yet
    SPDLOG_WARN("Legacy mode particle loading not yet implemented in refactor.");
  }

  // Running the simulation
  Simulation simulation(simConfig);
  try {
    simulation.run();
  } catch (const std::exception& e) {
    SPDLOG_ERROR("Simulation runtime error: {}", e.what());
    return 1;
  }

  // Ending message
  SPDLOG_INFO("Simulation ended.");

  return 0;
}
