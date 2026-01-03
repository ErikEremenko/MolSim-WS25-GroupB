#include <iostream>

#include "io/ArgParser.h"
#include "Simulation.h"

#ifndef SPDLOG_ACTIVE_LEVEL
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG  // TODO: Make this a global define using CMake or even remove completely?
#endif                                          // SPDLOG_ACTIVE_LEVEL
#include "spdlog/spdlog.h"

int main(const int argc, char* argv[]) {
  // Greeting
  spdlog::set_level(spdlog::level::info);
  SPDLOG_INFO("Hello from MolSim for PSE!");

  // Argument parsing and logger configuration
  const std::optional<RunConfig> configOpt = ArgParser::parseArgs(argc, argv);
  if (!configOpt) return 1;
  const RunConfig& config = configOpt.value();
  ArgParser::setLogLevel(config.logLevel);

  // Running the simulation
  try {
    if (config.isYaml) {
      auto parallelMode = config.useParallelization ?
                          YAMLSimulation::Parallelization::ON :
                          YAMLSimulation::Parallelization::OFF;

      // We can safely dereference *config.containerKind because the parser guarantees it exists in YAML mode
      YAMLSimulation simulation(
          config.filename,
          config.simulationMode,
          *config.containerKind,
          parallelMode
      );
      simulation.run();
    }
    else {  // Legacy Mode
      // The parser guarantees tEnd and deltaT exist if isYamlMode is false
      if (config.useParallelization) {
        CollisionSimulationParallel simulation(
            config.filename,
            config.tEnd.value(),
            config.deltaT.value(),
            config.simulationMode
        );
        simulation.run();
      } else {
        CollisionSimulation simulation(
            config.filename,
            *config.tEnd,
            *config.deltaT,
            config.simulationMode
        );
        simulation.run();
      }
    }
  } catch (const std::exception& e) {
    SPDLOG_ERROR("Simulation runtime error: {}", e.what());
    return 1;
  }

  return 0;
}
