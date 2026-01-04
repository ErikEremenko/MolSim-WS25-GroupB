#pragma once

#include <string>
#include <optional>

#include "Simulation.h"

enum class LogLevelConfig { OFF, ERROR, WARN, INFO, DEBUG, TRACE };

struct RunConfig {
  std::string filename;
  SimulationMode simulationMode;
  LogLevelConfig logLevel;

  bool useParallelization;
  bool isYaml;

  std::optional<double> tEnd;
  std::optional<double> deltaT;
  std::optional<YAMLSimulation::ContainerKind> containerKind;
};

class ArgParser {
private:
  static void printUsage();
  static LogLevelConfig parseLogLevel(const std::string& logLevelStr);
  static SimulationMode parseSimulationMode(const std::string& simModeStr);
public:
  static std::optional<RunConfig> parseArgs(int argc, char* argv[]);
  static void setLogLevel(LogLevelConfig logLevel);
};
