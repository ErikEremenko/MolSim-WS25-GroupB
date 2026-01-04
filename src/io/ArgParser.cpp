#include "io/ArgParser.h"

#include <spdlog/spdlog.h>

void ArgParser::printUsage() {
  SPDLOG_ERROR("Usage:");
  SPDLOG_ERROR(
      "  YAML mode: ./MolSim filename [file | benchmark] [off | error | debug | trace | info] [linked | direct]");
  SPDLOG_ERROR(
      "  Legacy mode: ./MolSim filename t_end delta_t [file | benchmark] [off | error | debug | trace | info] [P:OFF | "
      "P:ON]");
}

LogLevelConfig ArgParser::parseLogLevel(const std::string& logLevelStr) {
  if (logLevelStr == "off")
    return LogLevelConfig::OFF;
  if (logLevelStr == "error")
    return LogLevelConfig::ERROR;
  if (logLevelStr == "warn")
    return LogLevelConfig::WARN;
  if (logLevelStr == "info")
    return LogLevelConfig::INFO;
  if (logLevelStr == "debug")
    return LogLevelConfig::DEBUG;
  if (logLevelStr == "trace")
    return LogLevelConfig::TRACE;
  throw std::invalid_argument("Invalid log level: " + logLevelStr);
}

SimulationMode ArgParser::parseSimulationMode(const std::string& simModeStr) {
  if (simModeStr == "file")
    return SimulationMode::FILE_OUTPUT;
  if (simModeStr == "benchmark")
    return SimulationMode::BENCHMARK;
  throw std::invalid_argument("Invalid simulation mode: " + simModeStr);
}

std::optional<RunConfig> ArgParser::parseArgs(int argc, char* argv[]) {
  if (argc < 2) {
    printUsage();
    return std::nullopt;
  }

  std::vector<std::string> args(argv, argv + argc);  // cpp-safe version
  RunConfig config;

  // Set defaults
  config.filename = args[1];
  config.simulationMode = SimulationMode::FILE_OUTPUT;
  config.logLevel = LogLevelConfig::INFO;
  config.useParallelization = false;

  // Detect YAML extension
  config.isYaml =
      (config.filename.find(".yaml") != std::string::npos || config.filename.find(".yml") != std::string::npos);

  try {
    if (config.isYaml) {
      // Indices: 0-prog, 1-file, 2-mode, 3-log, 4-container, 5-parallel

      if (args.size() > 2)
        config.simulationMode = parseSimulationMode(args[2]);
      if (args.size() > 3)
        config.logLevel = parseLogLevel(args[3]);

      config.containerKind = YAMLSimulation::ContainerKind::LINKED;  // Default
      if (args.size() > 4) {
        if (args[4] == "direct")
          config.containerKind = YAMLSimulation::ContainerKind::DIRECT;
        else if (args[4] == "linked")
          config.containerKind = YAMLSimulation::ContainerKind::LINKED;
        else
          throw std::invalid_argument("Invalid container kind: " + args[4]);
      }

      if (args.size() > 5) {
        if (args[5] == "P:ON")
          config.useParallelization = true;
        else if (args[5] == "P:OFF")
          config.useParallelization = false;
        else
          throw std::invalid_argument("Invalid parallel option: " + args[5]);
      }
    } else {
      // Legacy mode parsing
      if (argc != 7) {
        SPDLOG_ERROR("Legacy mode requires exactly 7 arguments.");
        printUsage();
        return std::nullopt;
      }

      // Indices: 0:prog, 1:file, 2:t_end, 3:dt, 4:mode, 5:log, 6:parallel
      config.tEnd = std::stod(args[2]);
      config.deltaT = std::stod(args[3]);
      config.simulationMode = parseSimulationMode(args[4]);
      config.logLevel = parseLogLevel(args[5]);

      if (args[6] == "P:ON")
        config.useParallelization = true;
      else if (args[6] == "P:OFF")
        config.useParallelization = false;
      else
        throw std::invalid_argument("Invalid parallel option: " + args[6]);
    }
  } catch (const std::exception& e) {
    SPDLOG_ERROR("Argument parsing error: {}", e.what());
    return std::nullopt;
  }

  return config;
}

void ArgParser::setLogLevel(LogLevelConfig logLevel) {
  switch (logLevel) {
    case LogLevelConfig::OFF:
      spdlog::set_level(spdlog::level::off);
      break;
    case LogLevelConfig::ERROR:
      spdlog::set_level(spdlog::level::err);
      break;
    case LogLevelConfig::WARN:
      spdlog::set_level(spdlog::level::warn);
      break;
    case LogLevelConfig::INFO:
      spdlog::set_level(spdlog::level::info);
      break;
    case LogLevelConfig::DEBUG:
      spdlog::set_level(spdlog::level::debug);
      break;
    case LogLevelConfig::TRACE:
      spdlog::set_level(spdlog::level::trace);
      break;
  }
}
