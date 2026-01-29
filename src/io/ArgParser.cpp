#include "io/ArgParser.h"

#include <spdlog/spdlog.h>

ArgParser::ArgParser(const int argc, char* argv[]) : argc(argc), args(argv, argv + argc) {}

ArgParser::~ArgParser() = default;

void ArgParser::printUsage() {
  SPDLOG_INFO("Usage:");
  SPDLOG_INFO(
    "  YAML mode: ./MolSim filename [file | benchmark] [off | error | debug | trace | info] [linked | direct] [P:ON | P:OFF] [S:COLORING | S:TASKBASED]");
  SPDLOG_INFO(
    "  Legacy mode: ./MolSim filename t_end delta_t [file | benchmark] [off | error | debug | trace | info] [P:OFF | "
    "P:ON] [S:COLORING | S:TASKBASED]");
}

LogLevelConfig ArgParser::parseLogLevel(const std::string& logLevelStr) {
  // TODO: Replace this with compiler flags
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

ContainerType ArgParser::parseContainerType(const std::string& containerTypeStr) {
  // TODO: This function has a duplicate in YAMLFileReader
  // TODO: This function should be removed anyway, as specifying the container type in YAML is more natural
  if (containerTypeStr == "direct")
    return ContainerType::DIRECT;
  if (containerTypeStr == "linked")
    return ContainerType::LINKED;
  throw std::invalid_argument("Invalid container type: " + containerTypeStr);
}

bool ArgParser::parseParallelization(const std::string& parallelStr) {
  if (parallelStr == "P:ON")
    return true;
  if (parallelStr == "P:OFF")
    return false;
  throw std::invalid_argument("Invalid parallel option: " + parallelStr);
}

ParallelStrategy ArgParser::parseParallelStrategy(const std::string& strategyStr) {
  if (strategyStr == "S:COLORING" || strategyStr == "S:coloring")
    return ParallelStrategy::COLORING;
  if (strategyStr == "S:TASKBASED" || strategyStr == "S:taskbased")
    return ParallelStrategy::TASKBASED;
  throw std::invalid_argument("Invalid parallel strategy: " + strategyStr);
}

std::optional<CLIConfig> ArgParser::parse() const {
  if (argc < 2) {
    printUsage();
    return std::nullopt;
  }

  CLIConfig config;
  config.filename = args[1];

  // Detect YAML vs Legacy based on extension
  config.isYaml =
      (config.filename.find(".yaml") != std::string::npos || config.filename.find(".yml") != std::string::npos);

  try {
    if (config.isYaml) {
      // YAML Mode, sparse parsing
      // We only set values if arguments are present

      // Index 2: Mode (Optional)
      if (args.size() > 2)
        config.simulationMode = parseSimulationMode(args[2]);

      // Index 3: Log Level (Optional)
      if (args.size() > 3)
        config.logLevel = parseLogLevel(args[3]);

      /* TODO: Remove container and parallel options from the CLI
      * - Parallel should be moved to a compiler flag
      * - Container should be removed completely and only specified in YAML
      */
      // Index 4: Container (Optional)
      if (args.size() > 4)
        config.containerType = parseContainerType(args[4]);

      // Index 5: Parallel (Optional)
      if (args.size() > 5)
        config.useParallelization = parseParallelization(args[5]);

      // Index 6: Strategy (Optional)
      if (args.size() > 6)
        config.parallelStrategy = parseParallelStrategy(args[6]);

    } else {
      // Legacy mode, full parsing
      // Requires 7 or 8 arguments (program + 6 args + optional strategy)
      if (argc != 7 && argc != 8) {
        SPDLOG_ERROR("Legacy mode requires 7 arguments (or 8 with strategy).");
        printUsage();
        return std::nullopt;
      }

      // Default container for legacy is direct
      config.tEnd = std::stod(args[2]);
      config.deltaT = std::stod(args[3]);
      config.simulationMode = parseSimulationMode(args[4]);
      config.logLevel = parseLogLevel(args[5]);
      config.useParallelization = parseParallelization(args[6]);

      if (argc == 8)
        config.parallelStrategy = parseParallelStrategy(args[7]);
    }
  } catch (const std::exception& e) {
    SPDLOG_ERROR("Argument parsing error: {}", e.what());
    printUsage();
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
