#pragma once

#include <optional>
#include <string>

#include "simulation/SimulationConfig.h"

/**
 * @enum LogLevelConfig
 * @brief Used to set the logging level of the spdlog-logger.
 */
enum class LogLevelConfig { OFF, ERROR, WARN, INFO, DEBUG, TRACE };

/**
 * @struct CLIConfig
 * @brief Holds strictly what was parsed from the command line.
 */
struct CLIConfig {
  std::string filename;
  bool isYaml;
  std::optional<LogLevelConfig> logLevel;

  // CLI overrides the parameters below even if they are defined in YAML
  std::optional<SimulationMode> simulationMode;
  std::optional<ContainerType> containerType;
  std::optional<bool> useParallelization;

  // Legacy Parameters
  std::optional<double> tEnd;
  std::optional<double> deltaT;
};

constexpr int FILENAME_INDEX = 1;

class ArgParser {
 private:
  /**
   * @brief Number of command line arguments passed to the program.
   */
  int argc;

  /**
   * @brief Vector containing the command line arguments as strings.
   */
  std::vector<std::string> args;

  // Helper functions below
  /**
   * @brief Prints the correct usage syntax.
   * Used when argument validation fails.
   */
  static void printUsage();

  /**
   * @brief Parses a string into a LogLevelConfig enum.
   * @param logLevelStr The string representation of the log level.
   * @return The corresponding LogLevelConfig enum value.
   */
  static LogLevelConfig parseLogLevel(const std::string& logLevelStr);

  /**
   * @brief Parses a string into a SimulationMode enum.
   * @param simModeStr The string representation of the simulation mode (file or benchmark).
   * @return The corresponding SimulationMode enum value.
   */
  static SimulationMode parseSimulationMode(const std::string& simModeStr);

  /**
   * @brief Parses a string into a ContainerType enum.
   * @param containerTypeStr The string representation of the container strategy (direct or linked).
   * @return The corresponding ContainerType enum value.
   */
  static ContainerType parseContainerType(const std::string& containerTypeStr);

  /**
   * @brief Parses a string into a boolean flag for parallelization.
   * @param parallelStr The string representation of the parallel option ("P:ON" or "P:OFF").
   * @return true if P:ON, false if P:OFF.
   */
  static bool parseParallelization(const std::string& parallelStr);

 public:
  ArgParser(int argc, char* argv[]);
  ~ArgParser();

  /**
   * @brief Parses the command line arguments.
   * @return A config struct if parsing was successful, otherwise an empty option.
   */
  std::optional<CLIConfig> parse() const;

  /**
   * @brief Sets the logging level of the spdlog-logger.
   * @param logLevel Logging level to be used by spdlog
   */
  static void setLogLevel(LogLevelConfig logLevel);
};
