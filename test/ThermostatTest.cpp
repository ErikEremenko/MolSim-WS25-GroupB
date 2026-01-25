#include <functional>
#include <memory>

#include <gtest/gtest.h>

#include "../include/simulation/Simulation.h"
#include "physics/ParticleGenerator.h"

// Simulation defines
constexpr double SIM_END_TIME = 0.2;
constexpr double SIM_DT = 0.002;      // ~100 simulation iterations
constexpr int THERMO_FREQUENCY = 20;  // ~5 thermostat applications

constexpr double SIM_FORCE_EPSILON = 5.0;
constexpr double SIM_FORCE_SIGMA = 1.0;

constexpr double SIM_FORCE_CUTOFF_RADIUS = 3.0;

// Temperature test defines
constexpr double THERMO_TEMP_TOLERANCE = 0.5;  // how far off is the thermostat allowed to be?

constexpr double THERMO_CUBOID_TEMP = 27.0;

// Temperature holding test defines
constexpr double THERMO_HOLDING_TEMP = 50.0;

// Heating test defines
constexpr double THERMO_HEATING_START_TEMP = 10;
constexpr double THERMO_HEATING_TARGET_TEMP = 60;
constexpr double THERMO_HEATING_TEMP_DELTA = 5;

// Cooling test defines
constexpr double THERMO_COOLING_START_TEMP = 150;
constexpr double THERMO_COOLING_TARGET_TEMP = 50;
constexpr double THERMO_COOLING_TEMP_DELTA = 10;

class ThermostatTestingSimulation : public Simulation {
 private:
  Thermostat& thermostat;
  std::function<void(bool)> temperatureChecker;  // called after updating the temperature

 protected:
  void setupSimulation() override { /* empty override for compilation */ }
  void runFileOutput() override { /* empty as it will not be called */ }
  void runBenchmark() override {
    // No benchmarking, just simulate and regularly apply thermostat
    double current_time = 0;
    int iteration = 0;
    const int thermostatFrequency = thermostat.getUpdateFrequency();
    while (current_time < endTime) {
      ForceCalc::calculateX(*particles, dt);
      for (auto& p : *particles) {
        p.setOldF(p.getF());
        p.setF({});
      }
      for (const auto& force : forces) {
        force->calculateF();
      }
      ForceCalc::calculateV(*particles, dt);

      iteration++;
      if (iteration % thermostatFrequency == 0) {
        temperatureChecker(true);  // before update
        thermostat.updateTemperature();
        temperatureChecker(false);  // after update
      }
      current_time += dt;
    }
  }

 public:
  ThermostatTestingSimulation(std::unique_ptr<ParticleContainer> particles, Thermostat& thermostat,
                              std::function<void(bool)> temperatureChecker, SimulationConfig& config)
      : Simulation(config), thermostat(thermostat), temperatureChecker(std::move(temperatureChecker)) {

    this->particles = std::move(particles);
    this->forces.push_back(std::make_unique<LennardJonesForce>(*this->particles, SIM_FORCE_EPSILON, SIM_FORCE_SIGMA,
                                                               SIM_FORCE_CUTOFF_RADIUS));
  }
  ~ThermostatTestingSimulation() override = default;
};

class ThermostatTest : public testing::Test {
 protected:
  std::unique_ptr<ParticleContainer> particles;  // ownership transferred to simulation later on
  ParticleGenerator generator;

  ThermostatTest() : particles(std::make_unique<ParticleContainer>()) {}
  ~ThermostatTest() override = default;
};

/**
 * @brief Check that the thermostat initializes the temperature correctly.
 */
TEST_F(ThermostatTest, CheckThermostatInitialTemperature) {
  // TODO
}

/**
 * @brief Check that the thermostat correctly calculates the system's current temperature from its kinetic energy.
 */
TEST_F(ThermostatTest, CheckThermostatTemperatureCalculation) {
  // Set up 5 particles
  constexpr double sigma = 1;
  constexpr double epsilon = 1;
  particles->addParticle({0, 0, 0}, {1, 0, 0}, 2, sigma, epsilon);    // m*v^2 = 2*1*1 = 2
  particles->addParticle({1, 1, 1}, {0, 2, 0}, 1, sigma, epsilon);    // m*v^2 = 1*2*2 = 4
  particles->addParticle({50, 0, 0}, {0, 0, 4}, 1, sigma, epsilon);   // m*v^2 = 1*4*4 = 16
  particles->addParticle({0, 100, 0}, {3, 4, 0}, 2, sigma, epsilon);  // m*v^2 = 2*5*5 = 50
  particles->addParticle({0, 0, 150}, {1, 4, 8}, 3, sigma, epsilon);  // m*v^2 = 3*9*9 = 243

  // Initialize thermostat for temperature calculation
  Thermostat thermostat(*particles, 1, 1, 1);  // don't care values except 'particles'

  ASSERT_EQ(21, thermostat.calculateCurrentTemperature());  // compare with hand-calculated temperature value
}

/**
 * @brief Check that the particle generator initializes the temperature correctly.
 */
TEST_F(ThermostatTest, CheckGeneratorInitialTemperature) {
  // Initialize cuboid
  // Initialize cuboid
  generator.queueCuboid({0.0, 0.0, 0.0},  // don't care
                        {0.0, 0.0, 0.0},  // cuboid stationary
                        {100, 50, 1},     // 5000 particles
                        1.0,              // don't care
                        1.0,              // mass is 1 for easier calculations
                        THERMO_CUBOID_TEMP,
                        1.0,  // don't care
                        5.0   // don't care
  );
  generator.generate(*particles);

  // Initialize thermostat for temperature calculation
  Thermostat thermostat(*particles, 1, 1, 1);  // don't care values except 'particles'

  ASSERT_NEAR(thermostat.calculateCurrentTemperature(), THERMO_CUBOID_TEMP, 1);
}

// TODO: Uncomment before committing
TEST_F(ThermostatTest, HoldingTemperature) {
  // Set up simulation and thermostat
  Thermostat thermostat(*particles, THERMO_FREQUENCY, THERMO_HOLDING_TEMP);

  // Add objects (particles) to simulation
  // Add objects (particles) to simulation
  generator.queueCuboid(  // cuboid from collision8000.yaml
      {10.0, 10.0, 0.0}, {0.0, 0.0, 0.0}, {120, 60, 1}, 1.1225, 1.0, THERMO_HOLDING_TEMP, 1.0, 5.0);
  generator.generate(*particles);

  // Run simulation and check final temperature
  SimulationConfig config;
  config.tEnd = SIM_END_TIME;
  config.deltaT = SIM_DT;
  config.simulationMode = SimulationMode::BENCHMARK;
  config.domainSize = std::array<double, 3>{200.0, 100.0, 10.0};
  config.linkedCellCutoff = SIM_FORCE_CUTOFF_RADIUS;
  config.boundaryTypes =
      std::array<BoundaryType, 6>{BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
                                  BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW};

  ThermostatTestingSimulation simulation(std::move(particles), thermostat, [](bool) {}, config);
  simulation.run();

  ASSERT_NEAR(thermostat.calculateCurrentTemperature(), THERMO_HOLDING_TEMP, THERMO_TEMP_TOLERANCE);
}

TEST_F(ThermostatTest, CoolingDown) {
  // Set up simulation and thermostat
  Thermostat thermostat(*particles, THERMO_FREQUENCY, THERMO_COOLING_TARGET_TEMP, THERMO_COOLING_TEMP_DELTA);

  // Add objects (particles) to simulation
  // Add objects (particles) to simulation
  generator.queueCuboid(  // cuboid from collision8000.yaml
      {10.0, 10.0, 0.0}, {0.0, 0.0, 0.0}, {120, 60, 1}, 1.1225, 1.0, THERMO_COOLING_START_TEMP, 1.0, 5.0);
  generator.generate(*particles);

  // Define temperature checking function - check if every temperature jump is in the 'tempDelta' range
  double lastTemp = 0;
  auto heatingChecker = [&](bool beforeUpdate) {
    double currentTemp = thermostat.calculateCurrentTemperature();
    if (beforeUpdate) {
      // Called before thermostat application, capture last temperature of the system
      lastTemp = currentTemp;
    } else {
      // Called after thermostat application, run checks
      EXPECT_LE(currentTemp, lastTemp) << "Temperature should drop during cooling!";  // did we cool down?
      EXPECT_LE(lastTemp - currentTemp,
                THERMO_COOLING_TEMP_DELTA +
                    0.01);  // is the jump less than or equal to the max delta (accounting for rounding errors)?
    }
  };

  // Run simulation with checks every temperature update using the lambda function
  SimulationConfig config;
  config.tEnd = SIM_END_TIME;
  config.deltaT = SIM_DT;
  config.simulationMode = SimulationMode::BENCHMARK;
  config.domainSize = std::array<double, 3>{200.0, 100.0, 10.0};
  config.linkedCellCutoff = SIM_FORCE_CUTOFF_RADIUS;
  config.boundaryTypes =
      std::array<BoundaryType, 6>{BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
                                  BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW};
  // Forces are now configured via forceConfigs, but for this test we add forces directly in the constructor

  ThermostatTestingSimulation simulation(std::move(particles), thermostat, heatingChecker, config);
  simulation.run();
}

TEST_F(ThermostatTest, HeatingUp) {
  // Set up simulation and thermostat
  Thermostat thermostat(*particles, THERMO_FREQUENCY, THERMO_HEATING_TARGET_TEMP, THERMO_HEATING_TEMP_DELTA);

  // Add objects (particles) to simulation
  // Add objects (particles) to simulation
  generator.queueCuboid(  // cuboid from collision8000.yaml
      {10.0, 10.0, 0.0}, {0.0, 0.0, 0.0}, {120, 60, 1}, 1.1225, 1.0, THERMO_HEATING_START_TEMP, 1.0, 5.0);
  generator.generate(*particles);

  // Define temperature checking function - check if every temperature jump is in the 'tempDelta' range
  double lastTemp = 0;
  auto heatingChecker = [&](bool beforeUpdate) {
    double currentTemp = thermostat.calculateCurrentTemperature();
    if (beforeUpdate) {
      // Called before thermostat application, capture last temperature of the system
      lastTemp = currentTemp;
    } else {
      // Called after thermostat application, run checks
      EXPECT_LE(lastTemp, currentTemp) << "Temperature should rise during heating!";  // did we heat up?
      EXPECT_LE(currentTemp - lastTemp,
                THERMO_HEATING_TEMP_DELTA +
                    0.01);  // is the jump less than or equal to the max delta (accounting for rounding errors)?
    }
  };

  // Run simulation with checks every temperature update using the lambda function
  SimulationConfig config;
  config.tEnd = SIM_END_TIME;
  config.deltaT = SIM_DT;
  config.simulationMode = SimulationMode::BENCHMARK;
  config.domainSize = std::array<double, 3>{200.0, 100.0, 10.0};
  config.linkedCellCutoff = SIM_FORCE_CUTOFF_RADIUS;
  config.boundaryTypes =
      std::array<BoundaryType, 6>{BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
                                  BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW};

  ThermostatTestingSimulation simulation(std::move(particles), thermostat, heatingChecker, config);
  simulation.run();
}
