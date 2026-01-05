#include <memory>

#include <gtest/gtest.h>

#include "ParticleGenerator.h"
#include "Simulation.h"

// Simulation defines
constexpr double SIM_END_TIME = 0.2;
constexpr double SIM_DT = 0.002;      // ~100 simulation iterations
constexpr int THERMO_FREQUENCY = 20;  // ~5 thermostat applications

constexpr double SIM_FORCE_EPSILON = 5.0;
constexpr double SIM_FORCE_SIGMA = 1.0;

constexpr double SIM_FORCE_CUTOFF_RADIUS = 3.0;

// Temperature test defines
constexpr double THERMO_TEMP_TOLERANCE = 1e-4;  // how far off is the thermostat allowed to be?

// Temperature holding test defines
constexpr double THERMO_HOLDING_TEMP = 50.0;

class ThermostatTestingSimulation : public BaseSimulation {
 private:
  Thermostat& thermostat;

 protected:
  void setupSimulation() override { /* empty override for compilation */ }
  void runFileOutput(int frequency, const std::string& outputBaseName) override { /* empty as it will not be called */ }
  void runBenchmark() override {
    // No benchmarking, just simulate and regularly apply thermostat
    double current_time = 0;
    int iteration = 0;
    const int thermostatFrequency = thermostat.getUpdateFrequency();
    while (current_time < end_time) {
      forceCalc->calculateX(dt);
      for (auto& p : *particles) {
        p.setOldF(p.getF());
      }
      forceCalc->calculateF();
      forceCalc->calculateV(dt);

      iteration++;
      if (iteration % thermostatFrequency == 0) {
        thermostat.updateTemperature();
      }
      current_time += dt;
    }
  }

 public:
  explicit ThermostatTestingSimulation(std::unique_ptr<ParticleContainer> particles, Thermostat& thermostat)
      : BaseSimulation(SIM_END_TIME, SIM_DT,
                       0,                         // don't care about file output
                       "",                        // don't care about file output
                       SimulationMode::BENCHMARK  // we want runBenchmark() to run
                       ),
        thermostat(thermostat) {

    this->particles = std::move(particles);
    forceCalc = std::make_unique<LennardJonesForce>(*this->particles, SIM_FORCE_EPSILON, SIM_FORCE_SIGMA,
                                                    SIM_FORCE_CUTOFF_RADIUS);
  }
  ~ThermostatTestingSimulation() override = default;
};

class ThermostatTest : public testing::Test {
 protected:
  std::unique_ptr<ParticleContainer> particles;  // ownership transferred to simulation later on
  ParticleGenerator generator;

  ThermostatTest() : particles(std::make_unique<ParticleContainer>()), generator(*particles) {}
  ~ThermostatTest() override = default;
};

TEST_F(ThermostatTest, CheckGeneratorInitialTemperature) {
  // TODO: Check that the particle generator initializes the temperature correctly
}

TEST_F(ThermostatTest, CheckThermostatInitialTemperature) {
  // TODO: Check that the thermostat initializes the temperature correctly
}

TEST_F(ThermostatTest, CheckThermostatTemperatureCalculation) {
  // TODO: Check the correctness of the current temperature calculation from kinetic energy
}

TEST_F(ThermostatTest, Holding) {
  // TODO: This test somehow succeeds although the internal temperature values are wrong for the first ~3 thermostat applications
  // Set up simulation and thermostat
  Thermostat thermostat(*particles, THERMO_FREQUENCY, THERMO_HOLDING_TEMP);

  // Add objects (particles) to simulation
  generator.generateCuboid(  // cuboid from collision8000.yaml
      {10.0, 10.0, 0.0}, {0.0, 0.0, 0.0}, {120, 60, 1}, 1.1225, 1.0, THERMO_HOLDING_TEMP, 1.0, 5.0);

  // Check initial temperature being correct
  ThermostatTestingSimulation simulation(std::move(particles), thermostat);
  simulation.run();

  double finalTemp = thermostat.calculateCurrentTemperature();

  std::cout << "\n[   INFO   ] Final Temperature: " << finalTemp << " | Target Holding Temp: " << THERMO_HOLDING_TEMP
            << "\n"
            << std::endl;
  // -------------------------

  ASSERT_NEAR(thermostat.calculateCurrentTemperature(), THERMO_HOLDING_TEMP, THERMO_TEMP_TOLERANCE);
}

TEST_F(ThermostatTest, Cooling) {
  // TODO
}

TEST_F(ThermostatTest, Heating) {
  // TODO
}
