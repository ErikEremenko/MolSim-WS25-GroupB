#!/bin/bash

# Stop the script if any command fails
set -e

# 1. Create and enter build directory
if [ ! -d "build" ]; then
    echo "📂 Creating build directory..."
    mkdir build
fi
cd build

# 2. Configure CMake (only if needed, or forced)
echo "⚙️  Configuring CMake..."
cmake .. -DCMAKE_BUILD_TYPE=Release

# 3. Compile ONLY the test executable
# -j$(nproc) uses all available CPU cores for faster compilation
echo "🔨 Building MolSimTests..."
cmake --build . --target MolSimTests -- -j$(nproc)

# 4. Run the Thermostat tests
echo "🧪 Running Thermostat Tests..."
echo "---------------------------------------------------"

# This filter runs all tests in the ThermostatTest suite
./MolSimTests --gtest_filter="ThermostatTest.*"

echo "---------------------------------------------------"
echo "✅ Done."