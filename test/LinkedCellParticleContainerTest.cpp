#include <gtest/gtest.h>

#include "LinkedCellParticleContainer.h"

class LinkedCellParticleContainerTest : public ::testing::Test {
 protected:
  std::array<double, 3> domainDims = {10.0, 10.0, 10.0};
  double cutoffRadius = 2.5;
  std::array<LinkedCellParticleContainer::BoundaryType, 6> outflowBoundaries = {
      LinkedCellParticleContainer::BoundaryType::OUTFLOW, LinkedCellParticleContainer::BoundaryType::OUTFLOW,
      LinkedCellParticleContainer::BoundaryType::OUTFLOW, LinkedCellParticleContainer::BoundaryType::OUTFLOW,
      LinkedCellParticleContainer::BoundaryType::OUTFLOW, LinkedCellParticleContainer::BoundaryType::OUTFLOW};
};

// Test that the container initializes with correct domain dimensions
TEST_F(LinkedCellParticleContainerTest, ConstructorInitializesDomain) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  EXPECT_EQ(lpc.domain_dims()[0], 10.0);
  EXPECT_EQ(lpc.domain_dims()[1], 10.0);
  EXPECT_EQ(lpc.domain_dims()[2], 10.0);
  EXPECT_EQ(lpc.cutoff_radius(), 2.5);
}

// Test that invalid cutoff radius throws exception
TEST_F(LinkedCellParticleContainerTest, InvalidCutoffThrows) {
  EXPECT_THROW(LinkedCellParticleContainer(domainDims, 0.0, outflowBoundaries), std::invalid_argument);
  EXPECT_THROW(LinkedCellParticleContainer(domainDims, -1.0, outflowBoundaries), std::invalid_argument);
}

// Test that the cell grid is initialized correctly
TEST_F(LinkedCellParticleContainerTest, CellGridInitialization) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  // With domain 10, cutoff  2.5, we expect floor(10/2.5) = 4 inner cells per dim + 2 halo cells = 6 cells per dim
  auto numCells = lpc.num_cells();
  EXPECT_EQ(numCells[0], 6);
  EXPECT_EQ(numCells[1], 6);
  EXPECT_EQ(numCells[2], 6);
}

// Test that container is initially empty
TEST_F(LinkedCellParticleContainerTest, InitiallyEmpty) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);
  EXPECT_EQ(lpc.size(), 0);
}

// Test adding particles
TEST_F(LinkedCellParticleContainerTest, AddParticle) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  EXPECT_EQ(lpc.size(), 1);

  lpc.addParticle({2.0, 2.0, 2.0}, {1.0, 1.0, 1.0}, 2.0);
  EXPECT_EQ(lpc.size(), 2);
}

// Test that particles can be accessed after addition
TEST_F(LinkedCellParticleContainerTest, AccessParticleAfterAdd) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  lpc.addParticle({5.0, 5.0, 5.0}, {1.0, 2.0, 3.0}, 1.5);

  EXPECT_DOUBLE_EQ(lpc[0].getX()[0], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[1], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getX()[2], 5.0);
  EXPECT_DOUBLE_EQ(lpc[0].getM(), 1.5);
}

// Test pair iteration counts correct number of pairs
TEST_F(LinkedCellParticleContainerTest, IteratePairsCountsCorrectly) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  // Add 3 particles close together (within cutoff)
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  lpc.addParticle({5.5, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  lpc.addParticle({5.0, 5.5, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  int pairCount = 0;
  lpc.iteratePairs([&pairCount](Particle&, Particle&) { pairCount++; });

  // 3 particles -> 3 pairs
  EXPECT_EQ(pairCount, 3);
}

// Test outflow boundary removes particles outside domain
TEST_F(LinkedCellParticleContainerTest, OutflowRemovesParticlesOutside) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  // Add particle inside domain
  lpc.addParticle({5.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);
  // Add particle outside domain (will be placed in halo initially)
  lpc.addParticle({-1.0, 5.0, 5.0}, {0.0, 0.0, 0.0}, 1.0);

  EXPECT_EQ(lpc.size(), 2);

  lpc.handleOutflowBoundaries();

  // Particle outside should be removed
  EXPECT_EQ(lpc.size(), 1);
}

// Test that particles at boundaries are handled correctly
TEST_F(LinkedCellParticleContainerTest, ParticleAtBoundary) {
  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, outflowBoundaries);

  // Particle exactly at origin (edge of domain)
  lpc.addParticle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, 1.0);
  EXPECT_EQ(lpc.size(), 1);

  lpc.handleOutflowBoundaries();
  // Should still exist (not outside)
  EXPECT_EQ(lpc.size(), 1);
}

// Test reflective boundaries setup
TEST_F(LinkedCellParticleContainerTest, ReflectiveBoundaryTypes) {
  std::array<LinkedCellParticleContainer::BoundaryType, 6> reflectiveBoundaries = {
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE,
      LinkedCellParticleContainer::BoundaryType::REFLECTIVE, LinkedCellParticleContainer::BoundaryType::REFLECTIVE};

  LinkedCellParticleContainer lpc(domainDims, cutoffRadius, reflectiveBoundaries);

  auto bt = lpc.boundary_types();
  for (int i = 0; i < 6; ++i) {
    EXPECT_EQ(bt[i], LinkedCellParticleContainer::BoundaryType::REFLECTIVE);
  }
}
