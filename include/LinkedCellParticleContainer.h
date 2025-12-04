/**
* @file LinkedCellParticleContainer.h
 *
 *
 */
#pragma once
#include <array>
#include <functional>
#include <vector>

#include "ParticleContainer.h"

/**
 * @class LinkedCellParticleContainer
 * @brief Iterable container class which extends ParticleContainer used for storing particles for a simulation.
 * It offers methods for adding, accessing particles and iterating through a set of particles
 * in an easy and efficient manner.
 * Additionally, it implements boundary conditions and the Linked Cell Algorithm.
 *
 */
class LinkedCellParticleContainer : public ParticleContainer {
 public:
  /**
   * @enum BoundaryType
   * @brief Defines behavior of particles at  boundaries.
   */
  enum class BoundaryType {
    OUTFLOW,    ///< Particles leave the domain and are deleted
    REFLECTIVE, ///< Particles bounce back (velocity reflection)
    PERIODIC    ///< Particles wrap around to opposite boundary
  };

  /**
  * @enum CellType
  * @brief Defines type of particles in the cells.
  */
  enum class CellType {
    INNER,
    BOUNDARY,
    HALO,
  };


  /**
   * @brief Constructor for LinkedCellParticleContainer
   * @param domainDims Size of the simulation domain (x, y, z)
   * @param cutoffRadius The cutoff radius for interactions
   * @param boundaryTypes Boundary types for each face (left, right, bottom, top, back, front)
   */
  LinkedCellParticleContainer(const std::array<double, 3>& domain_dims, double cutoff_radius,
                              const std::array<BoundaryType, 6>& boundary_types)
      : domainDims(domain_dims), cutoffRadius(cutoff_radius), boundaryTypes(boundary_types) {}

  ~LinkedCellParticleContainer() = default;

  // override addParticle to place particle in correct cell
  void addParticle(std::array<double, 3> x, std::array<double, 3> v, double m);
  void addParticle(const Particle* p);

  /**
   * @brief Update cell assignments after particle positions change
   */
  void updateCells();

  /**
   * @brief Apply boundary conditions to all particles
   */
  void applyBoundaryConditions();

  /**
   * @brief Iterate over all distinct particle pairs within cutoff distance
   * @param pairFunc Function to apply to each pair (p1, p2)
   */
  void iteratePairs(const std::function<void(Particle&, Particle&)>& pairFunc);

  /**
   * @brief Iterate over all particles in a cell and its neighbors
   * @param cellIndex Index of the central cell
   * @param func Function to apply
   */
  void iterateCellNeighbors(size_t cellIndex,
                            const std::function<void(Particle&, Particle&)>& func);

  // getters
[[nodiscard]] std::array<double, 3> domain_dims() const { return domainDims; }
  [[nodiscard]] std::array<double, 3> domain_origin() const { return domainOrigin; }
  [[nodiscard]] double cutoff_radius() const { return cutoffRadius; }
  [[nodiscard]] std::array<double, 3> cell_size() const { return cellSize; }
  [[nodiscard]] std::array<int, 3> num_cells() const { return numCells; }
  [[nodiscard]] std::array<BoundaryType, 6> boundary_types() const { return boundaryTypes; }

 private:
    std::array<double, 3> domainDims;
    std::array<double, 3> domainOrigin{0.,0.,0.};
    double cutoffRadius;
    std::array<double, 3> cellSize; // size of ech cell
    std::array<int, 3> numCells; // number of cells in each dimension (including halo)
    std::array<BoundaryType, 6> boundaryTypes = {BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
      BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW}; // boundary types for 6 faces (left, right, bottom, top, back, front)
    std::vector<std::vector<Particle*>> cells; // 1D array of cells;
    std::vector<Particle> particles; // stores all particles

  /**
* @brief Initialize cell grid
*/
    void initCells();

  /**
 * @brief Get all neighbor cell indices for a cell
 */
  [[nodiscard]] std::vector<size_t> getNeighborCellIndices(int cellIndex) const;

  /**
 * @brief Handle outflow boundary: delete particles outside domain
 */
  void handleOutflow();

  /**
   * @brief Handle reflective boundary using ghost particles
   */
  void handleReflective();

  /**
 * @brief Check if position is inside domain
 */
  [[nodiscard]] bool isInsideDomain(const std::array<double, 3>& pos) const;

};