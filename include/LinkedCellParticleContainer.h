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
  LinkedCellParticleContainer( const std::array<double, 3>& domain_dims,
    double cutoff_radius,
    const std::array<BoundaryType, 6>& boundary_types);

  ~LinkedCellParticleContainer() = default;

  using ParticleContainer::addParticle;
  // override addParticle to place particle in correct cell
  void addParticle(std::array<double, 3> x, std::array<double, 3> v, double m);
  void addParticle(const Particle* p);

  /**
   * @brief Update cell assignments after particle positions change
   */
  void updateCells();

  /**
   * @brief Handle outflow boundaries by removing particles outside domain and updating cells
   */
  void handleOutflowBoundaries();

  /**
   * @brief Iterate over all distinct particle pairs within cutoff distance
   * @param pairFunc Function to apply to each pair (p1, p2)
   */
  void iteratePairs(const std::function<void(Particle&, Particle&)>& pairFunc) const;

  /**
   * @brief Iterate over all particles in a cell and its neighbors
   * @param cdx Index of the central cell
   * @param func Function to apply
   */
  void iterateCellNeighbors(size_t cdx,
                            const std::function<void(Particle&, Particle&)>& func) const;

  // getters
[[nodiscard]] std::array<double, 3> domain_dims() const { return domainDims; }
  [[nodiscard]] std::array<double, 3> domain_origin() const { return domainOrigin; }
  [[nodiscard]] double cutoff_radius() const { return cutoffRadius; }
  [[nodiscard]] std::array<double, 3> cell_size() const { return cellSize; }
  [[nodiscard]] std::array<int, 3> num_cells() const { return numCells; }
  [[nodiscard]] std::array<BoundaryType, 6> boundary_types() const { return boundaryTypes; }

 private:
    std::array<double, 3> domainDims;             ///< Dimensions of the simulation domain (x, y, z)
    std::array<double, 3> domainOrigin{0., 0., 0.}; ///< Coordinates of the domain origin (bottom-left-front corner)
    double cutoffRadius;                          ///< Cutoff radius for particle-particle interactions
    std::array<double, 3> cellSize;               ///< Dimensions of a single cell
    CellType cellType;                            ///< Type of the current cell
    std::array<int, 3> numCells;                  ///< Number of cells in each dimension (including halo layer)
    std::array<BoundaryType, 6> boundaryTypes = {
        BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
        BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW}; ///< Boundary conditions for the 6 faces (x-, x+, y-, y+, z-, z+)
    std::vector<std::vector<Particle*>> cells;    ///< Linearized vector of cells, where each cell contains pointers to particles

  /**
* @brief Initialize cell grid
*/
    void initCells();

  /**
 * @brief Get all 3*3*3=27 neighbor cell indices, including the centered cell itself
 */
  [[nodiscard]] std::vector<size_t> getNeighborCellIndices(int cdx) const;

  /**
* @brief Returns the CellType of the given cell index
*/
  [[nodiscard]] CellType getCellType(size_t cdx) const;

  /**
 * @brief Handle outflow boundary: delete particles outside domain
 */
  void handleOutflow();

  /**
 * @brief Check if position is inside domain
 */
  [[nodiscard]] bool isInsideDomain(const std::array<double, 3>& pos) const;

  /**
* @brief Map the position of a particle's position to the index of the respective cell
*/
  [[nodiscard]] int getCellIndex(const std::array<double, 3>& x) const;

};