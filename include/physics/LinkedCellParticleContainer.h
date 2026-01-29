/**
* @file LinkedCellParticleContainer.h
 *
 *
 */
#pragma once
#include <array>
#include <functional>
#include <vector>

#include "physics/ParticleContainer.h"

/**
 * @enum BoundaryType
 * @brief Defines behavior of particles at  boundaries.
 */
enum class BoundaryType {
  OUTFLOW,     ///< Particles leave the domain and are deleted
  REFLECTIVE,  ///< Particles bounce back (velocity reflection)
  PERIODIC     ///< Particles wrap around to opposite boundary
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
   * @brief Constructor for LinkedCellParticleContainer
   * @param domainDims Size of the simulation domain (x, y, z)
   * @param cutoffRadius The cutoff radius for interactions
   * @param boundaryTypes Boundary types for each face (left, right, bottom, top, back, front)
   */
  LinkedCellParticleContainer(const std::array<double, 3>& domain_dims, double cutoff_radius,
                              const std::array<BoundaryType, 6>& boundary_types);

  ~LinkedCellParticleContainer() = default;

  using ParticleContainer::addParticle;
  // Override addParticle to place particle in correct cell
  Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m) override;
  Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, double sigma,
                        double epsilon) override;
  Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, int type, double sigma,
                        double epsilon) override;
  Particle* addParticle(const Particle* p) override;
  Particle* addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, std::array<double, 3> f,
                        std::array<double, 3> oldF, int type, double sigma, double epsilon) override;

  /**
   * @brief Update cell assignments after particle positions change
   */
  void updateCells();

  /**
   * @brief Handle outflow boundaries by removing particles outside domain and updating cells
   */
  void handleOutflowBoundaries();

  /**
   * @brief Iterate over all distinct particle pairs within cutoff distance (std::function version)
   * @param pairFunc Function to apply to each pair (p1, p2)
   * @note This version has overhead from std::function type erasure. Use the template version for performeance-critical code.
   */
  void iteratePairs(const std::function<void(Particle&, Particle&)>& pairFunc) const;

  /**
   * @brief Iterate over all distinct particle pairs within cutoff distance (template version)
   * @param pairFunc Function to apply to each pair (p1, p2)
   * @note This template version is designed to eliminate std::function overhead
   */
  template <typename F>
  void iteratePairsTemplate(F&& pairFunc) const {
    const int nx = numCells[0];
    const int ny = numCells[1];
    const int nz = numCells[2];

    const int layerSize = nx * ny;

    // offsets for 13 forward neighbors
    static constexpr std::array<std::array<int, 3>, 13> neighborOffsets = {{{1, 0, 0},
                                                                            {1, 1, 0},
                                                                            {0, 1, 0},
                                                                            {-1, 1, 0},
                                                                            {-1, -1, 1},
                                                                            {0, -1, 1},
                                                                            {1, -1, 1},
                                                                            {-1, 0, 1},
                                                                            {0, 0, 1},
                                                                            {1, 0, 1},
                                                                            {-1, 1, 1},
                                                                            {0, 1, 1},
                                                                            {1, 1, 1}}};

    // Compute interactions between particles in the same cell, ignoring halo cells
    for (size_t cdx = 0; cdx < cells.size(); ++cdx) {
      if (getCellType(cdx) == CellType::HALO)
        continue;

      const auto& cell = cells[cdx];
      const size_t cellSize = cell.size();

      // Intra-cell pairs
      for (size_t i = 0; i < cellSize; ++i) {
        Particle* __restrict__ pi = cell[i];
        for (size_t j = i + 1; j < cellSize; ++j) {
          Particle* __restrict__ pj = cell[j];
          pairFunc(*pi, *pj);
        }
      }

      // Convert cell index to 3D coordinates
      const int idx = static_cast<int>(cdx);
      const int iz = idx / layerSize;
      const int remainder = idx % layerSize;
      const int iy = remainder / nx;
      const int ix = remainder % nx;

      // Compute interactions with forward neighbor cells (avoid double counting)
      for (const auto& off : neighborOffsets) {
        const int nix = ix + off[0];
        const int niy = iy + off[1];
        const int niz = iz + off[2];

        // If neighbor is out of bounds, skip it
        if (nix < 0 || nix >= nx || niy < 0 || niy >= ny || niz < 0 || niz >= nz) {
          continue;
        }

        // Convert neighbor coordinates back to 1D index
        const int nIndex = (niz * layerSize) + (niy * nx) + nix;
        const auto& ncell = cells[static_cast<size_t>(nIndex)];
        const size_t ncellSize = ncell.size();

        // Inter-cell pairs: loops over all pairs
        for (size_t i = 0; i < cellSize; ++i) {
          Particle* __restrict__ pi = cell[i];
          for (size_t j = 0; j < ncellSize; ++j) {
            Particle* __restrict__ pj = ncell[j];
            pairFunc(*pi, *pj);
          }
        }
      }
    }
  }

  /**
   * @brief Iterate over all particles in a cell and its neighbors
   * @param cdx Index of the central cell
   * @param func Function to apply
   */
  void iterateCellNeighbors(size_t cdx, const std::function<void(Particle&, Particle&)>& func) const;

  /**
   * @brief Access cell by grid coordinates.
   * @param cx x-index
   * @param cy y-index
   * @param cz z-index
   * @return Reference to vector of particles in the cell
   */
  std::vector<Particle*>& cell_at(int cx, int cy, int cz);

  // Getters
  [[nodiscard]] std::array<double, 3> domain_dims() const { return domainDims; }
  [[nodiscard]] std::array<double, 3> domain_origin() const { return domainOrigin; }
  [[nodiscard]] double cutoff_radius() const { return cutoffRadius; }
  [[nodiscard]] std::array<double, 3> cell_size() const { return cellSize; }
  [[nodiscard]] std::array<int, 3> num_cells() const { return numCells; }
  [[nodiscard]] std::array<BoundaryType, 6> boundary_types() const { return boundaryTypes; }

 private:
  std::array<double, 3> domainDims;                ///< Dimensions of the simulation domain (x, y, z)
  std::array<double, 3> domainOrigin{0., 0., 0.};  ///< Coordinates of the domain origin (bottom-left-front corner)
  double cutoffRadius;                             ///< Cutoff radius for particle-particle interactions
  std::array<double, 3> cellSize;                  ///< Dimensions of a single cell
  CellType cellType;                               ///< Type of the current cell
  std::array<int, 3> numCells;                     ///< Number of cells in each dimension (including halo layer)
  std::array<BoundaryType, 6> boundaryTypes = {
      BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW, BoundaryType::OUTFLOW,
      BoundaryType::OUTFLOW, BoundaryType::OUTFLOW};  ///< Boundary conditions for the 6 faces (x-, x+, y-, y+, z-, z+)
  std::vector<std::vector<Particle*>>
      cells;  ///< Linearized vector of cells, where each cell contains pointers to particles

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
  /**
 * @brief Handle outflow boundary: delete particles outside domain
 */
  void handleOutflow();

  /**
   * @brief Handle periodic boundaries: wrap particles around domain
   */
  void handlePeriodicBoundaries();

  /**
 * @brief Check if position is inside domain
 */
  [[nodiscard]] bool isInsideDomain(const std::array<double, 3>& pos) const;

  /**
* @brief Map the position of a particle's position to the index of the respective cell
*/
  [[nodiscard]] int getCellIndex(const std::array<double, 3>& x) const;
};