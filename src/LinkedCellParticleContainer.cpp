#include "../include/LinkedCellParticleContainer.h"

#include <cmath>

#include "spdlog/fmt/bundled/format.h"

LinkedCellParticleContainer::LinkedCellParticleContainer(const std::array<double, 3>& domain_dims, double cutoff_radius,
                                                         const std::array<BoundaryType, 6>& boundary_types)
    : domainDims(domain_dims), cutoffRadius(cutoff_radius), boundaryTypes(boundary_types) {
  if (cutoffRadius <= 0.0) {
    throw std::invalid_argument("Cutoff radius must be positive");
  }
  initCells();
}

void LinkedCellParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m) {
  ParticleContainer::addParticle(x, v, m);
  Particle& p = (*this)[this->size() - 1];  // the newly added particle
  const int cdx = getCellIndex(p.getX());
  cells[cdx].push_back(&p);
}

void LinkedCellParticleContainer::addParticle(std::array<double, 3> x, std::array<double, 3> v, double m, double sigma,
                                              double epsilon) {
  ParticleContainer::addParticle(x, v, m, sigma, epsilon);
  Particle& p = (*this)[this->size() - 1];  // the newly added particle
  const int cdx = getCellIndex(p.getX());
  cells[cdx].push_back(&p);
}

void LinkedCellParticleContainer::addParticle(const Particle* p) {
  addParticle(p->getX(), p->getV(), p->getM(), p->getSigma(), p->getEpsilon());
}

void LinkedCellParticleContainer::updateCells() {
  for (auto& cell : cells) {
    cell.clear();
  }

  const std::size_t n = this->size();
  for (std::size_t i = 0; i < n; ++i) {
    Particle& p = (*this)[i];
    const std::size_t c = getCellIndex(p.getX());
    cells[c].push_back(&p);
  }
}

LinkedCellParticleContainer::CellType LinkedCellParticleContainer::getCellType(size_t cdx) const {
  const int nx = numCells[0];
  const int ny = numCells[1];
  const int nz = numCells[2];

  const int idx = static_cast<int>(cdx);
  const int layerSize = nx * ny;

  const int iz = idx / layerSize;
  const int remainder = idx % layerSize;
  const int iy = remainder / nx;
  const int ix = remainder % nx;

  // outermost layer in any dimension -> Halo
  if (ix == 0 || ix == nx - 1 || iy == 0 || iy == ny - 1 || iz == 0 || iz == nz - 1) {
    return CellType::HALO;
  }

  // next layer  from halo -> Boundary
  if (ix == 1 || ix == nx - 2 || iy == 1 || iy == ny - 2 || iz == 1 || iz == nz - 2) {
    return CellType::BOUNDARY;
  }

  return CellType::INNER;
}

void LinkedCellParticleContainer::handleOutflowBoundaries() {
  handleOutflow();
  updateCells();
}
void LinkedCellParticleContainer::iteratePairs(const std::function<void(Particle&, Particle&)>& pairFunc) const {
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

    auto& cell = cells[cdx];
    for (std::size_t i = 0; i < cell.size(); ++i) {
      for (std::size_t j = i + 1; j < cell.size(); ++j) {
        pairFunc(*cell[i], *cell[j]);
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
      auto& ncell = cells[static_cast<size_t>(nIndex)];

      // Compute interactions between current cell and neighbor
      for (auto* pi : cell) {
        for (auto* pj : ncell) {
          pairFunc(*pi, *pj);
        }
      }
    }
  }
}
void LinkedCellParticleContainer::iterateCellNeighbors(size_t cdx,
                                                       const std::function<void(Particle&, Particle&)>& func) const {
  auto& centerCell = cells[cdx];  // the currently selected cell at the center ouf the surrounding 26 neighbors
  const auto neighbors = getNeighborCellIndices(static_cast<int>(cdx));
  for (auto nidx : neighbors) {
    auto& ncell = cells[nidx];
    for (auto* pi : centerCell) {
      for (auto* pj : ncell) {
        func(*pi, *pj);
      }
    }
  }
}
void LinkedCellParticleContainer::initCells() {
  // Calculate cell size and number of cells in each dimension
  for (int i = 0; i < 3; ++i) {
    int inner = static_cast<int>(std::floor(domain_dims()[i] / cutoff_radius()));
    if (inner < 1)
      inner = 1;
    // cell size is chosen so that all inner cell exactly cover the domain in dim i
    cellSize[i] = domainDims[i] / static_cast<double>(inner);
    numCells[i] = inner + 2;  // +2 accounts for halo layer on each side
  }
  const int nx = numCells[0];
  const int ny = numCells[1];
  const int nz = numCells[2];

  // Allocate cell grid (linearized 3D array)
  cells.assign(nx * ny * nz, {});
}
std::vector<size_t> LinkedCellParticleContainer::getNeighborCellIndices(int cdx) const {
  std::vector<std::size_t> neighbors;
  const int nx = numCells[0];
  const int ny = numCells[1];
  const int nz = numCells[2];
  const int layerSize = nx * ny;

  // Convert 1D index to 3D coordinates
  const int idx = cdx;
  const int iz = idx / layerSize;
  const int remainder = idx % layerSize;
  const int iy = remainder / nx;
  const int ix = remainder % nx;

  // Iterate over all 27 neighboring cells (including center cell)
  for (int dz = -1; dz <= 1; ++dz) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        const int nix = ix + dx;
        const int niy = iy + dy;
        const int niz = iz + dz;

        // Skip f neighbor is out of bounds
        if (nix < 0 || nix >= nx || niy < 0 || niy >= ny || niz < 0 || niz >= nz) {
          continue;
        }
        // Convert 3D coordinates back to 1D index
        const int nIndex = (niz * layerSize) + (niy * nx) + nix;
        neighbors.push_back(nIndex);
      }
    }
  }
  return neighbors;
}
void LinkedCellParticleContainer::handleOutflow() {
  std::vector<size_t> indicesToRemove;
  // Check each particle if it's outside the domain
  for (size_t i = 0; i < this->size(); ++i) {
    Particle& p = (*this)[i];
    auto x = p.getX();
    bool outside = false;

    for (int d = 0; d < 3; ++d) {
      // lower boundary
      if (boundaryTypes[2 * d] == BoundaryType::OUTFLOW && x[d] < domainOrigin[d]) {
        outside = true;
        break;
      }
      // upper boundary
      if (boundaryTypes[2 * d + 1] == BoundaryType::OUTFLOW && x[d] > domainOrigin[d] + domainDims[d]) {
        outside = true;
        break;
      }
    }
    if (outside) {
      indicesToRemove.push_back(i);
    }
  }
  // Removing the particles in reverse order to maintain valid indexation
  for (auto idx = indicesToRemove.rbegin(); idx != indicesToRemove.rend(); ++idx) {
    this->removeParticle(*idx);
  }
}

bool LinkedCellParticleContainer::isInsideDomain(const std::array<double, 3>& pos) const {
  for (int i = 0; i < 3; ++i) {
    if (pos[i] < domainOrigin[i] || pos[i] > domainOrigin[i] + domainDims[i]) {
      return false;
    }
  }
  return true;
}
int LinkedCellParticleContainer::getCellIndex(const std::array<double, 3>& x) const {
  std::array<int, 3> idx{};
  // Calculate cell coordinates from particle position
  for (int i = 0; i < 3; ++i) {
    const double rel_dim = x[i] - domainOrigin[i];
    int cdx = static_cast<int>(std::floor(rel_dim / cell_size()[i])) + 1;  // +1 as to skip the halo cells
    // Adjust to valid cell range
    if (cdx < 0)
      cdx = 0;
    if (cdx >= numCells[i])
      cdx = numCells[i] - 1;  // Limiting inner domain to indices 1...numCells[d]-2
    idx[i] = cdx;
  }

  // Flatten 3D coordinates into 1D
  const int nx = numCells[0];
  const int ny = numCells[1];
  // reference: https://stackoverflow.com/questions/7367770/
  return (idx[2] * nx * ny) + (idx[1] * nx) + idx[0];
}
