#pragma once

#include "physics/ForceCalc.h"

#include <cmath>
#include <vector>

/**
 * @class SmoothedLJForce
 * @brief Smoothed Lennard-Jones potential with continuous force at cutoff
 *
 * Uses the standard Lennard-Jones force up to the smoothing radius, then
 * tapers to zero at the cutoff radius. This avoids discontinuities
 * in the force. Mixed particle types use Lorentz-Berthelot rules.
 */
class SmoothedLJForce final : public ForceCalc {
 private:
  /// Global Lennard-Jones epsilon parameter
  const double epsilon;
  /// Global Lennard-Jones sigma parameter
  const double sigma;
  /// Cutoff radius (r_c) beyond which interactions are zero
  const double cutoffRadius;
  /// Smoothing start radius (r_l) below which standard LJ is used
  const double smoothingRadius;

  // Precomputed constants
  double cutoffRadiusSq;     ///< r_c^2 for squared distance comparisons
  double smoothingRadiusSq;  ///< r_l^2 for squared distance comparisons
  double rcMinusRlCubed;     ///< (r_c - r_l)^3 for smoothing function

  /// Lookup tables for pair interactions (indexed by type_i * tableWidth + type_j)
  std::vector<double> pairEpsilon;  ///< Mixed epsilon values
  std::vector<double> pairSigma6;   ///< sigma_ij^6
  std::vector<double> pairSigma12;  ///< sigma_ij^12
  int tableWidth = 0;

 public:
  /**
   * @param particles ParticleContainer that stores the particles
   * @param epsilon Epsilon in the Lennard-Jones potential formula
   * @param sigma Sigma in the Lennard-Jones potential formula
   * @param cutoffRadius Cutoff radius r_c beyond which interactions are ignored
   * @param smoothingRadius Smoothing radius r_l below which standard LJ is used
   */
  SmoothedLJForce(ParticleContainer& particles, double epsilon, double sigma, double cutoffRadius,
                  double smoothingRadius);

  void calculateF() override;
  void precomputeConstants() override;

 private:
  void calculateFDirectSum() const;
  void calculateFLinkedCell();
  void calcFPeriodicPair(Particle* p1, Particle* p2) const;
  void applyPeriodicBoundaries(LinkedCellParticleContainer* lc) const;

  /**
   * @brief Compute and apply smoothed LJ pair force (must be defined in header for proper inlining)
   */
  __attribute__((always_inline)) inline void applySmoothedPairForce(Particle& p_i, Particle& p_j) const {
    const auto& xi = p_i.getX();
    const auto& xj = p_j.getX();

    const double dx = xj[0] - xi[0];
    const double dy = xj[1] - xi[1];
    const double dz = xj[2] - xi[2];
    const double distSq = dx * dx + dy * dy + dz * dz;

    if (distSq >= cutoffRadiusSq) {
      return;
    }

    // Use precomputed lookup tables
    const int idx = p_i.getType() * tableWidth + p_j.getType();
    const double epsilon_ij = pairEpsilon[idx];
    const double sigma6 = pairSigma6[idx];
    const double sigma12 = pairSigma12[idx];

    const double inv_distSq = 1.0 / distSq;
    const double inv_distSq3 = inv_distSq * inv_distSq * inv_distSq;
    const double sigma6_d6 = sigma6 * inv_distSq3;
    const double sigma12_d12 = sigma12 * inv_distSq3 * inv_distSq3;

    const double U_LJ = 4.0 * epsilon_ij * (sigma12_d12 - sigma6_d6);
    const double F_LJ_scalar = 24.0 * epsilon_ij * inv_distSq * (sigma6_d6 - 2.0 * sigma12_d12);

    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;

    if (distSq <= smoothingRadiusSq) {
      fx = F_LJ_scalar * dx;
      fy = F_LJ_scalar * dy;
      fz = F_LJ_scalar * dz;
    } else {
      const double d = std::sqrt(distSq);
      const double dMinusRl = d - smoothingRadius;
      const double dMinusRl2 = dMinusRl * dMinusRl;
      const double S = 1.0 - dMinusRl2 * (3.0 * cutoffRadius - smoothingRadius - 2.0 * d) / rcMinusRlCubed;
      const double dS_dd = -6.0 * dMinusRl * (cutoffRadius - d) / rcMinusRlCubed;

      const double F_total_scalar = S * F_LJ_scalar + U_LJ * dS_dd / d;
      fx = F_total_scalar * dx;
      fy = F_total_scalar * dy;
      fz = F_total_scalar * dz;
    }

    auto& fi = p_i.getF();
    auto& fj = p_j.getF();
    fi[0] += fx;
    fi[1] += fy;
    fi[2] += fz;
    fj[0] -= fx;
    fj[1] -= fy;
    fj[2] -= fz;
  }
};
