#pragma once

/**
 * @file PhysicsConstants.h
 * @brief Precomputed mathematical constants for physics calculations
 */

namespace PhysicsConstants {

/// sqrt(2)
inline constexpr double SQRT_2 = 1.4142135623730950488016887242096980785696718753769480731766797379;

/// 1 / sqrt(2)
inline constexpr double INV_SQRT_2 = 0.7071067811865475244008443621048490392848359376884740365883398689;

// LJ-potential constants

/// 2^(1/6) equilibrium distance factor for LJ-potential
/// sigma_eq = 2^(1/6) * sigma, where F(sigma_eq) = 0
inline constexpr double TWO_POW_1_6 = 1.1224620483093729814335330496791795162324111106139867534404505950;

/// 2^(1/3) = (2^(1/6))^2 - for squared distance comparisons
/// Used to avoid sqrt: distSq < (2^(1/6) * sigma)^2 = 2^(1/3) * sigma^2
inline constexpr double TWO_POW_1_3 = 1.2599210498948731647672106072782283505702514647015079800819751121;

/// 2^(2/3) = (2^(1/3))^2
inline constexpr double TWO_POW_2_3 = 1.5874010519681994747517056392723082603914933278998530098082857618;

}  // namespace PhysicsConstants
