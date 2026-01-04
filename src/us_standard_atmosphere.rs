//! This module provides functions and data structures to model the 1976 US Standard Atmosphere.

use uom::si::f64::{Pressure, Temperature, Length, Density, Velocity};

/// The Boltzmann constant.
/// 
/// Symbol: `k`
/// 
/// Value: `1.380649 × 10⁻²³ N*m/K`
pub const BOLTZMANN_CONSTANT: f64 = 1.380649e-23; // J/K

/// Represents the US Standard Atmosphere model.
pub struct UsStandardAtmosphere {
}

