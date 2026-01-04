
// pub mod us_standard_atmosphere;
// pub use us_standard_atmosphere::*;

pub mod intl_standard_atmos;
pub use intl_standard_atmos::InternationalStandardAtmosphere;

pub mod macros;
pub use macros::*;

/// Precision for floating point comparisons in tests. 
/// 
/// Represents 0.05%.
pub const PRECISION: f64 = 0.0005;
