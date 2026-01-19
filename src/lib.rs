//! # Atmospheric and Aeronautical Calculations Library
//! 
//! This library provides functions and calculations based on standard atmospheric models.
//! 
//! Numerous atmospheric models have been published. This crate seeks to replicate the International Standard Atmosphere (ISA) published by the International Civil Aviation Organization (ICAO) via Doc 7488/3.
//! 
//! It is intended that this crate will support other atmospheric models:
//! 
//! - [x] ICAO ISA, Doc 7488/3 (1993)
//! - [ ] US Standard Atmosphere (1976)
//! - [ ] ISO Standard Atmosphere, ISO 2533, (1975)
//! 
//! # Examples
//! 
//! ## Get the ISA Standard Pressure at a Given Geopotential Altitude
//! 
//! ```rust
//! use aero_atmos::InternationalStandardAtmosphere;
//! 
//! // Altitude units
//! use uom::si::length::meter;
//! use uom::si::f64::Length;
//! 
//! // Pressure units
//! use uom::si::pressure::pascal;
//! use uom::si::f64::Pressure;
//! 
//! // Create a new uom Length representing 1,000 meters of Geopotential Altitude. 
//! // Any input unit can be used, so long as it is within the range -5 to 80 km.
//! let altitude = uom::si::f64::Length::new::<uom::si::length::meter>(1_000.0);
//! 
//! // Calculate the ISA standard pressure at the given geopotential altitude.
//! let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude);
//! 
//! match pressure {
//!     Ok(p) => println!("ISA Standard Pressure at 1,000 meters: {} Pa", p.get::<uom::si::pressure::pascal>()),
//!     Err(e) => println!("Error calculating pressure: {:?}", e),
//! }
//! ```
//! 
//! ## Convert between Geopotential and Geometric Altitudes
//! 
//! ```rust
//! use aero_atmos::InternationalStandardAtmosphere;
//! 
//! // Altitude units
//! use uom::si::length::{meter, foot};
//! use uom::si::f64::Length;
//! 
//! // Create a new uom Length representing 10,000 feet of Geometric Altitude.
//! let geometric_altitude = Length::new::<foot>(10_000.0);
//! 
//! // Convert to Geopotential Altitude.
//! // This example shocases the `uom` crate's ability to handle multiple units.
//! let geopotential_altitude = InternationalStandardAtmosphere::altitude_geometric_to_geopotential(geometric_altitude);
//! 
//! match geopotential_altitude {
//!     Ok(alt) => {
//!         println!("Geopotential Altitude at 10,000 feet Geometric: {} meters", alt.get::<meter>());
//!         println!("Geopotential Altitude at 10,000 feet Geometric: {} feet", alt.get::<foot>());
//!     },
//!     Err(e) => println!("Error converting altitude: {:?}", e),
//! }
//!```

pub mod intl_standard_atmos;
pub use intl_standard_atmos::InternationalStandardAtmosphere;

pub mod macros;

// pub use crate::assert_eq_precision;
// pub use self::assert_eq_sigfigs;

/// Precision for floating point comparisons in tests. 
/// 
/// Represents 0.05%.
pub const PRECISION: f64 = 0.0005;
