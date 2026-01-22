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
//! // uom Altitude units
//! use uom::si::length::meter;
//! use uom::si::f64::Length;
//! 
//! // uom Pressure units
//! use uom::si::pressure::pascal;
//! use uom::si::f64::Pressure;
//! 
//! // Create a new uom Length representing 16,000 meters of Geopotential Altitude. 
//! // Any input unit can be used, so long as it is within the range -5 to 80 km.
//! let altitude = uom::si::f64::Length::new::<uom::si::length::meter>(16_000.0);
//! 
//! // Calculate the ISA standard pressure at the given geopotential altitude.
//! let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude);
//! 
//! match pressure {
//!     // prints ~ "p = 101325.0 Pa"
//!     Ok(p) => println!("p = {} Pa", p.get::<uom::si::pressure::pascal>()),
//!     // should not print an error message, since the altitude is within the valid range for the ISA model
//!     Err(e) => println!("Error calculating pressure: {:?}", e),
//! }
//! ```
//! 
//! ## Convert between Geopotential and Geometric Altitudes
//! 
//! The ISA uses an appoximation formula to convert between geopotential and geometric altitudes.
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
//!         // prints something like "Geopotential Altitude at 10,000 feet Geometric: 9995.0 feet"
//!         println!("Geopotential Altitude at 10,000 feet Geometric: {} feet", alt.get::<foot>());
//!     },
//!     Err(e) => println!("Error converting altitude: {:?}", e),
//! }
//!```
//! 
//! ## Get the Geometric altitude for a given pressure.
//! 
//! ```rust
//! use aero_atmos::InternationalStandardAtmosphere;
//! use aero_atmos::dbg_sci;
//! 
//! // Altitude units
//! use uom::si::length::{meter, foot};
//! use uom::si::f64::Length;
//! 
//! // Pressure units
//! use uom::si::pressure::pascal;
//! use uom::si::f64::Pressure;
//! 
//! // Create a new uom Pressure representing 50,000 Pa.
//! let pressure = Pressure::new::<pascal>(50_000.0);
//! 
//! // First, we calculate the Geopotential Altitude for the given pressure.
//! // In Doc 7488/3 Table 7, 50_000 pa is about 5_574.0 meters.
//! let geopotential_altitude = InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap();
//! dbg_sci!(geopotential_altitude.get::<meter>(), 4);
//! // prints something like "[src/main.rs:31] geopotential_altitude.get::<meter>() = 5.5744e3"
//! 
//! // Then, we convert the Geopotential Altitude to Geometric Altitude.
//! let geometric_altitude = InternationalStandardAtmosphere::altitude_geopotential_to_geometric(geopotential_altitude).unwrap();
//! dbg_sci!(geometric_altitude.get::<meter>(), 4);
//! // prints something like "[src/main.rs:36] geometric_altitude.get::<meter>() = 5.5793e3"
//! ```

pub mod intl_standard_atmos;
pub use intl_standard_atmos::InternationalStandardAtmosphere;

pub mod macros;

/// Precision for testing. 0.05%. 
/// 
/// Most or all tests comnpare the output of these functions to the
/// tabular data in the tables of the ICAO Doc 7488/3 (1993). Most or all
/// functions are tested to this precision.
pub const PRECISION: f64 = 0.0005;
