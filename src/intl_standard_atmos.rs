//! The International Standard Atmosphere (ISA) module.
//! 
//! This module provides functions and data structures to model the 1993 International Standard Atmosphere (ISA) as defined by ICAO Doc 7488/3.

use uom::si::length::*;
use uom::si::pressure::pascal;



/// Struuct for the ISA
pub struct InternationalStandardAtmosphere {

}

/// Error types for ISA calculations
#[derive(Debug, PartialEq)]
pub enum IsaError {
    /// Input value is out of the valid range
    InputOutOfRange,
    ComputationError,
}

impl InternationalStandardAtmosphere {

    /// The standard temperature for a geopotential altitude
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::ThermodynamicTemperature;
    /// use uom::si::length::foot;
    /// use uom::si::thermodynamic_temperature::kelvin;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let altitude = Length::new::<foot>(5_000.0);
    /// let temperature = InternationalStandardAtmosphere::altitude_to_temperature(altitude).unwrap().value;
    /// assert_eq_precision!(temperature, 278.244, PRECISION);
    /// ```
    pub fn altitude_to_temperature (geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::ThermodynamicTemperature, IsaError> {
        // match altitude ranges to determine temperature
        let h = geopotential_altitude.get::<uom::si::length::kilometer>();

        // check if `h` is out of range
        if h < -5.0 || 80.0 < h {
            return Err(IsaError::InputOutOfRange);
        }

        let delta_h:f64;
        let lapse_rate:f64;
        let t0:f64;

        // Match the height ranges to the constants used for the calculation.
        match geopotential_altitude.get::<uom::si::length::kilometer>() {
            h if -5.0 <= h && h <= 11.0 => {
                // Troposphere
                t0 = 288.15; // K
                lapse_rate = -6.5; // K/km
                delta_h = h - 0.0; // km
            },
            h if 11.0 < h && h <= 20.0 => {
                // Lower Stratosphere
                t0 = 216.65; // K
                lapse_rate = 0.0; // K/km
                delta_h = h - 11.0; // km
            },
            h if 20.0 < h && h <= 32.0 => {
                // Middle Stratosphere
                t0 = 216.65; // K
                lapse_rate = 1.0; // K/km
                delta_h = h - 20.0; // km
            },
            h if 32.0 < h && h <= 47.0 => {
                // Upper Stratosphere
                t0 = 228.65; // K
                lapse_rate = 2.8; // K/km
                delta_h = h - 32.0; // km
                
            },
            h if 47.0 < h && h <= 51.0 => {
                // Stratopause
                t0 = 270.65; // K
                lapse_rate = 0.0; // K/km
                delta_h = h - 47.0; // km
            },
            h if 51.0 < h && h <= 71.0 => {
                // Lower Mesosphere
                t0 = 270.65; // K
                lapse_rate = -2.8; // K/km
                delta_h = h - 51.0; // km
            },
            h if 71.0 < h && h <= 80.0 => {
                // Upper Mesosphere
                t0 = 214.65; // K
                lapse_rate = -2.0; // K/km
                delta_h = h - 71.0; // km
            },
            _ => {
                return Err(IsaError::InputOutOfRange);
            }
        }
        Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(t0 + lapse_rate * delta_h))
    }

    /// Base temperature for a geopotential altitude layer
    /// 
    /// THis function returns the base temperature for the layer that contains the given geopotential altitude. If
    /// the given altitude is equal to the bottom of a layer, the base temperature for the next lower layer is returned.
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::ThermodynamicTemperature;
    /// use uom::si::length::kilometer;
    /// use uom::si::thermodynamic_temperature::kelvin;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let altitude = Length::new::<kilometer>(11.0);
    /// let base_temperature = InternationalStandardAtmosphere::altitude_to_base_temperature(altitude).unwrap().value;
    /// // 11km is at the bottom of a layer, so the base temperature of the next lower layer is returned.
    /// assert_eq_precision!(base_temperature, 288.15, PRECISION);
    /// ```
    pub fn altitude_to_base_temperature(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::ThermodynamicTemperature, IsaError> {
        let h = geopotential_altitude.get::<uom::si::length::kilometer>();

        // check if `h` is out of range
        if h < -5.0 || 80.0 < h {
            return Err(IsaError::InputOutOfRange);
        }

        match geopotential_altitude.get::<uom::si::length::kilometer>() {
            h if -5.0 <= h && h <= 11.0 => {
                // Troposphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(288.15))
            },
            h if 11.0 < h && h <= 20.0 => {
                // Lower Stratosphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(216.65))
            },
            h if 20.0 < h && h <= 32.0 => {
                // Middle Stratosphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(216.65))
            },
            h if 32.0 < h && h <= 47.0 => {
                // Upper Stratosphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(228.65))
            },
            h if 47.0 < h && h <= 51.0 => {
                // Stratopause
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(270.65))
            },
            h if 51.0 < h && h <= 71.0 => {
                // Lower Mesosphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(270.65))
            },
            h if 71.0 < h && h <= 80.0 => {
                // Upper Mesosphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(214.65))
            },
            _ => {
                Err(IsaError::InputOutOfRange)
            }
        }
    }

    /// Base geopotential altitude for a geopotential altitude layer
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// If the input H is equal to an H_b, the H_b for the next lower layer is returned.
    /// 
    /// The layers are:
    /// - Troposphere: -5 km to 11 km
    /// - Lower Stratosphere: 11 km to 20 km
    /// - Middle Stratosphere: 20 km to 32 km
    /// - Upper Stratosphere: 32 km to 47 km
    /// - Stratopause: 47 km to 51 km
    /// - Lower Mesosphere: 51 km to 71 km
    /// - Upper Mesosphere: 71 km to 80 km
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::length::kilometer;
    /// 
    /// let altitude = Length::new::<kilometer>(15.0);
    /// let base_altitude = InternationalStandardAtmosphere::altitude_to_base_geopotential_altitude(altitude).unwrap();
    /// assert_eq!(base_altitude.get::<kilometer>(), 11.0);
    /// ```
    pub fn altitude_to_base_geopotential_altitude(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Length, IsaError> {
        let h = geopotential_altitude.get::<uom::si::length::kilometer>();

        // check if `h` is out of range
        if h < -5.0 || 80.0 < h {
            return Err(IsaError::InputOutOfRange);
        }

        match geopotential_altitude.get::<uom::si::length::kilometer>() {
            h if -5.0 <= h && h <= 11.0 => {
                // Troposphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(0.0))
            },
            h if 11.0 < h && h <= 20.0 => {
                // Lower Stratosphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(11.0))
            },
            h if 20.0 < h && h <= 32.0 => {
                // Middle Stratosphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(20.0))
            },
            h if 32.0 < h && h <= 47.0 => {
                // Upper Stratosphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(32.0))
            },
            h if 47.0 < h && h <= 51.0 => {
                // Stratopause
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(47.0))
            },
            h if 51.0 < h && h <= 71.0 => {
                // Lower Mesosphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(51.0))
            },
            h if 71.0 < h && h <= 80.0 => {
                // Upper Mesosphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(71.0))
            },
            _ => {
                Err(IsaError::InputOutOfRange)
            }
        }
    }

    /// Temperature lapse rate for a geopotential altitude layer.
    /// 
    /// Note, this function returns the lapse rate in K/km.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::length::kilometer;
    /// 
    /// let altitude = Length::new::<kilometer>(5.0);
    /// let lapse_rate = InternationalStandardAtmosphere::altitude_to_lapse_rate(altitude).unwrap();
    /// assert_eq!(lapse_rate, -6.5);
    /// ```
    pub fn altitude_to_lapse_rate(geopotential_altitude: uom::si::f64::Length) -> Result<f64, IsaError> {
        let h = geopotential_altitude.get::<uom::si::length::kilometer>();

        // check if `h` is out of range
        if h < -5.0 || 80.0 < h {
            return Err(IsaError::InputOutOfRange);
        }

        // match altitude ranges to determine lapse rate
        match geopotential_altitude.get::<uom::si::length::kilometer>() {
            h if -5.0 <= h && h <= 11.0 => {
                // Troposphere
                Ok(-6.5) // K/km
            },
            h if 11.0 < h && h <= 20.0 => {
                // Lower Stratosphere
                Ok(0.0) // K/km
            },
            h if 20.0 < h && h <= 32.0 => {
                // Middle Stratosphere
                Ok(1.0) // K/km
            },
            h if 32.0 < h && h <= 47.0 => {
                // Upper Stratosphere
                Ok(2.8) // K/km
            },
            h if 47.0 < h && h <= 51.0 => {
                // Stratopause
                Ok(0.0) // K/km
            },
            h if 51.0 < h && h <= 71.0 => {
                // Lower Mesosphere
                Ok(-2.8) // K/km
            },
            h if 71.0 < h && h <= 80.0 => {
                // Upper Mesosphere
                Ok(-2.0) // K/km
            },
            _ => {
                Err(IsaError::InputOutOfRange)
            }
        }
    }

    /// Convert geopotential to geometric altitude
    /// 
    /// This function applies the formula from ICAO Doc 7488/3 section 2.3 to convert H to h.
    /// 
    /// Note that this method is an estimation. The computation that uses lattitude and longitude to determine the precise earth radius at a given location is not implemented in the ISA nor here.
    /// 
    /// # Arguements
    /// `geopotential_altitude` - Geopotential altitude (H) in uom length
    /// 
    /// # Returns
    /// geometric altitude (h) in uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::length::meter;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let geopotential_alt = uom::si::f64::Length::new::<meter>(10_000.0);
    /// let geometric_alt = InternationalStandardAtmosphere::altitude_geopotential_to_geometric(geopotential_alt).unwrap();
    /// assert_eq_precision!(geometric_alt.get::<meter>(), 10015.71, PRECISION);
    /// ```
    pub fn altitude_geopotential_to_geometric(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Length, IsaError> {
        let rad_meter = Self::constant_earth_radius().get::<meter>();
        let geop_alt_meter = geopotential_altitude.get::<meter>();
        let radius_minus_altitude = rad_meter - geop_alt_meter;
        // handle divide by zero
        if radius_minus_altitude == 0.0 {return Err(IsaError::ComputationError)}
        Ok(uom::si::f64::Length::new::<meter>((rad_meter * geop_alt_meter) / radius_minus_altitude))
    }

    /// Convert geometric (h) to geopotential (H) altitude
    /// 
    /// This function applies the formula from ICAO Doc 7488/3 section 2.3 to convert h to H.
    /// 
    /// Note that this method is an estimation. The computation that uses lattitude and longitude to determine the precise earth radius at a given location is not implemented in the ISA nor here.
    /// 
    /// # Arguements
    /// `geometric_altitude` - Geometric altitude (h) in uom length
    /// 
    /// # Returns
    /// Geopotential altitude (H) in uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::length::meter;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let geometric_alt = Length::new::<meter>(10_000.0);
    /// let geopotential_alt = InternationalStandardAtmosphere::altitude_geometric_to_geopotential(geometric_alt).unwrap();
    /// assert_eq_precision!(geopotential_alt.get::<meter>(), 9984.29, PRECISION);
    /// ```
    pub fn altitude_geometric_to_geopotential(geometric_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Length, IsaError> {
        let rad_meter = Self::constant_earth_radius().get::<meter>();
        let geo_alt_meter = geometric_altitude.get::<meter>();
        let radius_plus_altitude = rad_meter + geo_alt_meter;
        // handle divide by zero. Probably impossible, but just in case.
        if radius_plus_altitude == 0.0 {return Err(IsaError::ComputationError)}
        Ok(uom::si::f64::Length::new::<meter>((rad_meter * geo_alt_meter) / radius_plus_altitude))
    }
    
    /// Pressure from geopotential altitude.
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// This function uses an iterative method to compute pressure if the Geopotential altitude is outside of the 
    /// range -5 to 11 km.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::Pressure;
    /// use uom::si::length::kilometer;
    /// use uom::si::pressure::pascal;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let altitude = Length::new::<kilometer>(5.0);
    /// let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap();
    /// assert_eq_precision!(pressure.get::<pascal>(), 54_020.0, PRECISION);
    /// ```
    pub
    fn altitude_to_pressure (geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Pressure, IsaError> {
        let t_base: f64 = Self::altitude_to_base_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // Kelvin
        let h: f64 = geopotential_altitude.get::<uom::si::length::kilometer>(); // km
        let h_base: f64 = Self::altitude_to_base_geopotential_altitude(geopotential_altitude)?.get::<uom::si::length::kilometer>(); // km
        let m: f64 = Self::constant_molar_mass_air().get::<uom::si::molar_mass::gram_per_mole>(); // g/mol
        let g: f64 = Self::constant_gravity_sl().get::<uom::si::acceleration::meter_per_second_squared>(); // m/s^2
        let r_star: f64 = Self::constant_universal_gas_constant(); // J/(kmol·K)
        let beta: f64 = Self::altitude_to_lapse_rate(geopotential_altitude)?; // K/km
        
        let p_base: f64 = match h {
            h if -5.0 <= h && h <= 11.0 => {
                101_325.0 // SL pressure in Pa
            },
            _ => {
                // recursive call to get base pressure for higher layers
                InternationalStandardAtmosphere::altitude_to_pressure(
                    uom::si::f64::Length::new::<uom::si::length::kilometer>(h_base)
                )?.get::<pascal>()
            }
        };

        // Is `H`` out of range?
        if h < -5.0 || h > 80.0 { return Err(IsaError::InputOutOfRange); }

        // match lapse rate to determine pressure
        let p = match beta {
            beta if beta == 0.0 => {
                // isothermal layer
                p_base * (- ( (m * g * (h - h_base)) / (r_star * t_base) ) ).exp()
            },
            _ => {
                // gradient layer
                p_base * (t_base / (t_base + (beta * (h - h_base) )) ).powf( (m * g) / (r_star * beta))
            }
        };

        return Ok(uom::si::f64::Pressure::new::<pascal>(p));

    }

    /// Pressure ratio
    /// 
    /// Returns the pressure ratio (p/P₀) for a given geopotential altitude.
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::{InternationalStandardAtmosphere, IsaError};
    /// use uom::si::f64::Length;
    /// use uom::si::length::kilometer;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let altitude = Length::new::<kilometer>(5.0);
    /// let pressure_ratio = InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    /// assert_eq_precision!(pressure_ratio, 0.5329, PRECISION);
    /// 
    /// // out of range example
    /// let altitude_out_of_range = Length::new::<kilometer>(100.0); // 100 km is too high
    /// let pressure_ratio_result = InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude_out_of_range);
    /// assert_eq!(pressure_ratio_result, Err(IsaError::InputOutOfRange));
    /// ```
    pub fn altitude_to_pressure_ratio(geopotential_altitude: uom::si::f64::Length) -> Result<f64, IsaError> {
        let p = Self::altitude_to_pressure(geopotential_altitude)?.get::<pascal>();
        let p0 = Self::constant_sea_level_standard_pressure().get::<pascal>();
        
        // unlikely to be divide by zero, but just in case
        if p0 == 0.0 {return Err(IsaError::ComputationError);}
        
        Ok(p / p0)
    }

    /// Density from geopotential altitude.
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::MassDensity;
    /// use uom::si::length::foot;
    /// use uom::si::mass_density::kilogram_per_cubic_meter;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // H = 14,000 ft => 7.96523e-1 kg/m³ (Doc 7488/3 Table 4)
    /// let altitude = uom::si::f64::Length::new::<foot>(14_000.0);
    /// let density = InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap();
    /// assert_eq_precision!(density.get::<kilogram_per_cubic_meter>(), 7.96523e-1, PRECISION);
    /// ```
    pub fn altitude_to_density (geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::MassDensity, IsaError> {
        let p = Self::altitude_to_pressure(geopotential_altitude)?.get::<pascal>(); // Pa
        let t = Self::altitude_to_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K
        let r_specific = Self::constant_specific_gas_constant_air(); // J/(kg·K)

        // handle divide by zero
        let divisor = r_specific * t;
        if divisor == 0.0 {return Err(IsaError::ComputationError);}

        let rho = p / divisor; // kg/m^3

        Ok(uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(rho))
    }

    /// Altitude to Density ratio
    /// 
    /// Returns the density ratio (ρ/ρ₀) for a given geopotential altitude.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::{InternationalStandardAtmosphere, IsaError};
    /// use uom::si::f64::Length;
    /// use uom::si::length::meter;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 2: H = 5,000 meter => ρ/ρ₀ = 6.00911e-1
    /// let altitude = Length::new::<meter>(5_000.0);
    /// let density_ratio = InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
    /// assert_eq_precision!(density_ratio, 6.00911e-1, PRECISION);
    /// ```
    pub fn altitude_to_density_ratio(geopotential_altitude: uom::si::f64::Length) -> Result<f64, IsaError> {
        let rho = Self::altitude_to_density(geopotential_altitude)?.get::<uom::si::mass_density::kilogram_per_cubic_meter>();
        let rho_0 = Self::altitude_to_density(uom::si::f64::Length::new::<uom::si::length::meter>(0.0))?.get::<uom::si::mass_density::kilogram_per_cubic_meter>(); // kg/m^3 at sea level

        // unlikely to be divide by zero, but just in case
        if rho_0 == 0.0 {return Err(IsaError::ComputationError);}

        Ok(rho / rho_0)
    }

    /// Square root of density ratio
    /// 
    /// Returns the square root of the density ratio (√(ρ/ρ₀)) for a given geopotential altitude.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::{InternationalStandardAtmosphere, IsaError};
    /// use uom::si::f64::Length;
    /// use uom::si::length::meter;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 2: H = 5,000 meter => √(ρ/ρ₀) = 7.75184e-1
    /// let altitude = Length::new::<meter>(5_000.0);
    /// let sqrt_density_ratio = InternationalStandardAtmosphere::altitude_to_sqrt_density_ratio(altitude).unwrap();
    /// assert_eq_precision!(sqrt_density_ratio, 7.75184e-1, PRECISION);
    /// ```
    pub fn altitude_to_sqrt_density_ratio(geopotential_altitude: uom::si::f64::Length) -> Result<f64, IsaError> {
        let density_ratio = Self::altitude_to_density_ratio(geopotential_altitude)?;
        Ok(density_ratio.sqrt())
    }

    /// Speed of sound from geopotential altitude.
    /// 
    /// The symbol for speed of sound is `a`.
    /// 
    /// "This expression ... presents the speed of
    /// propagation of an infinitesimal perturbation in a gas. That
    /// is why this formula may not be used for calculation, for
    /// example, of the speed of propagation of shock waves
    /// induced by blast, detonation, body motion in the air at
    /// supersonic speed, etc." (ICAO Doc 7488/3, section 2.14)
    /// 
    /// Tabular values for speed of sound are given in ICAO Doc 7488/3 Table 5.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::Velocity;
    /// use uom::si::length::foot;
    /// use uom::si::velocity::meter_per_second;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 5: 
    /// // H = 16_000.0 ft => a = 321.031 m/s
    /// let altitude = Length::new::<foot>(16_000.0);
    /// let speed_of_sound = InternationalStandardAtmosphere::altitude_to_speed_of_sound(altitude).unwrap();
    /// assert_eq_precision!(speed_of_sound.get::<meter_per_second>(), 321.031, PRECISION);
    /// ```
    pub fn altitude_to_speed_of_sound(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Velocity, IsaError> {
        let gamma: f64 = Self::constant_ratio_specific_heats_air(); // dimensionless
        let r_specific: f64 = Self::constant_specific_gas_constant_air(); // J/(kg·K)
        let t: f64 = Self::altitude_to_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K

        let a = (gamma * r_specific * t).sqrt(); // m/s

        Ok(uom::si::f64::Velocity::new::<uom::si::velocity::meter_per_second>(a))
    }

    /// Dynamic viscosity from geopotential altitude.
    /// 
    /// The symbol for dynamic viscosity is `μ`.
    /// 
    /// Tabular values for dynamic viscosity are given in ICAO Doc 7488/3 Table 5.
    /// The table lists `H` in feet.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::DynamicViscosity;
    /// use uom::si::length::foot;
    /// use uom::si::dynamic_viscosity::pascal_second;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 5: 
    /// // H = 16_000.0 ft => μ = 1.6322e-5 Pa·s
    /// let altitude = Length::new::<foot>(16_000.0);
    /// let dynamic_viscosity = InternationalStandardAtmosphere::altitude_to_dynamic_viscosity(altitude).unwrap();
    /// assert_eq_precision!(dynamic_viscosity.get::<pascal_second>(), 1.6322e-5, PRECISION);
    /// ```
    pub fn altitude_to_dynamic_viscosity(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::DynamicViscosity, IsaError> {
        let t: f64 = Self::altitude_to_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K
        let beta_s: f64 = Self::constant_sutherland_beta_s(); // kg/(m·s·K^0.5)
        let s: f64 = Self::constant_sutherland_s().get::<uom::si::thermodynamic_temperature::kelvin>(); // K

        let mu = (beta_s * t.powf(3.0/2.0)) / (t + s);

        Ok(uom::si::f64::DynamicViscosity::new::<uom::si::dynamic_viscosity::pascal_second>(mu))
    }

    /// The earth's radius.
    /// 
    /// This function returns the constant earth radius used in the ISA model.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::length::kilometer;
    ///
    /// let earth_radius = InternationalStandardAtmosphere::constant_earth_radius();
    /// assert_eq!(earth_radius.get::<kilometer>(), 6_356.766);
    /// ```
    pub fn constant_earth_radius() -> uom::si::f64::Length {
        return uom::si::f64::Length::new::<uom::si::length::meter>(6_356_766.0);
    }

    /// Standard Gravity constant at sea level.
    /// 
    /// Returns the standard gravity constant used in the ISA model.
    /// 
    /// The symbol is `G`.
    /// 
    /// "standard acceleration due to gravity. It conforms 
    /// with latitude n = 45°32'33"N using
    /// Lambert’s equation of the acceleration due
    /// to gravity as a function of latitude" (ICAO Doc 7488/3, section 2.2)
    pub fn constant_gravity_sl() -> uom::si::f64::Acceleration {
        uom::si::f64::Acceleration::new::<uom::si::acceleration::meter_per_second_squared>(9.806_65)
    }

    /// Molar mass of Earth's air.
    /// 
    /// Returns the molar mass of Earth's air used in the ISA model.
    /// 
    /// The symbol is `M_0`.
    pub fn constant_molar_mass_air() -> uom::si::f64::MolarMass {
        uom::si::f64::MolarMass::new::<uom::si::molar_mass::gram_per_mole>(28_964.420)
    }

    /// Avogadro's number.
    /// 
    /// Returns Avogadro's number used in the ISA model.
    pub fn constant_avogadro_number() -> f64 {
        602.257e24 // kmol^-1
    }

    /// Specific gas constant for air.
    /// 
    /// Returns the specific gas constant for air used in the ISA model.
    pub fn constant_specific_gas_constant_air() -> f64 {
        287.052_870 // J/(kg·K)
    }

    /// Universal gas constant.
    /// 
    /// Returns the universal gas constant used in the ISA model in J/(kmol·K).
    /// 
    /// The symbol is `R*`.
    pub fn constant_universal_gas_constant() -> f64 {
        8_314.32 // J/(kmol·K)
    }

    /// Sutherland constant.
    /// 
    /// Returns the Sutherland constant used in the ISA model.
    /// 
    /// The symbol is `S`.
    pub fn constant_sutherland_s() -> uom::si::f64::ThermodynamicTemperature {
        uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(110.4)
    }

    /// Sutherland's constant (β_s)
    /// 
    /// Returns Sutherland's constant (β_s) used in the ISA model.
    /// 
    /// The symbol is `β_s`.
    pub fn constant_sutherland_beta_s() -> f64 {
        1.458e-6 // kg/(m·s·K^0.5)
    }

    /// Ratio of specific heats for air.
    /// 
    /// Returns the ratio of specific heats for air used in the ISA model.
    /// 
    /// The symbol is `γ`.
    pub fn constant_ratio_specific_heats_air() -> f64 {
        1.4
    }

    /// effective collision diameter of an air molecule
    /// 
    /// Returns the effective collision diameter of an air molecule used in the ISA model.
    /// 
    /// The symbol is `σ` (sigma).
    pub fn constant_effective_collision_diameter_air_molecule() -> uom::si::f64::Length {
        uom::si::f64::Length::new::<uom::si::length::meter>(3.621e-10)
    }

    /// Sea level standard pressure.
    /// 
    /// Returns the sea level standard pressure used in the ISA model.
    /// 
    /// The symbol is `P_0`.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Pressure;
    /// use uom::si::pressure::hectopascal;
    ///
    /// let sl_pressure = InternationalStandardAtmosphere::constant_sea_level_standard_pressure();
    /// assert_eq!(sl_pressure.get::<hectopascal>(), 1013.25);
    /// ```
    pub fn constant_sea_level_standard_pressure() -> uom::si::f64::Pressure {
        uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(101_325.0)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::length::kilometer;
    // use uom::si::thermodynamic_temperature::kelvin;
    use uom::si::pressure::pascal;
    use crate::PRECISION;
    use crate::assert_eq_precision;

    /// Test altitude to temperature conversion
    #[test]
    fn test_altitude_to_pressure() {
        let input_alt: f64 = 11.466;
        let desired_pressure = 210.3 * 100.0; // Pa (Doc 7488/3 Table 7 gives hPa)
        let altitude = uom::si::f64::Length::new::<kilometer>(input_alt);
        let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
        println!("Pressure at {} km: Desired {:.1} Pa, Calculated {:.1} Pa", input_alt, desired_pressure, pressure);
        assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<pascal>(desired_pressure).value, PRECISION);

        let input_alt:f64 = 31.985; // kilometers
        let desired_pressure = 8.70 * 100.0; // Pa
        let altitude = uom::si::f64::Length::new::<kilometer>(input_alt);
        let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
        println!("Pressure at {} km: Desired {:.2} Pa, Calculated {:.2} Pa", input_alt, desired_pressure, pressure);
        assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<pascal>(desired_pressure).value, PRECISION);
    }

}