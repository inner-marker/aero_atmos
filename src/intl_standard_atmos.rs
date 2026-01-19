//! The International Standard Atmosphere (ISA) module.
//! 
//! This module provides functions and data structures to model the 1993 International Standard Atmosphere (ISA) as defined by ICAO Doc 7488/3.

/// Error types for ISA calculations
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum IsaError {
    /// Input value is out of the valid range
    InputOutOfRange,

    /// A computation error occurred (e.g., division by zero)
    ComputationError,
}

/// Struuct for the International Standard Atmosphere
/// 
/// This struct contains methods to compute various atmospheric properties based on the International Standard Atmosphere (ISA) model.
/// 
/// The source is ICAO Doc 7488/3, "Manual of the ICAO Standard Atmosphere", Third Edition, 1993.
/// 
/// Computations are based on the values and equations provided in the document. Some more accurate or precise 
/// methods and values may exist elsewhere, but this implementation aims to closely follow the ISA specification.
/// 
/// For all methods where physical quantities are involved, the `uom` crate is used to ensure type safety and unit correctness.
/// Some values are unitless, so are represented as `f64`.
pub struct InternationalStandardAtmosphere {
}

impl InternationalStandardAtmosphere {

    /// List the Base Altitudes
    /// 
    /// This function returns a list of the ISA base altitudes for layers as defined in Doc 7488/3, Table D.
    ///
    /// # Returns
    /// `Vec<uom::si::f64::Length>`
    pub fn list_base_altitudes () -> Vec<uom::si::f64::Length> {
        let out: Vec<uom::si::f64::Length> = vec![
            uom::si::f64::Length::new::<uom::si::length::kilometer>(-5.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(0.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(11.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(20.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(32.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(47.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(51.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(71.0),
            uom::si::f64::Length::new::<uom::si::length::kilometer>(80.0),
        ];

        out
    }

    /// The standard temperature for a geopotential altitude
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Returns
    /// Temperature in uom thermodynamic temperature.
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range (-5 to 80 km).
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
            h if -5.0 <= h && h <= 0.0 => {
                // Troposphere
                Ok(uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(320.65))
            },
            h if 0.0 < h && h <= 11.0 => {
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
    /// # Returns
    /// Base geopotential altitude for the layer in uom length.
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range (-5 to 80 km).
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

        match h {
            h if -5.0 <= h && h <= 0.0 => {
                // Troposphere
                Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(-5.0))
            },
            h if 0.0 < h && h <= 11.0 => {
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
    /// Note, this function returns the lapse rate in K/km. Also, if the input geopotential
    /// altitude is exactly equal to the base altitude of a layer, the lapse rate for the next lower layer is returned.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length.
    /// 
    /// # Returns
    /// Lapse rate in K/km as f64.
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range (-5 to 80 km).
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
        if h < -5.0 || 80.0 < h { return Err(IsaError::InputOutOfRange); }

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
    /// # Errors
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range (-5 to 80 km).
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
        let geopotential_alt_meter = geopotential_altitude.get::<uom::si::length::meter>();

        // check for valid range
        if geopotential_alt_meter < -5.0 || geopotential_alt_meter > 80_000.0 { return Err(IsaError::InputOutOfRange); }

        let rad_meter = Self::constant_earth_radius().get::<uom::si::length::meter>();
        let radius_minus_altitude = rad_meter - geopotential_alt_meter;

        // handle divide by zero
        if radius_minus_altitude == 0.0 {return Err(IsaError::ComputationError)}

        Ok(uom::si::f64::Length::new::<uom::si::length::meter>((rad_meter * geopotential_alt_meter) / radius_minus_altitude))
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
    /// # Errors
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// - Returns `IsaError::InputOutOfRange` if the geometric altitude is negative.
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
        let geometric_alt_meter = geometric_altitude.get::<uom::si::length::meter>();
        
        // check for valid range 
        if geometric_alt_meter < -5.0 || geometric_alt_meter > 80_000.0 { return Err(IsaError::InputOutOfRange); }
        
        let rad_meter = Self::constant_earth_radius().get::<uom::si::length::meter>();
        let radius_plus_altitude = rad_meter + geometric_alt_meter;

        // handle divide by zero. Probably impossible, but just in case.
        if radius_plus_altitude == 0.0 {return Err(IsaError::ComputationError)}
        
        Ok(uom::si::f64::Length::new::<uom::si::length::meter>((rad_meter * geometric_alt_meter) / radius_plus_altitude))
    }
    
    /// Pressure from geopotential altitude.
    /// 
    /// Tabular values for pressure are given in ICAO Doc 7488/3 Table 4.
    /// In Table 4, values for `H` are listed in feet, and pressure in pascals.
    /// 
    /// This function uses a recursive method to compute pressure if the Geopotential altitude is outside of the 
    /// lowest layer range of -5 to 11 km.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range (-5 to 80 km).
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::Pressure;
    /// use uom::si::length::{kilometer, foot};
    /// use uom::si::pressure::{pascal, hectopascal};
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let altitude = Length::new::<kilometer>(5.0);
    /// let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap();
    /// assert_eq_precision!(pressure.get::<pascal>(), 54_020.0, PRECISION);
    /// 
    /// let altitude = Length::new::<foot>(-16_250.0);
    /// let pressure = InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap();
    /// assert_eq_precision!(pressure.get::<hectopascal>(), 1.76799e3, PRECISION);
    /// ```
    pub
    fn altitude_to_pressure (geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Pressure, IsaError> {
        let t_base: f64 = Self::altitude_to_base_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // Kelvin
        let h: f64 = geopotential_altitude.get::<uom::si::length::kilometer>(); // km
        let h_base: f64 = Self::altitude_to_base_geopotential_altitude(geopotential_altitude)?.get::<uom::si::length::kilometer>(); // km
        let m: f64 = Self::constant_molar_mass_dry_air().get::<uom::si::molar_mass::gram_per_mole>(); // g/mol
        let g: f64 = Self::constant_gravity_sealevel().get::<uom::si::acceleration::meter_per_second_squared>(); // m/s^2
        let r_star: f64 = Self::constant_universal_gas_constant(); // J/(kmol·K)
        let beta: f64 = Self::altitude_to_lapse_rate(geopotential_altitude)?; // K/km
        
        
        let p_base: f64 = match h {
            h if -5.0 <= h && h <= 0.0 => {
                1.75364e0 * 101_325.0 // -5 km pressure in Pa
            },
            h if 0.0 <= h && h <= 11.0 => {
                101_325.0 // SL pressure in Pa
            },
            _ => {
                // recursive call to get base pressure for higher layers
                InternationalStandardAtmosphere::altitude_to_pressure(
                    uom::si::f64::Length::new::<uom::si::length::kilometer>(h_base)
                )?.get::<uom::si::pressure::pascal>()
            }
        };
        
        // Is `H`` out of range?
        if h < -5.0 || h > 80.0 { return Err(IsaError::InputOutOfRange); }
        
        // dbg!(t_base, h, h_base, p_base, m, g, r_star, beta);

        // match lapse rate to determine pressure
        let p = match beta {
            beta if beta == 0.0 => {
                // isothermal layer
                p_base * (- ( (m * g * (h - h_base)) / (r_star * t_base) ) ).exp()
            },
            _ => {
                // non-isothermal layer
                p_base * (t_base / (t_base + (beta * (h - h_base) )) ).powf( (m * g) / (r_star * beta))
            }
        };

        return Ok(uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(p));

    }

    /// Pressure ratio
    /// 
    /// Returns the pressure ratio (p/P₀) for a given geopotential altitude.
    /// 
    /// Tabular values for p/P₀ are given in ICAO Doc 7488/3 Table 2.
    /// In Table 2, values for `H` are listed in meters.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range (-5 to 80 km).
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
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
        let p = Self::altitude_to_pressure(geopotential_altitude)?.get::<uom::si::pressure::pascal>();
        let p0 = Self::constant_sea_level_standard_pressure().get::<uom::si::pressure::pascal>();
        
        // unlikely to be divide by zero, but just in case
        if p0 == 0.0 {return Err(IsaError::ComputationError);}
        
        Ok(p / p0)
    }

    /// Density from geopotential altitude.
    /// 
    /// The valid range for geopotential altitude is -5 km to 80 km.
    /// 
    /// Tabular values for density are given in ICAO Doc 7488/3 Table 4.
    /// In Table 4, values for `H` are listed in feet.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range.
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
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
        let p = Self::altitude_to_pressure(geopotential_altitude)?.get::<uom::si::pressure::pascal>(); // Pa
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
    /// Tabular values for ρ/ρ₀ are given in ICAO Doc 7488/3 Table 5.
    /// Values for `H` are listed in feet.
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range.
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero). 
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
    /// Tabular values for √(ρ/ρ₀) are given in ICAO Doc 7488/3 Table 2.
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range.
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// 
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

    /// Kinematic viscosity from geopotential altitude.
    /// 
    /// The symbol for kinematic viscosity is `ν` (nu).
    /// 
    /// Tabular values for kinematic viscosity are given in ICAO Doc 7488/3 Table 5.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::KinematicViscosity;
    /// use uom::si::length::foot;
    /// use uom::si::kinematic_viscosity::square_meter_per_second;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 5: 
    /// // H = 16_000.0 ft => ν = 2.1880e-5 m²/s
    /// let altitude = Length::new::<foot>(16_000.0);
    /// let kinematic_viscosity = InternationalStandardAtmosphere::altitude_to_kinematic_viscosity(altitude).unwrap();
    /// assert_eq_precision!(kinematic_viscosity.get::<square_meter_per_second>(), 2.1880e-5, PRECISION);
    /// ```
    pub fn altitude_to_kinematic_viscosity(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::KinematicViscosity, IsaError> {
        let mu = Self::altitude_to_dynamic_viscosity(geopotential_altitude)?.get::<uom::si::dynamic_viscosity::pascal_second>(); // Pa·s
        let rho = Self::altitude_to_density(geopotential_altitude)?.get::<uom::si::mass_density::kilogram_per_cubic_meter>(); // kg/m^3

        // handle divide by zero
        if rho == 0.0 {return Err(IsaError::ComputationError);}

        let nu = mu / rho; // m^2/s

        Ok(uom::si::f64::KinematicViscosity::new::<uom::si::kinematic_viscosity::square_meter_per_second>(nu))
    }

    /// Thermal conductivity from geopotential altitude.
    /// 
    /// The symbol for thermal conductivity is `λ` (lambda).
    /// 
    /// Tabular values for thermal conductivity are given in ICAO Doc 7488/3 Table 5.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::ThermalConductivity;
    /// use uom::si::length::foot;
    /// use uom::si::thermal_conductivity::watt_per_meter_kelvin;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 5: 
    /// // H = 16_000.0 ft => λ = 2.2810e-2 W/(m·K)
    /// let altitude = Length::new::<foot>(16_000.0);
    /// let thermal_conductivity = InternationalStandardAtmosphere::altitude_to_thermal_conductivity(altitude).unwrap();
    /// assert_eq_precision!(thermal_conductivity.get::<watt_per_meter_kelvin>(), 2.2810e-2, PRECISION);
    /// ```
    pub fn altitude_to_thermal_conductivity(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::ThermalConductivity, IsaError> {
        let t: f64 = Self::altitude_to_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K

        let dividend = 2.648_151e-3 * t.powf(3.0/2.0);
        let divisor = t + (245.4 * 10.0_f64.powf(-(12.0/t)));

        // handle divide by zero, unlikely but just in case
        if divisor == 0.0 {return Err(IsaError::ComputationError);}

        let lambda = dividend / divisor; // W/(m·K)

        Ok(uom::si::f64::ThermalConductivity::new::<uom::si::thermal_conductivity::watt_per_meter_kelvin>(lambda))
    }

    /// Number Density from geopotential altitude.
    /// 
    /// The symbol for number density is _n_. 
    /// 
    /// The unit is m⁻³ (molecules per cubic meter). Therefore, `uom` is not used here.
    /// 
    /// 
    /// Tabular values for number density are given in ICAO Doc 7488/3 Table 3.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 3: 
    /// // H = 3.0 km => n = 1.8903e25 m⁻³
    /// let altitude = Length::new::<uom::si::length::kilometer>(3.0);
    /// let number_density = InternationalStandardAtmosphere::altitude_to_number_density(altitude).unwrap();
    /// assert_eq_precision!(number_density, 1.8903e25, PRECISION);
    /// ```
    pub fn altitude_to_number_density(geopotential_altitude: uom::si::f64::Length) -> Result<f64, IsaError> {
        let p = Self::altitude_to_pressure(geopotential_altitude)?.get::<uom::si::pressure::pascal>(); // Pa
        let t = Self::altitude_to_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K
        let r_star = Self::constant_universal_gas_constant(); // J/(kmol·K)
        let n_a = Self::constant_avogadro_number(); // kmol^-1

        // handle divide by zero
        let divisor = r_star * t;
        if divisor == 0.0 {return Err(IsaError::ComputationError);}

        let n = (p * n_a) / divisor; // m^-3

        Ok(n)
    }

    /// Mean particle speed
    /// 
    /// Returns the mean particle speed ($\bar{v}$) for a given geopotential altitude.
    /// 
    /// Tabular values for mean particle speed are given in ICAO Doc 7488/3 Table 3.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Errors
    /// Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range of -5 km to 80 km.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::Velocity;
    /// use uom::si::length::meter;
    /// use uom::si::velocity::meter_per_second;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 3: 
    /// // H = 3_000.0 m => v̄ = 443.14 m/s
    /// let altitude = Length::new::<meter>(3_000.0);
    /// let mean_particle_speed = InternationalStandardAtmosphere::altitude_to_mean_particle_speed(altitude).unwrap();
    /// assert_eq_precision!(mean_particle_speed.get::<meter_per_second>(), 443.14, PRECISION);
    /// ```
    pub fn altitude_to_mean_particle_speed(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Velocity, IsaError> {
        let t: f64 = Self::altitude_to_temperature(geopotential_altitude)?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K
        let r: f64 = Self::constant_specific_gas_constant_air(); // J/(kmol·K)

        let constant = 1.595_769;

        let v_bar = constant * (r * t).sqrt(); // m/s

        Ok(uom::si::f64::Velocity::new::<uom::si::velocity::meter_per_second>(v_bar))
    }

    /// Mean Free Path
    /// 
    /// Returns the mean free path (_`l`_) for a given geopotential altitude.
    /// 
    /// Tabular values for mean free path are given in ICAO Doc 7488/3 Table 3.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Returns
    /// Mean free path in uom length
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range of -5 km to 80 km.
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
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
    /// // ICAO Doc 7488/3 Table 3: 
    /// // H = 3_000.0 m => l = 8.9374e-8 m
    /// let altitude = Length::new::<meter>(3_000.0);
    /// let mean_free_path = InternationalStandardAtmosphere::altitude_to_mean_free_path(altitude).unwrap();
    /// assert_eq_precision!(mean_free_path.get::<meter>(), 8.9374e-8, PRECISION);
    /// ```
    pub fn altitude_to_mean_free_path(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Length, IsaError> {
        let pi = std::f64::consts::PI;
        let sigma = Self::constant_effective_collision_diameter_air_molecule().get::<uom::si::length::meter>(); // m
        let n = Self::altitude_to_number_density(geopotential_altitude)?; // m^-3

        let divisor = (2.0_f64).sqrt() * pi * sigma.powf(2.0_f64) * n;
        // handle divide by zero
        if divisor == 0.0 { return Err(IsaError::ComputationError)}

        Ok(uom::si::f64::Length::new::<uom::si::length::meter>(1.0_f64 / divisor))
    }

    /// Collision Frequency
    /// 
    /// Returns the collision frequency (_`ω`_) for a given geopotential altitude.
    /// 
    /// Tabular values for collision frequency are given in ICAO Doc 7488/3 Table 3.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Errors
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range of -5 km to 80 km.
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::Frequency;
    /// use uom::si::length::meter;
    /// use uom::si::frequency::hertz;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 3: 
    /// // H = 3_000.0 m => ω = 4.9583e9 s⁻¹
    /// let altitude = Length::new::<meter>(3_000.0);
    /// let collision_frequency = InternationalStandardAtmosphere::altitude_to_collision_frequency(altitude).unwrap();
    /// assert_eq_precision!(collision_frequency.get::<hertz>(), 4.9583e9, PRECISION);
    /// ```
    pub fn altitude_to_collision_frequency(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Frequency, IsaError> {
        let v_bar = Self::altitude_to_mean_particle_speed(geopotential_altitude)?.get::<uom::si::velocity::meter_per_second>(); // m/s
        let l = Self::altitude_to_mean_free_path(geopotential_altitude)?.get::<uom::si::length::meter>(); // m

        // handle divide by zero
        if l == 0.0 {return Err(IsaError::ComputationError);}

        let omega = v_bar / l; // s^-1

        Ok(uom::si::f64::Frequency::new::<uom::si::frequency::hertz>(omega))
    }

    /// Gravitational Acceleration from geopotential altitude.
    /// 
    /// The symbol for gravitational acceleration is _`g`_.
    /// 
    /// Tabular values for gravitational acceleration are given in ICAO Doc 7488/3 Table 4.
    /// Geopotential altitude in Table 4 is given in feet.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Errors
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range of -5 km to 80 km.
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::f64::Acceleration;
    /// use uom::si::length::foot;
    /// use uom::si::acceleration::meter_per_second_squared;
    /// use aero_atmos::assert_eq_precision;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// // ICAO Doc 7488/3 Table 4: 
    /// // H = 16_000.0 ft => g = 9.7916 m/s²
    /// let altitude = Length::new::<foot>(16_000.0);
    /// let gravitational_acceleration = InternationalStandardAtmosphere::altitude_to_gravitational_acceleration(altitude).unwrap();
    /// assert_eq_precision!(gravitational_acceleration.get::<meter_per_second_squared>(), 9.7916, PRECISION);
    /// ```
    pub fn altitude_to_gravitational_acceleration(geopotential_altitude: uom::si::f64::Length) -> Result<uom::si::f64::Acceleration, IsaError> {
        let h = geopotential_altitude.get::<uom::si::length::meter>(); // m
        
        // Is `H` out of range?
        if h < -5_000.0 || h > 80_000.0 { return Err(IsaError::InputOutOfRange); }

        let g0 = Self::constant_gravity_sealevel().get::<uom::si::acceleration::meter_per_second_squared>(); // m/s^2
        let r = Self::constant_earth_radius().get::<uom::si::length::meter>(); // m

        let divisor = r + h;
        // handle divide by zero
        if divisor == 0.0 {return Err(IsaError::ComputationError);}

        let g = g0 * (r / divisor).powf(2.0);

        Ok(uom::si::f64::Acceleration::new::<uom::si::acceleration::meter_per_second_squared>(g))
    }

    /// Specific Weight from geopotential altitude.
    /// 
    /// The symbol for specific weight is _`γ`_ (gamma). The unit is N/m³.
    /// 
    /// Tabular values for specific weight are given in ICAO Doc 7488/3 Table 3.
    /// Values for geopotential altitude are given in meters.
    /// 
    /// # Arguments
    /// * `geopotential_altitude` - Geopotential altitude using uom length
    /// 
    /// # Returns
    /// Specific weight in N/m³ as `f64`. NOTE: uom has no quantity for specific weight, so it is not used.
    /// 
    /// # Errors
    /// - Returns `IsaError::ComputationError` if a computation error occurs (e.g., division by zero).
    /// - Returns `IsaError::InputOutOfRange` if the geopotential altitude is outside the valid range of -5 km to 80 km.
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
    /// // ICAO Doc 7488/3 Table 3: 
    /// // H = 3_000.0 m => γ = 8.9070 N/m³
    /// let altitude = Length::new::<meter>(3_000.0);
    /// let specific_weight = InternationalStandardAtmosphere::altitude_to_specific_weight(altitude).unwrap();
    /// assert_eq_precision!(specific_weight, 8.9070, PRECISION);
    /// ```
    pub fn altitude_to_specific_weight(geopotential_altitude: uom::si::f64::Length) -> Result<f64, IsaError> {
        let rho = Self::altitude_to_density(geopotential_altitude)?.get::<uom::si::mass_density::kilogram_per_cubic_meter>(); // kg/m^3
        let g = Self::altitude_to_gravitational_acceleration(geopotential_altitude)?.get::<uom::si::acceleration::meter_per_second_squared>(); // m/s^2

        let specific_weight = rho * g; // N/m^3

        Ok(specific_weight)
    }

    /// Geopotential Altitude(s) from Pressure
    /// 
    /// Given a pressure, this function returns the corresponding geopotential altitude(s) based on the ISA model.
    /// 
    /// # Arguments
    /// * `pressure` - Pressure using uom pressure
    /// 
    /// # Returns
    /// `Result<uom::si::f64::Length, IsaError>`
    /// 
    /// # Examples
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use aero_atmos::{assert_eq_sigfigs, assert_eq_precision};
    /// use uom::si::f64::Pressure;
    /// use uom::si::pressure::hectopascal;
    /// use uom::si::f64::Length;
    /// use uom::si::length::foot;
    /// 
    /// const PRECISION: f64 = 0.0005; // 0.05%
    /// 
    /// let pressure = Pressure::new::<hectopascal>(1.01325e3);
    /// let altitude = InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap().get::<uom::si::length::foot>();
    /// let altitude_check = 0.0; // feet
    /// assert_eq_precision!(altitude, altitude_check, PRECISION);
    /// 
    /// let pressure = Pressure::new::<hectopascal>(5.491_52e2);
    /// let altitude = InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap().get::<uom::si::length::foot>();
    /// let altitude_check = 16_000.0; // feet
    /// assert_eq_precision!(altitude, altitude_check, PRECISION);
    /// 
    /// 
    /// let pressure = Pressure::new::<hectopascal>(9.08455e-3);
    /// let altitude = InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap().get::<uom::si::length::foot>();
    /// let altitude_check = 262_000.0; // feet
    /// assert_eq_precision!(altitude, altitude_check, PRECISION);
    /// 
    /// let pressure = Pressure::new::<hectopascal>(1.767_99e3);
    /// let altitude = InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap().get::<uom::si::length::foot>();
    /// let altitude_check = -16_250.0; // feet
    /// assert_eq_precision!(altitude, altitude_check, PRECISION);
    /// ```
    pub fn altitude_from_pressure (pressure: uom::si::f64::Pressure) -> Result<uom::si::f64::Length, IsaError> {
        
        // input pressure in pascals
        let p_pa: f64 = pressure.get::<uom::si::pressure::pascal>();
        
        // sea-level standard pressure
        let p_0: f64 = 101325.0; // pascal
        
        // dbg!(p_pa, p_0;
        // set H_base and p_base for each ISA atmospheric level
        let h_b_and_p_b = vec![
            (-5.0, 1.75364e0 * p_0), 
            (0.0, p_0), 
            (11.0, 2.23361e-1*p_0), 
            (20.0, 5.40328e-2*p_0), 
            (32.0, 8.56664e-3*p_0), 
            (47.0, 1.09455e-3*p_0), 
            (51.0, 6.60631e-4*p_0), 
            (71.0, 3.90465e-5*p_0), 
            (80.0, 8.74682e-6*p_0), 

        ];

        // check if p is out of bounds.
        // If input pressure is less than the lowest or greater than the highest pressure in the table, return an error
        if p_pa > h_b_and_p_b[0].1 || p_pa < h_b_and_p_b[h_b_and_p_b.len()-1].1 { return Err(IsaError::InputOutOfRange) }

        let mut h_b: f64 = -10.0; // set a waaay low value as a place holder
        let mut p_b: f64 = -10.0; // set a waaay low value as a place holder

        // Figure out h_b and p_b?
        for layer_index in 0..h_b_and_p_b.len() - 1 {
            // if the input pressure is between the pressures for the lower and upper 
            // geopotential altitude bounds for the current layer...
            if p_pa <= h_b_and_p_b[layer_index].1 && p_pa > h_b_and_p_b[layer_index+1].1 {
                // set h_b and p_b for later use
                (h_b, p_b) = (h_b_and_p_b[layer_index].0, h_b_and_p_b[layer_index].1);
                // dbg_sci!(h_b, 5);
                // dbg_sci!(p_b, 5);
                break;
            }
        }
        
        // remaining constants/variables
        let r_star: f64 = Self::constant_universal_gas_constant(); // J/(kmol·K)
        let m: f64 = Self::constant_molar_mass_dry_air().get::<uom::si::molar_mass::gram_per_mole>(); // g/mol
        let g: f64 = Self::constant_gravity_sealevel().get::<uom::si::acceleration::meter_per_second_squared>(); // m/s^2
        let h_b_false: f64 = h_b + 1.0;
        // dbg!(h_b_false);
        let t_b: f64 = Self::altitude_to_base_temperature(uom::si::f64::Length::new::<uom::si::length::kilometer>(h_b_false))?.get::<uom::si::thermodynamic_temperature::kelvin>(); // K. Adds a small offset to avoid issues at layer boundaries
        let beta = Self::altitude_to_lapse_rate(uom::si::f64::Length::new::<uom::si::length::kilometer>(h_b_false)).unwrap(); // K/km. Adds a small offset to avoid issues at layer boundaries
        
        // dbg!(h_b, p_b, r_star, m, g, t_b, beta);

        if beta == 0.0 {
            // Calculate H for an Isothermal layer
            let h = h_b + (- (r_star * t_b * (p_pa / p_b).ln())) / (m * g);
            // dbg!(h);
            Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(h)) // return H as a uom length representing the geopotential altitude
        } else {
            // Calculate H for a Non-isothermal layer
            let h = h_b + (t_b / (beta * (p_pa/p_b).powf((r_star*beta)/(m*g)))) - (t_b / beta);
            // dbg!(h);
            Ok(uom::si::f64::Length::new::<uom::si::length::kilometer>(h)) // return H as a uom length representing the geopotential altitude
        }
    }

    /// Geopotential Altitude from density and temperature
    /// 
    /// This function returns the geopotential altitude that corresponds to the input density. Temperature is a 
    /// normal component in the various formulas to estimate air density. Although it may have been used in computing
    /// the density that is used as input to this formula, it is still required as an input to this function.
    /// 
    /// The formula for density in ICAO Doc 7488/3 is $\rho = \frac{p}{RT}$. $T$ is a required variable. This function 
    /// calculates the pressure from the input density and temperature, and then calls `altitude_from_pressure` 
    /// to get the corresponding geopotential altitude.
    /// 
    /// # Arguments
    /// * `rho` - Density (`uom::si::f64::MassDensity`)
    /// * `t` - Temperature (`uom::si::f64::ThermodynamicTemperature`)
    /// 
    /// # Returns
    /// Geopotential altitude (`uom::si::f64::Length`)
    /// 
    /// # Example
    /// ```
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use aero_atmos::{assert_eq_sigfigs, assert_eq_precision};
    /// 
    /// use uom::si::f64::{MassDensity, ThermodynamicTemperature, Length};
    /// use uom::si::mass_density::kilogram_per_cubic_meter;
    /// use uom::si::thermodynamic_temperature::kelvin;
    /// use uom::si::length::foot;
    ///
    /// let rho = MassDensity::new::<kilogram_per_cubic_meter>(7.45979e-1);
    /// let t = ThermodynamicTemperature::new::<kelvin>(256.451);
    /// let alt = InternationalStandardAtmosphere::altitude_from_density_and_temperature(rho, t).unwrap();
    /// assert_eq_sigfigs!(alt.get::<foot>(), 16_000.0, 5);
    /// ```
    pub fn altitude_from_density_and_temperature (rho: uom::si::f64::MassDensity, t: uom::si::f64::ThermodynamicTemperature) -> Result<uom::si::f64::Length, IsaError> {
        // get the values.
        let r = Self::constant_specific_gas_constant_air();
        let rho = rho.get::<uom::si::mass_density::kilogram_per_cubic_meter>();
        let t = t.get::<uom::si::thermodynamic_temperature::kelvin>();

        // calculate the pressure
        let p = uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(rho * r * t);

        // get height based on the pressure
        let alt = Self::altitude_from_pressure(p)?;

        Ok(alt)
    }

    /// The earth's radius.
    /// 
    /// This function returns the constant earth radius used in the ISA model.
    /// 
    /// /// The symbol is `R_e`.
    /// 
    /// # Returns
    /// Earth radius as uom length.
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
    /// The symbol is `G` or _g_0_.
    /// 
    /// "standard acceleration due to gravity. It conforms 
    /// with latitude n = 45°32'33"N using
    /// Lambert’s equation of the acceleration due
    /// to gravity as a function of latitude" (ICAO Doc 7488/3, section 2.2)
    pub fn constant_gravity_sealevel() -> uom::si::f64::Acceleration {
        uom::si::f64::Acceleration::new::<uom::si::acceleration::meter_per_second_squared>(9.806_65)
    }

    /// Molar mass of Earth's air.
    /// 
    /// Returns the molar mass of Earth's air used in the ISA model.
    /// 
    /// The symbol is `M_0`.
    pub fn constant_molar_mass_dry_air() -> uom::si::f64::MolarMass {
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
    /// 
    /// The symbol is `R`.
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
    /// The symbol is _β_s_.
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
    /// 
    /// # Example
    /// ```rust
    /// use aero_atmos::intl_standard_atmos::InternationalStandardAtmosphere;
    /// use uom::si::f64::Length;
    /// use uom::si::length::meter;
    /// 
    /// let collision_diameter = InternationalStandardAtmosphere::constant_effective_collision_diameter_air_molecule();
    /// assert_eq!(collision_diameter.get::<meter>(), 0.365e-9);
    /// ```
    pub fn constant_effective_collision_diameter_air_molecule() -> uom::si::f64::Length {
        uom::si::f64::Length::new::<uom::si::length::meter>(0.365e-9)
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

