//! These are tests for the ISA atmospheric model

use aero_atmos::intl_standard_atmos::IsaError;
use uom::si::{length::*, pressure::hectopascal};

// use aero_atmos::{InternationalStandardAtmosphere};
use aero_atmos::{assert_eq_precision, PRECISION};


/// Test for Standard temperature.
#[test]
fn standard_temperature_tests () {
    // test data (geopotential altitude in feet, temperature in kelvin)
    // Table 4 from Doc7488
    let test_data = [
        (-16_250.0,  320.325),    // very low in the range
        (      0.0,  288.150),    // sea level
        ( 16_000.0,  256.451),    // mid range
        ( 50_000.0,  216.650),    // mid range
        (100_000.0,  227.130),    // mid range
        (262_000.0,  196.935),    // very high in the range
    ];

    for (h_ft, temp_expected) in test_data {
        let length = uom::si::f64::Length::new::<uom::si::length::foot>(h_ft);
        let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
        assert_eq_precision!(temp, temp_expected, PRECISION);
    }

    // out of bounds low
    let length = uom::si::f64::Length::new::<kilometer>(-6.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length);
    assert_eq!(temp, Err(IsaError::InputOutOfRange));
    
    // out of bounds high
    let length = uom::si::f64::Length::new::<kilometer>(85.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length);
    assert_eq!(temp, Err(IsaError::InputOutOfRange));
}

/// Test for Standard pressure.
#[test]
fn standard_pressure_tests () {
    // test data base on Table 4: 
    // vec of (geopotential altitude in feet, pressure in hectopascals)
    let test_data = [
        (-16_250.0, 1.76799e3 ),    // very low in the range
        (      0.0, 1.01325e3 ),   // sea level
        ( 16_000.0, 5.49152e2 ),    // mid range
        ( 50_000.0, 1.15972e2 ),    // mid range
        (100_000.0, 1.09015e1 ),   // mid range
        (262_000.0, 9.08455e-3), // very high in the range
    ];
    
    for (h_ft, pressure_hpa_expected) in test_data {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().get::<hectopascal>();
        assert_eq_precision!(pressure, pressure_hpa_expected, PRECISION);
    }
    
    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude);
    assert_eq!(pressure, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");
    
    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude);
    assert_eq!(pressure, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
    
}

/// Test pressure ratios from altitude.
#[test]
fn standard_pressure_ratio_tests () {
    // test data based on Table 5:
    // (geopotential altitude in feet, pressure ratio)
    let test_data = [
        (-16_000.0, 1.73074e0),    // very low in the range
        (      0.0, 1.0),   // sea level
        ( 16_000.0, 5.41971e-1),   // mid range
        ( 20_000.0, 4.59543e-1),   // mid range
        ( 50_000.0, 1.14456e-1),   // mid range
        (100_000.0, 1.07590e-2),   // mid range
        (262_000.0, 8.96576e-6),   // very high in the range
    ];

    for (h_ft, pressure_ratio_expected) in test_data {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
        assert_eq_precision!(pressure_ratio, pressure_ratio_expected, PRECISION);
    }
    
    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude);
    assert_eq!(pressure_ratio, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude);
    assert_eq!(pressure_ratio, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Density tests
#[test]
fn standard_density_tests () {
    // test data base on Table 4: 
    // vec of (geopotential altitude in feet, density in kg/m^3)
    let test_data = [
        (-16_250.0, 1.92265e0 ),    // very low in the range
        (      0.0, 1.22500e0 ),   // sea level
        ( 16_000.0, 7.45979e-1 ),    // mid range
        ( 50_000.0, 1.86481e-1 ),    // mid range
        (100_000.0, 1.67206e-2 ),   // mid range
        (262_000.0, 1.60701e-5 ), // very high in the range
    ];

    for (h_ft, density_expected) in test_data {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().get::<uom::si::mass_density::kilogram_per_cubic_meter>();
        assert_eq_precision!(density, density_expected, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude);
    assert_eq!(density, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude);
    assert_eq!(density, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Density ratio tests
#[test]
fn standard_density_ratio_tests () {
    // (alt in km, density ratio)
    // data from Doc7488, Table 2
    let test_data = vec! [
        (-4.95,      1.56911),
        (0.0,        1.0),
        (10.0,       3.36903e-1),
        (20.0,       7.18649e-2),
        (79.8,       1.32424e-5),
    ];

    for (h_km, density_ratio_expected) in test_data {
        let altitude = uom::si::f64::Length::new::<kilometer>(h_km);
        let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
        assert_eq_precision!(density_ratio, density_ratio_expected, PRECISION);
    }

    // Too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude);
    assert_eq!(density_ratio, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // Too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude);
    assert_eq!(density_ratio, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Speed of sound tests
#[test]
fn standard_speed_of_sound_tests () {

    // Test values from Doc7488, Table 2, given as ('H' in ft, 'a' in m/s)
    let test_values = vec![
        (-16_000.0, 358.523),      // very low in the range
        (0.0,       340.294),      // sea level
        (16_000.0,  321.031),
        (20_000.0,  316.032),
        (100_000.0, 302.122),
        (262_000.0, 281.324),      // very high in the range
    ];

    for (h_ft, a_m_per_s) in test_values {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let speed_of_sound = aero_atmos::InternationalStandardAtmosphere::altitude_to_speed_of_sound(altitude).unwrap().value;
        assert_eq_precision!(speed_of_sound, a_m_per_s, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let speed_of_sound = aero_atmos::InternationalStandardAtmosphere::altitude_to_speed_of_sound(altitude);
    assert_eq!(speed_of_sound, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let speed_of_sound = aero_atmos::InternationalStandardAtmosphere::altitude_to_speed_of_sound(altitude);
    assert_eq!(speed_of_sound, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
    
}

/// Dynamic viscosity tests
#[test]
fn standard_dynamic_viscosity_tests () {
    // Test values from Doc7488, Table 5, given as ('H' in ft, 'mu' in microPascal seconds)
    let test_values = vec![
        (-16_000.0,  1.9385e-5),     // very low in the range
        (0.0,        1.7894e-5),     // sea level
        (16_000.0,   1.6322e-5),
        (20_000.0,   1.5915e-5),
        (100_000.0,  1.4786e-5),
        (262_000.0,  1.3111e-5),     // very high in the range
    ];

    for (h_ft, mu_pascal_second) in test_values {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let dynamic_viscosity = aero_atmos::InternationalStandardAtmosphere::altitude_to_dynamic_viscosity(altitude).unwrap().get::<uom::si::dynamic_viscosity::pascal_second>();
        assert_eq_precision!(dynamic_viscosity, mu_pascal_second, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let dynamic_viscosity = aero_atmos::InternationalStandardAtmosphere::altitude_to_dynamic_viscosity(altitude);
    assert_eq!(dynamic_viscosity, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let dynamic_viscosity = aero_atmos::InternationalStandardAtmosphere::altitude_to_dynamic_viscosity(altitude);
    assert_eq!(dynamic_viscosity, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Thermal conductivity tests
#[test]
fn standard_thermal_conductivity_tests () {
    // Test values from Doc7488, Table 5, given as ('H' in ft, 'lambda' in W/(m·K))
    let test_values = vec![
        (-16_000.0,  2.7798e-2),     // very low in the range
        (0.0,        2.5343e-2),     // sea level
        (16_000.0,   2.2810e-2),     // mid range
        (20_000.0,   2.2164e-2),     // mid range
        (100_000.0,  2.0397e-2),     // mid range
        (262_000.0,  1.7841e-2),     // very high in the range
    ];

    for (h_ft, lambda_w_m_k) in test_values {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let thermal_conductivity = aero_atmos::InternationalStandardAtmosphere::altitude_to_thermal_conductivity(altitude).unwrap().get::<uom::si::thermal_conductivity::watt_per_meter_kelvin>();
        assert_eq_precision!(thermal_conductivity, lambda_w_m_k, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let thermal_conductivity = aero_atmos::InternationalStandardAtmosphere::altitude_to_thermal_conductivity(altitude);
    assert_eq!(thermal_conductivity, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let thermal_conductivity = aero_atmos::InternationalStandardAtmosphere::altitude_to_thermal_conductivity(altitude);
    assert_eq!(thermal_conductivity, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Number Density tests
#[test]
fn standard_number_density_tests () {
    // test data (geopotential altitude in kilometers, number density in molecules per cubic meter)
    // Table 3 from Doc7488
    let test_values = vec![
        (-4.95,  3.9967e25),     // very low in the range
        (0.0  ,  2.5471e25),     // sea level
        (3.0  ,  1.8903e25),     // mid range
        (15.0 ,  4.0270e24),     // mid range
        (29.0 ,  4.3752e23),     // mid range
        (40.0 ,  8.0074e22),     // mid range
        (79.8 ,  3.3730e20),     // very high in the range
    ];

    for (h_km, n_per_m3) in test_values {
        let altitude = uom::si::f64::Length::new::<kilometer>(h_km);
        let number_density: f64 = aero_atmos::InternationalStandardAtmosphere::altitude_to_number_density(altitude).unwrap();
        assert_eq_precision!(number_density, n_per_m3, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let number_density = aero_atmos::InternationalStandardAtmosphere::altitude_to_number_density(altitude);
    assert_eq!(number_density, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let number_density = aero_atmos::InternationalStandardAtmosphere::altitude_to_number_density(altitude);
    assert_eq!(number_density, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}