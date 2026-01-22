//! These are tests for the ISA atmospheric model

use aero_atmos::{assert_eq_sigfigs, intl_standard_atmos::IsaError};
use uom::si::{length::*, pressure::hectopascal};

// use aero_atmos::{InternationalStandardAtmosphere};
use aero_atmos::{assert_eq_precision, PRECISION};


/// Test for Standard temperature.
#[test]
fn test_standard_temperature_tests () {
    // test data (geopotential altitude in feet, temperature in kelvin)
    // Table 4 from Doc7488
    let test_data = [
        (-16_250.0,  320.345),    // very low in the range
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
        assert_eq_sigfigs!(temp, temp_expected, 6);
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
fn test_altitude_to_pressure () {
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
fn test_standard_pressure_ratio () {
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
fn test_standard_density () {
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
fn test_standard_density_ratio () {
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
fn test_standard_speed_of_sound () {

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
fn test_standard_dynamic_viscosity () {
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
fn test_standard_thermal_conductivity () {
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
fn test_standard_number_density () {
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

/// Mean particle speed tests
#[test]
fn test_standard_mean_particle_speed () {
    // Test values from Doc7488, ('H' in meters, 'v_bar' in m/s)
    // Table 3 from Doc7488
    let test_values = vec![
        (-4.95,  483.89),     // very low in the range
        (0.0  ,  458.94),     // sea level
        (3.0  ,  443.14),     // mid range
        (15.0 ,  397.95),     // mid range
        (29.0 ,  406.13),     // mid range
        (40.0 ,  428.38),     // mid range
        (79.8 ,  379.52),     // very high in the range
    ];

    for (h_km, v_bar_m_per_s) in test_values {
        let altitude = uom::si::f64::Length::new::<kilometer>(h_km);
        let mean_particle_speed = aero_atmos::InternationalStandardAtmosphere::altitude_to_mean_particle_speed(altitude).unwrap().value;
        assert_eq_precision!(mean_particle_speed, v_bar_m_per_s, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let mean_particle_speed = aero_atmos::InternationalStandardAtmosphere::altitude_to_mean_particle_speed(altitude);
    assert_eq!(mean_particle_speed, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let mean_particle_speed = aero_atmos::InternationalStandardAtmosphere::altitude_to_mean_particle_speed(altitude);
    assert_eq!(mean_particle_speed, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Mean free path tests
#[test]
fn test_standard_mean_free_path () {
    // Test values from Doc7488, ('H' in meters, 'l' in meters)
    // Table 3 from Doc7488
    let test_values = vec![
        (-4.95,  4.2271e-8),     // very low in the range
        (0.0  ,  6.6328e-8),     // sea level
        (3.0  ,  8.9374e-8),     // mid range
        (15.0 ,  4.1953e-7),     // mid range
        (29.0 ,  3.8614e-6),     // mid range
        (40.0 ,  2.1099e-5),     // mid range
        (79.8 ,  5.0088e-3),     // very high in the range
    ];

    for (h_km, l_m) in test_values {
        let altitude = uom::si::f64::Length::new::<kilometer>(h_km);
        let mean_free_path = aero_atmos::InternationalStandardAtmosphere::altitude_to_mean_free_path(altitude).unwrap().get::<uom::si::length::meter>();
        assert_eq_precision!(mean_free_path, l_m, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let mean_free_path = aero_atmos::InternationalStandardAtmosphere::altitude_to_mean_free_path(altitude);
    assert_eq!(mean_free_path, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");
    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let mean_free_path = aero_atmos::InternationalStandardAtmosphere::altitude_to_mean_free_path(altitude);
    assert_eq!(mean_free_path, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Collision frequency tests
#[test]
fn test_standard_collision_frequency () {
    // Test values from Doc7488, ('H' in meters, 'omega' in Hz)
    // Table 3 from Doc7488
    let test_values = vec![
        (-4.95,  1.1447e10),     // very low in the range
        (0.0  ,  6.9193e9),     // sea level
        (3.0  ,  4.9583e9),     // mid range
        (15.0 ,  9.4857e8),     // mid range
        (29.0 ,  1.0518e8),     // mid range
        (40.0 ,  2.0303e7),     // mid range
        (79.8 ,  7.5772e4),     // very high in the range
    ];
    
    for (h_km, omega_hz) in test_values {
        let altitude = uom::si::f64::Length::new::<kilometer>(h_km);
        let collision_frequency = aero_atmos::InternationalStandardAtmosphere::altitude_to_collision_frequency(altitude).unwrap().value;
        assert_eq_precision!(collision_frequency, omega_hz, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let collision_frequency = aero_atmos::InternationalStandardAtmosphere::altitude_to_collision_frequency(altitude);
    assert_eq!(collision_frequency, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let collision_frequency = aero_atmos::InternationalStandardAtmosphere::altitude_to_collision_frequency(altitude);
    assert_eq!(collision_frequency, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Test for gravity at altitude
/// Tabular values for gravitational acceleration are given in ICAO Doc 7488/3 Table 4.
/// Geopotential altitude in Table 4 is given in feet.
#[test]
fn test_standard_gravity () {
    // test data (geopotential altitude in feet, gravity in m/s^2)
    let test_data = [
        (-16_250.0, 9.8219),    // very low in the range
        (      0.0, 9.8067),   // sea level
        ( 16_000.0, 9.7916),    // mid range
        ( 50_000.0, 9.7597),    // mid range
        (100_000.0, 9.7128),    // mid range
        (262_000.0, 9.5618),    // very high in the range
    ];

    for (h_ft, g_expected) in test_data {
        let altitude = uom::si::f64::Length::new::<foot>(h_ft);
        let gravity = aero_atmos::InternationalStandardAtmosphere::altitude_to_gravitational_acceleration(altitude).unwrap().get::<uom::si::acceleration::meter_per_second_squared>();
        assert_eq_precision!(gravity, g_expected, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let gravity = aero_atmos::InternationalStandardAtmosphere::altitude_to_gravitational_acceleration(altitude);
    assert_eq!(gravity, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let gravity = aero_atmos::InternationalStandardAtmosphere::altitude_to_gravitational_acceleration(altitude);
    assert_eq!(gravity, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Specific weight tests
/// Tabular values for specific weight are given in ICAO Doc 7488/3 Table 3.
/// Values for geopotential altitude are given in meters.
#[test]
fn test_standard_specific_weight_tests () {
    // test data base on Table 4: 
    // vec of (geopotential altitude in kilometers, specific weight in N/m^3)
    let test_data = [
        (-4.95,  1.8879e1),     // very low in the range
        (0.0  ,  1.2013e1),     // sea level
        (3.0  ,  8.9070e0),     // mid range
        (15.0 ,  1.8903e0),     // mid range
        (29.0 ,  2.0447e-1),     // mid range
        (40.0 ,  3.7291e-2),     // mid range
        (79.8 ,  1.5511e-4),     // very high in the range
    ];

    for (h_km, specific_weight_expected) in test_data {
        let altitude = uom::si::f64::Length::new::<kilometer>(h_km);
        let specific_weight = aero_atmos::InternationalStandardAtmosphere::altitude_to_specific_weight(altitude).unwrap();
        assert_eq_precision!(specific_weight, specific_weight_expected, PRECISION);
    }

    // too low
    let altitude = uom::si::f64::Length::new::<kilometer>(-6.0);
    let specific_weight = aero_atmos::InternationalStandardAtmosphere::altitude_to_specific_weight(altitude);
    assert_eq!(specific_weight, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = uom::si::f64::Length::new::<kilometer>(85.0);
    let specific_weight = aero_atmos::InternationalStandardAtmosphere::altitude_to_specific_weight(altitude);
    assert_eq!(specific_weight, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Test altitude_from_pressure
/// Data from Table 4, H in feet
#[test]
fn test_altitude_from_pressure () {
    // (pressure in pascal, expected altitude in feet)
    let test_data = [
        (1.01325e3,        0.0),       // sea level
        (5.49152e2,    16_000.0),     
        (1.15972e2,    50_000.0),     
        (1.09015e1,   100_000.0),     
        (9.08455e-3,  262_000.0),     
    ];

    for (p_pa, h_expected_ft) in test_data {
        let pressure = uom::si::f64::Pressure::new::<uom::si::pressure::hectopascal>(p_pa);
        let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap();
        assert_eq_precision!(altitude.get::<uom::si::length::foot>(), h_expected_ft, PRECISION);
    }

    // test data for less than sea level
    // (pressure (hpa), altitude (feet), precision (f64))
    let test_data_below_sea_level = [
        (1.76799e3,  -16_250.0, PRECISION),
        (1.64246e3,  -14_000.0, PRECISION),
        (1.53702e3,  -12_000.0, PRECISION),
        (1.43714e3,  -10_000.0, PRECISION),
        (1.34258e3,  -8_000.0,  PRECISION),
        (1.25312e3,  -6_000.0,  PRECISION),
        (1.16855e3,  -4_000.0,  PRECISION),
        (1.08866e3,  -2_000.0,  PRECISION),
    ];

    for (p_pa, h_expected_ft, precision) in test_data_below_sea_level {
        let pressure = uom::si::f64::Pressure::new::<uom::si::pressure::hectopascal>(p_pa);
        let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_pressure(pressure).unwrap();
        assert_eq_precision!(altitude.get::<uom::si::length::foot>(), h_expected_ft, precision);
    }

    // too low
    let pressure = uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(0.0);
    let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_pressure(pressure);
    assert_eq!(altitude, Err(IsaError::InputOutOfRange), "Pressure is too low out of bounds.");

    // too high
    let pressure = uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(200_000.0);
    let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_pressure(pressure);
    assert_eq!(altitude, Err(IsaError::InputOutOfRange), "Pressure is too high out of bounds.");
}

/// Test altitude from density and temperature
#[test]
fn test_altitude_from_density_and_temperature () {
    // (density in kg/m^3, temperature in K, expected altitude in feet)
    let test_data = [
        (1.92265e0,   320.345,    -16_250.0),
        //(1.22500e0,   288.150,          0.0),
        (1.21785e0,   287.754,        200.0),
        (7.45979e-1,  256.451,     16_000.0),
        (1.86481e-1,  216.650,     50_000.0),
        (1.60701e-5,  196.935,    262_000.0),

    ];

    for (rho_kg_m3, t_k, h_expected_ft) in test_data {
        let density = uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(rho_kg_m3);
        let temperature = uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(t_k);
        let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_density_and_temperature(density, temperature).unwrap();
        assert_eq_precision!(altitude.get::<uom::si::length::foot>(), h_expected_ft, PRECISION);
    }

    // too low
    let density = uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(0.0);
    let temperature = uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(288.15);
    let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_density_and_temperature(density, temperature);
    assert_eq!(altitude, Err(IsaError::InputOutOfRange), "Density is too low out of bounds.");

    // too high
    let density = uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(10.0);
    let temperature = uom::si::f64::ThermodynamicTemperature::new::<uom::si::thermodynamic_temperature::kelvin>(288.15);
    let altitude = aero_atmos::InternationalStandardAtmosphere::altitude_from_density_and_temperature(density, temperature);
    assert_eq!(altitude, Err(IsaError::InputOutOfRange), "Density is too high out of bounds.");
}