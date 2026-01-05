//! These are tests for the ISA atmospheric model

use aero_atmos::intl_standard_atmos::IsaError;
use uom::si::{length::*};

// use aero_atmos::{InternationalStandardAtmosphere};
use aero_atmos::{assert_eq_precision, PRECISION};


/// Test for Standard temperature.
#[test]
fn standard_temperature_tests () {
    // sea level
    let length = uom::si::f64::Length::new::<kilometer>(0.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert_eq_precision!(temp, 288.15, PRECISION);

    // middle of the range
    let length = uom::si::f64::Length::new::<kilometer>(11.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert_eq_precision!(temp, 216.65, PRECISION);

    // test with feet as an input. Values from Doc7488, Table 4.
    let length = uom::si::f64::Length::new::<foot>(23_972.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert_eq_precision!(temp, 240.656, PRECISION);

    // test with feet as an input. Values from Doc7488, Table 4.
    let length = uom::si::f64::Length::new::<foot>(-14_009.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert_eq_precision!(temp, 315.905, PRECISION);

    // test with feet as an input, celcius as output. Values from Doc7488, Table 4.
    let length = uom::si::f64::Length::new::<foot>(-14_009.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().get::<uom::si::thermodynamic_temperature::degree_celsius>();
    assert_eq_precision!(temp, 42.755, PRECISION);

    // very high
    // H = 262_000 ft => 196.935 K
    let length = uom::si::f64::Length::new::<foot>(262_000.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert_eq_precision!(temp, 196.935, PRECISION);

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
    // sea level
    let altitude = uom::si::f64::Length::new::<kilometer>(0.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(101325.0).value, PRECISION);

    // middle of the range, right on a layer boundary
    let altitude = uom::si::f64::Length::new::<kilometer>(11.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(22632.1).value, PRECISION);

    // withing range, but near the top of the range (based on Doc7488 Table 7)
    let altitude = uom::si::f64::Length::new::<kilometer>(31.985);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(870.0).value, PRECISION);

    // within range, negative altitude (based on Doc7488 Table 7)
    let altitude = uom::si::f64::Length::new::<kilometer>(-4.990);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(177500.0).value, PRECISION);

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
    // sea level
    let altitude = uom::si::f64::Length::new::<kilometer>(0.0);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    assert_eq_precision!(pressure_ratio, 1.0, PRECISION);

    // mid range, not on a layer boundary
    // H = 20 km => 5.40328e–2
    let altitude = uom::si::f64::Length::new::<kilometer>(20.0);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    assert_eq_precision!(pressure_ratio, 5.40328e-2, PRECISION);

    // Low end of the range
    // Doc7488 Table 2 value at H = -4.950 km => 1.74431
    let altitude = uom::si::f64::Length::new::<kilometer>(-4.950);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    assert_eq_precision!(pressure_ratio, 1.74431, PRECISION);

    // High end of the range
    // Doc7488 Table 2 value at 79.800 km => 9.05575e-6
    let altitude = uom::si::f64::Length::new::<kilometer>(79.800);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    assert_eq_precision!(pressure_ratio, 9.05575e-6, PRECISION);

    // Bottom of the range
    // H = -5 km => 1.75364
    let altitude = uom::si::f64::Length::new::<kilometer>(-5.0);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    assert_eq_precision!(pressure_ratio, 1.75364, PRECISION);

    // Top of the range
    // H = 80.0 km => 8.74682e-6
    let altitude = uom::si::f64::Length::new::<kilometer>(80.0);
    let pressure_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure_ratio(altitude).unwrap();
    assert_eq_precision!(pressure_ratio, 8.74682e-6, PRECISION);
}

/// Density tests
#[test]
fn standard_density_tests () {
    // sea level
    let altitude = uom::si::f64::Length::new::<kilometer>(0.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().value;
    assert_eq_precision!(density, uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(1.225).value, PRECISION);

    // mid range, not on a layer boundary, in feet (Doc 7488 Table 4 gives H in feet)
    // H = 65,600 ft => 8.81057e-2 kg/m^3
    let altitude = uom::si::f64::Length::new::<foot>(65_600.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().value;
    assert_eq_precision!(density, uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(8.81057e-2).value, PRECISION);

    // mid range on a layer boundary (Doc 7488 Table 4 gives H in feet)
    // H = 11 km = 36_089.24 ft => 3.63918e-1 kg/m^3
    let altitude = uom::si::f64::Length::new::<foot>(36_089.24);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().value;
    assert_eq_precision!(density, uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(3.63918e-1).value, PRECISION);

    // midrange
    // H = 65_000 ft => 9.06835e-2 kg/m^3
    let altitude = uom::si::f64::Length::new::<foot>(65_000.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().value;
    assert_eq_precision!(density, uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(9.06835e-2).value, PRECISION);

    // negative altitude (Doc 7488 Table 4 gives H in feet)
    // H = –16_250 ft => 1.92265 kg/m^3
    let altitude = uom::si::f64::Length::new::<foot>(-16_250.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().value;
    assert_eq_precision!(density, uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(1.92265).value, PRECISION);

    // near the top of the range (Doc 7488 Table 7)
    // H = 262_000 ft => 1.60701e-05 kg/m^3
    let altitude = uom::si::f64::Length::new::<foot>(262_000.0);
    let density = aero_atmos::InternationalStandardAtmosphere::altitude_to_density(altitude).unwrap().value;
    assert_eq_precision!(density, uom::si::f64::MassDensity::new::<uom::si::mass_density::kilogram_per_cubic_meter>(1.60701e-5).value, PRECISION);

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
    // sea level
    let altitude = uom::si::f64::Length::new::<kilometer>(0.0);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
    assert_eq_precision!(density_ratio, 1.0, PRECISION);

    // mid range, not on a layer boundary
    // H = 20 km => 7.18649e-2
    let altitude = uom::si::f64::Length::new::<kilometer>(20.0);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
    assert_eq_precision!(density_ratio, 7.18649e-2, PRECISION);

    // mid range, not on a layer boundary
    // H = 10 km => 3.36903e-1
    let altitude = uom::si::f64::Length::new::<kilometer>(10.0);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
    assert_eq_precision!(density_ratio, 3.36903e-1, PRECISION);

    // Low end of the range
    // Doc7488 Table 2 value at H = -4.950 km => 1.56911
    let altitude = uom::si::f64::Length::new::<kilometer>(-4.950);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
    assert_eq_precision!(density_ratio, 1.56911, PRECISION);

    // High end of the range
    // Doc7488 Table 2 value at 79.800 km => 1.32424e-5
    let altitude = uom::si::f64::Length::new::<kilometer>(79.800);
    let density_ratio = aero_atmos::InternationalStandardAtmosphere::altitude_to_density_ratio(altitude).unwrap();
    assert_eq_precision!(density_ratio, 1.32424e-5, PRECISION);

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

    // Test values from Doc7488, Table 4, given as ('H' in ft, 'a' in m/s)
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
    // Test values from Doc7488, Table 4, given as ('H' in ft, 'mu' in microPascal seconds)
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