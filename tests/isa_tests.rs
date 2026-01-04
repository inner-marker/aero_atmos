//! These are tests for the ISA atmospheric model

use aero_atmos::intl_standard_atmos::IsaError;
use uom::si::f64::{Length, ThermodynamicTemperature};
use uom::si::{length::*};
use uom::si::thermodynamic_temperature::{kelvin, degree_celsius};

// use aero_atmos::{InternationalStandardAtmosphere};
use aero_atmos::{assert_eq_precision, PRECISION};


/// Test for Standard temperature.
#[test]
fn standard_temperature_tests () {
    // sea level
    let length = Length::new::<kilometer>(0.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length);
    assert_eq!(temp, Ok(ThermodynamicTemperature::new::<kelvin>(288.15)), "Temperature is not equal.");

    // middle of the range
    let length = Length::new::<kilometer>(11.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length);
    assert_eq!(temp, Ok(ThermodynamicTemperature::new::<kelvin>(216.65)), "Temperature is not equal.");

    // test with feet as an input. Values from Doc7488, Table 4.
    let length = Length::new::<foot>(23_972.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert!(((temp - 240.656)/temp) < aero_atmos::PRECISION, "Computed temperature is not close enough to the Doc 7488 Table IV value.");
    // test with feet as an input. Values from Doc7488, Table 4.
    let length = Length::new::<foot>(-14_009.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().value;
    assert!(((temp - 315.905)/temp) < aero_atmos::PRECISION, "Computed temperature is not close enough to the Doc 7488 Table IV value.");
    // test with feet as an input, celcius as output. Values from Doc7488, Table 4.
    let length = Length::new::<foot>(-14_009.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length).unwrap().get::<degree_celsius>();
    assert!(((temp - 42.755)/temp) < aero_atmos::PRECISION, "Computed temperature is not close enough to the Doc 7488 Table IV value.");

    // out of bounds low
    let length = Length::new::<kilometer>(-6.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length);
    assert_eq!(temp, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");
    
    // out of bounds high
    let length = Length::new::<kilometer>(85.0);
    let temp = aero_atmos::InternationalStandardAtmosphere::altitude_to_temperature(length);
    assert_eq!(temp, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");
}

/// Test for Standard pressure.
#[test]
fn standard_pressure_tests () {
    // sea level
    let altitude = Length::new::<kilometer>(0.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(101325.0).value, PRECISION);

    // middle of the range, right on a layer boundary
    let altitude = Length::new::<kilometer>(11.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(22632.1).value, PRECISION);

    // withing range, but near the top of the range (based on Doc7488 Table 7)
    let altitude = Length::new::<kilometer>(31.985);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(870.0).value, PRECISION);

    // within range, negative altitude (based on Doc7488 Table 7)
    let altitude = Length::new::<kilometer>(-4.990);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude).unwrap().value;
    assert_eq_precision!(pressure, uom::si::f64::Pressure::new::<uom::si::pressure::pascal>(177500.0).value, PRECISION);

    // too low
    let altitude = Length::new::<kilometer>(-6.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude);
    assert_eq!(pressure, Err(IsaError::InputOutOfRange), "Geopotential altitude is too low out of bounds.");

    // too high
    let altitude = Length::new::<kilometer>(85.0);
    let pressure = aero_atmos::InternationalStandardAtmosphere::altitude_to_pressure(altitude);
    assert_eq!(pressure, Err(IsaError::InputOutOfRange), "Geopotential altitude is too high out of bounds.");

}