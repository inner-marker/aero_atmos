//! These are tests for the ISA atmospheric model

use aero_atmos::intl_standard_atmos::IsaError;
use uom::si::f64::{Length, ThermodynamicTemperature};
use uom::si::{Quantity, length::*};
use uom::si::thermodynamic_temperature::{kelvin, degree_celsius};

use aero_atmos::{InternationalStandardAtmosphere};


/// Test for Standard sea level temperature.
#[test]
fn standard_pressure () {
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

/// Tests for earth's radius
#[test]
fn earth_radius() {
    let radius = aero_atmos::InternationalStandardAtmosphere::constant_earth_radius();
    assert_eq!(radius, Length::new::<kilometer>(6_356.766));
}
