//! Main file for testing the library.
#![allow(unused)]

use uom::si::f64::{Pressure, Length};
use uom::si::length::*;
use uom::si::pressure::*;
use crate::macros::*;

use aero_atmos::*;


fn main () {
    // println!("\nTesting pressure to altitude...");
    // let pressure = Pressure::new::<hectopascal>(1.76799e3);
    // let geop_alt = InternationalStandardAtmosphere::altitude_from_pressure(pressure);
    // dbg!(geop_alt.unwrap().get::<foot>());

    // println!("\nTesting altitude to pressure...");
    // let geop_alt = Length::new::<foot>(-16_250.0);
    // let pressure = InternationalStandardAtmosphere::altitude_to_pressure(geop_alt);
    // dbg_sci!(pressure.clone().unwrap().get::<hectopascal>(), 5);

    println!("\nTesting altitude to pressure...");
    let geop_alt = Length::new::<foot>(0.0);
    let pressure = InternationalStandardAtmosphere::altitude_to_pressure(geop_alt);
    dbg_sci!(pressure.clone().unwrap().get::<hectopascal>(), 5);
}