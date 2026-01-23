//! This macros module contains custom macros used primarily for testing.
//! 
//! The macros here are used primarily in testing to assert numerical equality within specified tolerances. However, they can be used in other contexts as needed.
//! 
//! # Examples
//! 
//! ## assert_eq_precision
//! 
//! This macro checks whether two numerical values are approximately equal within a specified precision (ratio). 
//! 
//! ```rust
//! use aero_atmos::assert_eq_precision;
//! let a: f64 = 100.0;
//! let b: f64 = 100.04;
//! let precision: f64 = 0.0005; // 0.05%
//! assert_eq_precision!(a, b, precision); // This will pass
//! 
//! let c: f64 = 100.1;
//! // assert_eq_precision!(a, c, precision); // This will panic
//! ```
//! 
//! ## assert_eq_sigfigs
//! 
//! This macro checks whether two numerical values are approximately equal when rounded to a specified number of significant figures. 
//! 
//! ```rust
//! use aero_atmos::assert_eq_sigfigs;
//! assert_eq_sigfigs!(1.11111, 1.11222, 3); // <- will pass     (1.11  == 1.11)
//! assert_eq_sigfigs!(1.1118, 1.1122, 4); // <- will pass    (1.112 == 1.112)
//! ```
//! 
//! ## dbg_sci
//! 
//! This macro prints a numerical value, the expression that generates it (variable), in scientific notation for debugging purposes.
//! Additionally, it prints the file name and line number where the macro is invoked. Notably, this macro is only effective in debug builds; it 
//! is not included in release builds.
//! 
//! ```rust
//! use aero_atmos::dbg_sci;
//! 
//! let x = 12345.6789;
//! dbg_sci!(x, 3); // prints something like "[src/main.rs:10] x = 1.234e4"
//! ```

/// Macro that checks whether two numerical values are approximately equal within a specified precision. Used in tests.
/// 
/// The specified precision is a ratio representing the maximum allowable relative difference between the two values.
/// 
/// # Arguments
/// - `left`: The first numerical value to compare.
/// - `right`: The second numerical value to compare.
/// - `precision`: The maximum allowable relative difference as a ratio (e.g., 0.01 for 1%).
/// 
/// # Examples
/// ```
/// use aero_atmos::assert_eq_precision;
/// let a: f64 = 100.0;
/// let b: f64 = 100.04;
/// let precision: f64 = 0.0005; // 0.05%
/// assert_eq_precision!(a, b, precision); // This will pass
/// 
/// let c: f64 = 100.1;
/// // assert_eq_precision!(a, c, precision); // This will panic
/// ```
#[macro_export]
macro_rules! assert_eq_precision {
    ($left:expr, $right:expr, $precision:expr) => {{
        let left_val = $left;
        let right_val = $right;
        let precision_val = $precision;

        if left_val == right_val {
            // Values are exactly equal
        } else {
            let relative_difference = ((left_val - right_val) / right_val).abs();
            // dbg!(relative_difference.clone());
            if relative_difference >= precision_val {
                panic!(
                    "Assertion failed, Values' relative difference is greater than the precision: `(left ≈ right)`\n  left: `{:?}`,\n right: `{:?}`,\n precision: `{:?}`,\n relative difference: `{:?}`",
                    left_val, right_val, precision_val, relative_difference
                );
            }
        }
    }};
}

/// Assert numbers rounded to a specified number of significant figures are equal. Used in tests.
/// 
/// This macro rounds both numbers to the specified number of significant figures
/// before comparing them for equality. The rounding will assume that the leftmost (greatest) 
/// digit is the first significant figure.
/// 
/// # Arguments
/// - `left`: The first numerical value to compare.
/// - `right`: The second numerical value to compare.
/// - `sigfigs`: The number of significant figures to which both values should be rounded.
/// 
/// # Example
/// ```rust
/// use rust_decimal::prelude::*;
/// use aero_atmos::assert_eq_sigfigs;
/// 
/// assert_eq_sigfigs!(1.11111, 1.11222, 3); // <- will pass     (1.11  == 1.11)
/// assert_eq_sigfigs!(1.1118, 1.1122, 4); // <- will pass    (1.112 == 1.112)
/// ```
#[macro_export]
macro_rules! assert_eq_sigfigs {
    ($left:expr, $right:expr, $sigfigs:expr) => {
        use rust_decimal::prelude::*;
        let left_rounded = rust_decimal::Decimal::from_f64($left).unwrap().round_sf_with_strategy($sigfigs, RoundingStrategy::MidpointAwayFromZero).unwrap();
        let right_rounded = rust_decimal::Decimal::from_f64($right).unwrap().round_sf_with_strategy($sigfigs, RoundingStrategy::MidpointAwayFromZero).unwrap();

        if left_rounded != right_rounded {
            panic!("Assertion Failed. Rounded values are not equal.
            left: {}
            right: {}
            left_rounded: {}
            right_rounded: {}",
            $left, $right, left_rounded, right_rounded
            )
        }
    };
}

/// Debug macro for numerical values in scientific notation
/// 
/// Similar to `dbg!()` but prints a numeric value in scientific notation.
/// 
/// Just as `dbg!()` prints the file name and line number along with the value, 
/// this macro prints the file name and line number along with the numeric value 
/// in scientific notation.
/// 
/// The functionality of this macro only works in debug builds. In release builds, the macro does nothing.
/// 
/// # Arguments
/// * `val` - The numeric value to print in scientific notation.
/// * `decimal_places` - The number of decimal places to include in the scientific notation.
/// 
/// # Examples
/// ```
/// use aero_atmos::dbg_sci;
/// 
/// let x = 12345.6789;
/// dbg_sci!(x, 3); // prints something like "[src/main.rs:10] x = 1.235e4"
/// ```
#[macro_export]
macro_rules! dbg_sci {
    ($val:expr, $decimal_places:expr) => {
        // prevent this from being included in release builds
        #[cfg(debug_assertions)]
        println!("[{}:{}] {} = {:.*e}", file!(), line!(), stringify!($val), $decimal_places, $val);
    };
}