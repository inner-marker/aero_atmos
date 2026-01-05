/// Macro that checks whether two numerical values are approximately equal within a specified precision.
/// 
/// The specified precision is a ratio representing the maximum allowable relative difference between the two values.
/// 
/// # Examples
/// ```
/// use aero_atmos::assert_eq_precision;
/// let a = 100.0;
/// let b = 100.04;
/// let precision = 0.0005; // 0.05%
/// // assert_eq_precision!(a, b, precision); // This will pass
/// 
/// let c = 100.1;
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
            if relative_difference >= precision_val {
                panic!(
                    "Assertion failed, Values' relative difference is greater than the precision: `(left ≈ right)`\n  left: `{:?}`,\n right: `{:?}`,\n precision: `{:?}`,\n relative difference: `{:?}`",
                    left_val, right_val, precision_val, relative_difference
                );
            }
        }
    }};
}

/// Assert numbers rounded to a specified number of significant figures are equal.
/// 
/// This macro rounds both numbers to the specified number of significant figures
/// before comparing them for equality. The rounding will assume that the leftmost (greatest) 
/// digit is the first significant figure.
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