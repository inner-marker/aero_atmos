/// Macro that checks whether two numerical values are approximately equal within a specified precision.
/// 
/// The specified precision is a ratio representing the maximum allowable relative difference between the two values.
/// 
/// Anything else after the first three arguments is not allowed.
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
    ($left:expr, $right:expr, $precision:expr $(,)?) => {{
        let left_val = $left;
        let right_val = $right;
        let precision_val = $precision;

        if left_val == right_val {
            // Values are exactly equal
        } else {
            let relative_difference = ((left_val - right_val) / right_val).abs();
            if relative_difference >= precision_val {
                panic!(
                    "assertion failed: `(left ≈ right)`\n  left: `{:?}`,\n right: `{:?}`,\n precision: `{:?}`,\n relative difference: `{:?}`",
                    left_val, right_val, precision_val, relative_difference
                );
            }
        }
    }};
}
