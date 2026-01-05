/// Macro that wraps around the `assert_eq!` macro to check whether two numerical values are approximately equal within a specified precision.
/// 
/// Anything else after the first three arguments is ignored.
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

