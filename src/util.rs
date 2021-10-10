///
/// A type that has a zero value.
///
/// Note that if the type also conforms to the [Add](std::ops::Add)
/// or [AddAssign](std::ops::AddAssign) trait the zero element should
/// act neutral (for all element: element + zero() == element)
///
pub trait Zeroed: Copy + PartialEq {
    /// The constant zero value for the type.
    fn zero() -> Self;

    /// A indicater to check a zero value.
    #[inline(always)]
    fn is_zero(&self) -> bool {
        self == &Self::zero()
    }
}

macro_rules! apply(
    ($name:ty, $zero:expr) => (
        impl Zeroed for $name {
            #[inline(always)]
            fn zero() -> Self {
                $zero
            }
        }
    );
    ($name:ty) => (
        apply!($name, 0);
    );
);

apply!(bool, false);

apply!(u8);
apply!(u16);
apply!(u32);
apply!(u64);
apply!(usize);

apply!(i8);
apply!(i16);
apply!(i32);
apply!(i64);
apply!(isize);

apply!(f32, 0.0);
apply!(f64, 0.0);

use num_complex::Complex32;
use num_complex::Complex64;

apply!(Complex32, Complex32::new(0.0, 0.0));
apply!(Complex64, Complex64::new(0.0, 0.0));
