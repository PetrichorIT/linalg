use std::{
    fmt::{Binary, Debug, Display},
    num::ParseIntError,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign},
};

use crate::num::Abs;

///
/// Interger with arithemtic module N aka. GF(N).
///
#[derive(Clone, Copy, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GF<const N: usize>(pub usize);

impl<const N: usize> Abs for GF<N> {
    fn abs(&self) -> Self {
        *self
    }
}

// == ADD ==

impl<const N: usize> Add<GF<N>> for GF<N> {
    type Output = GF<N>;
    fn add(self, rhs: GF<N>) -> Self::Output {
        GF((self.0 + rhs.0) % N)
    }
}

impl<'a, const N: usize> Add<GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn add(self, rhs: GF<N>) -> Self::Output {
        GF((self.0 + rhs.0) % N)
    }
}

impl<const N: usize> Add<&'_ GF<N>> for GF<N> {
    type Output = GF<N>;
    fn add(self, rhs: &'_ GF<N>) -> Self::Output {
        GF((self.0 + rhs.0) % N)
    }
}

impl<'a, const N: usize> Add<&'_ GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn add(self, rhs: &'_ GF<N>) -> Self::Output {
        GF((self.0 + rhs.0) % N)
    }
}

impl<const N: usize> AddAssign<GF<N>> for GF<N> {
    fn add_assign(&mut self, rhs: GF<N>) {
        *self = GF((self.0 + rhs.0) % N)
    }
}

impl<const N: usize> AddAssign<&'_ GF<N>> for GF<N> {
    fn add_assign(&mut self, rhs: &'_ GF<N>) {
        *self = GF((self.0 + rhs.0) % N)
    }
}

// FMT

impl<const N: usize> Debug for GF<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(&self.0, f)
    }
}

impl<const N: usize> Binary for GF<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Binary::fmt(&self.0, f)
    }
}

impl<const N: usize> Display for GF<N> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

// == DIV ==

impl<const N: usize> Div<GF<N>> for GF<N> {
    type Output = GF<N>;
    fn div(self, rhs: GF<N>) -> Self::Output {
        Div::div(&self, &rhs)
    }
}

impl<'a, const N: usize> Div<GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn div(self, rhs: GF<N>) -> Self::Output {
        Div::div(self, &rhs)
    }
}

impl<const N: usize> Div<&'_ GF<N>> for GF<N> {
    type Output = GF<N>;
    fn div(self, rhs: &'_ GF<N>) -> Self::Output {
        Div::div(&self, rhs)
    }
}

impl<'a, const N: usize> Div<&'_ GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn div(self, rhs: &'_ GF<N>) -> Self::Output {
        // Call: a / b = c
        // this means c * b = a
        // so that (a / b) * b = a

        for i in 0..N {
            if GF(i) * rhs == *self {
                return GF(i);
            }
        }

        panic!("Found no mutltiplicative inverse. Note that N must be a primepower")
    }
}

impl<const N: usize> DivAssign<GF<N>> for GF<N> {
    fn div_assign(&mut self, rhs: GF<N>) {
        *self = *self * rhs;
    }
}

impl<const N: usize> DivAssign<&'_ GF<N>> for GF<N> {
    fn div_assign(&mut self, rhs: &'_ GF<N>) {
        *self = *self * rhs;
    }
}

// == FROM ==

impl From<bool> for GF<2> {
    fn from(b: bool) -> Self {
        GF(if b { 1 } else { 0 })
    }
}

impl<const N: usize> From<usize> for GF<N> {
    fn from(v: usize) -> Self {
        GF(v % N)
    }
}

// == MUL ==

impl<const N: usize> Mul<GF<N>> for GF<N> {
    type Output = GF<N>;
    fn mul(self, rhs: GF<N>) -> Self::Output {
        GF((self.0 * rhs.0) % N)
    }
}

impl<'a, const N: usize> Mul<GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn mul(self, rhs: GF<N>) -> Self::Output {
        GF((self.0 * rhs.0) % N)
    }
}

impl<const N: usize> Mul<&'_ GF<N>> for GF<N> {
    type Output = GF<N>;
    fn mul(self, rhs: &'_ GF<N>) -> Self::Output {
        GF((self.0 * rhs.0) % N)
    }
}

impl<'a, const N: usize> Mul<&'_ GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn mul(self, rhs: &'_ GF<N>) -> Self::Output {
        GF((self.0 * rhs.0) % N)
    }
}

impl<const N: usize> MulAssign<GF<N>> for GF<N> {
    fn mul_assign(&mut self, rhs: GF<N>) {
        *self = GF((self.0 * rhs.0) % N);
    }
}

impl<const N: usize> MulAssign<&'_ GF<N>> for GF<N> {
    fn mul_assign(&mut self, rhs: &'_ GF<N>) {
        *self = GF((self.0 * rhs.0) % N);
    }
}

// SUB

impl<const N: usize> Sub<GF<N>> for GF<N> {
    type Output = GF<N>;
    fn sub(self, rhs: GF<N>) -> Self::Output {
        GF(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            N - (rhs.0 - self.0)
        })
    }
}

impl<'a, const N: usize> Sub<GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn sub(self, rhs: GF<N>) -> Self::Output {
        GF(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            N - (rhs.0 - self.0)
        })
    }
}

impl<const N: usize> Sub<&'_ GF<N>> for GF<N> {
    type Output = GF<N>;
    fn sub(self, rhs: &'_ GF<N>) -> Self::Output {
        GF(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            N - (rhs.0 - self.0)
        })
    }
}

impl<'a, const N: usize> Sub<&'_ GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn sub(self, rhs: &'_ GF<N>) -> Self::Output {
        GF(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            N - (rhs.0 - self.0)
        })
    }
}

impl<const N: usize> SubAssign<GF<N>> for GF<N> {
    fn sub_assign(&mut self, rhs: GF<N>) {
        *self = GF(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            N - (rhs.0 - self.0)
        });
    }
}

impl<const N: usize> SubAssign<&'_ GF<N>> for GF<N> {
    fn sub_assign(&mut self, rhs: &'_ GF<N>) {
        *self = GF(if self.0 >= rhs.0 {
            self.0 - rhs.0
        } else {
            N - (rhs.0 - self.0)
        });
    }
}

// == REM ==

impl<const N: usize> Rem<GF<N>> for GF<N> {
    type Output = GF<N>;
    fn rem(self, rhs: GF<N>) -> Self::Output {
        GF((self.0 % rhs.0) % N)
    }
}

impl<'a, const N: usize> Rem<GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn rem(self, rhs: GF<N>) -> Self::Output {
        GF((self.0 % rhs.0) % N)
    }
}

impl<const N: usize> Rem<&'_ GF<N>> for GF<N> {
    type Output = GF<N>;
    fn rem(self, rhs: &'_ GF<N>) -> Self::Output {
        GF((self.0 % rhs.0) % N)
    }
}

impl<'a, const N: usize> Rem<&'_ GF<N>> for &'a GF<N> {
    type Output = GF<N>;
    fn rem(self, rhs: &'_ GF<N>) -> Self::Output {
        GF((self.0 % rhs.0) % N)
    }
}

impl<const N: usize> RemAssign<GF<N>> for GF<N> {
    fn rem_assign(&mut self, rhs: GF<N>) {
        *self = GF((self.0 % rhs.0) % N);
    }
}

impl<const N: usize> RemAssign<&'_ GF<N>> for GF<N> {
    fn rem_assign(&mut self, rhs: &'_ GF<N>) {
        *self = GF((self.0 % rhs.0) % N);
    }
}

// == NUM ==

impl<const N: usize> num_traits::Zero for GF<N> {
    fn zero() -> Self {
        GF(0)
    }

    fn is_zero(&self) -> bool {
        self.0 == 0
    }
}

impl<const N: usize> num_traits::One for GF<N> {
    fn one() -> Self {
        GF(1)
    }
}

impl<const N: usize> num_traits::Num for GF<N> {
    type FromStrRadixErr = ParseIntError;

    fn from_str_radix(src: &str, radix: u32) -> Result<Self, Self::FromStrRadixErr> {
        Ok(GF(usize::from_str_radix(src, radix)? % N))
    }
}
