//! Primitves related to set theory, as well as usefull sets
//! like the real numbers, or integers with modulo arithmetic
//!

// == Type definitions ==

mod finite_integer;

pub use finite_integer::*;
use num_traits::{One, Zero};

/// Natural number from 0 to infinity in discrete steps of size 1.
pub type Nat = usize;
/// Signed natural numbers.
pub type Integer = isize;
/// Real floating point numbers using 32 bits.
pub type Real32 = f32;
/// Real floating point numbers using 64 bits.
pub type Real64 = f64;
/// Complex numbers using 32 bit for each component for a total of 64 bits.
pub type Complex32 = num_complex::Complex32;
/// Complex numbers using 64 bit for each component for a total of 128 bits.
pub type Complex64 = num_complex::Complex64;

// == Marker Traits ==

///
/// A marker trait for an additive Group (G, +) so that:
///
/// (G0) For each a,b of G: (a+b) is in G
///
/// (G1) For each a,b,c of G:
///         (a+b)+c = a+(b+c)
///
/// (G2) There is one and only one element n of G so that:
///         For each a of G:
///             a+n = n+a = a
///
/// (G3) Fo each a of G there exists only one (-a) of G so that:
///         a+(-a) = (-a)+a =e
///
pub trait AdditiveGroup
where
    Self: Zero + PartialEq,
{
}

///
/// A marker trait for an abelian / kommutativ Group (G, +) so that:
///
/// (GA0) (G, +) is a group
///
/// (GA1) For each a,b of G:
///     a+b = b+a
///
pub trait AbelianAdditiveGroup: AdditiveGroup {}

///
/// A marker trait for a multiplicative half group (G, *) so that:
///
/// (H0) For each a,b of G: (a*b) is in G
///
/// (H1) For each a,b,c of G:
///         (a*b)*c = a*(b*c)
///
pub trait MultiplicativeHalfGroup
where
    Self: One,
{
}

///
/// A matrker traif for a ring (G, +, *) so that:
///
/// (R1) (G, +) is an abelian group.
///
/// (R2) (G, *) is a half group (if ablian that the ring aswell).
///
/// (R3) For each a,b,c of G:
///     a*(b+c) = a*b + a*c
///     (a+b)*c = a*c + b*c
///
pub trait Ring: AbelianAdditiveGroup + MultiplicativeHalfGroup {}

///
/// A marker trait for a field (K, +, *) so that:
///
/// (K1) (K, +, *) is a ring.
///
/// (K2) (G, +) is an abelian group with the neural element e.
///
/// (K3) The neutral element for addition n is NOT equal to the
///      neutral element of mutiplication e.
///
pub trait Field: Ring {}

/// A field with copyable elements.
pub trait CField: Field + Copy {}
impl<T: Field + Copy> CField for T {}

///
/// A marker trait that indicates that a set has a finite size.
///
pub trait Finite
where
    Self: Sized,
{
    fn all() -> Vec<Self>;
}

// == Impls ===

macro_rules! impl_trait {
    ($t:ty => $($v:ty),*) => {
        $(
            impl $t for $v {}
        )*
    };
}

impl_trait!(AdditiveGroup => Integer, Real32, Real64, Complex32, Complex64);
impl_trait!(AbelianAdditiveGroup => Integer, Real32, Real64, Complex32, Complex64);
impl_trait!(MultiplicativeHalfGroup =>  Real32, Real64, Complex32, Complex64);
impl_trait!(Ring => Real32, Real64, Complex32, Complex64);
impl_trait!(Field => Real32, Real64, Complex32, Complex64);

// == Generic Impl ==

impl<const N: usize> AdditiveGroup for GF<N> {}
impl<const N: usize> AbelianAdditiveGroup for GF<N> {}
impl<const N: usize> MultiplicativeHalfGroup for GF<N> {}
impl<const N: usize> Ring for GF<N> {}
impl<const N: usize> Field for GF<N> {}
impl<const N: usize> Finite for GF<N> {
    fn all() -> Vec<Self> {
        let mut result = Vec::with_capacity(N);
        for i in 0..N {
            result.push(GF(i))
        }
        result
    }
}
