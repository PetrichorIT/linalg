//! A collection of function usfull for discrete mathematics.

use num_traits::{Num, Zero};

///
/// Computes the faculty of the given operand n.
///
/// # Panics
///
/// This function panics should the operand be negativ.
///
/// # Example
///
/// ```
///use linalg::prelude::*;
///
/// let f = fac(6);
/// assert_eq!(f, 6*5*4*3*2*1);
/// ```
///
#[inline(always)]
pub fn fac<T: Num + Copy + PartialOrd>(n: T) -> T {
    prt_under(n, n)
}

///
/// Computes n^k_(under) the decresing partial faculty.
///
/// # Panics
///
/// This function panics if k is greater than n, or either n or k are negativ.
///
pub fn prt_under<T: Num + Copy + PartialOrd>(n: T, k: T) -> T {
    assert!(n >= k);
    assert!(n >= T::zero() && k >= T::zero());

    let mut res = n;
    let mut fc = n;

    while fc > n - k + T::one() {
        fc = fc - T::one();
        res = res * fc;
    }

    res
}

///
/// Computes n^k_(over) the decresing partial faculty.
///
/// # Panics
///
/// This function panics if either n or k are negativ.
///
pub fn prt_over<T: Num + Copy + PartialOrd>(n: T, k: T) -> T {
    assert!(n >= T::zero() && k >= T::zero());

    let mut res = n;
    let mut fc = n;

    while fc < n + k {
        fc = fc + T::one();
        res = res * fc;
    }

    res
}

///
/// This function computes the number of permutations of the given type.
///
/// # Note
///
/// Thsi function assums typ to be a permutation cycle type of an abitrary
/// n-element set, where n is implicitly defined through the type.
/// It computes the number of different permutation which will match this type
/// in their cycle form.
///
pub fn permut_of_type(typ: &[(usize, usize)]) -> usize {
    let mut numerator = 0;
    let mut n = 0;
    for class in typ {
        numerator *= fac(class.1) * class.0.pow(class.1 as u32);
        n += class.0 * class.1;
    }
    fac(n) / numerator
}

///
/// This function calculates the binomalial coefficent.
///
/// This function calculates n over k,
/// the number if subset of size k of a superset size n,
/// using a recursiv (inlined iterativ) method, based on
/// pascal's triagnle using a memory efficent diagonalized representation.
///
/// # Panics
///
/// Note that this function panics if k is greater that n
///
/// # Example
///
/// ```
/// use linalg::prelude::*;
///
/// let b = binom(4, 2);
/// // Ther are 6 subsets of size 2
/// // {1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}.
/// assert_eq!(b, 6);
/// ```
///
pub fn binom(n: usize, k: usize) -> usize {
    assert!(k <= n);
    // There is a pascal-triag of size NxN
    // The only relevant parts of this triag are above and diagonal-left-up
    // from the anchor point (n, k)
    // that means there are n - k element per col
    // => k*d relevant cells

    if k.is_zero() {
        return 1;
    }

    let d: usize = (n - k) + 1;
    if d.is_zero() {
        return 1;
    }

    let k: usize = k;
    let size: usize = d * (k + 1);
    let mut table = Vec::with_capacity(size);

    // fill first col
    for _ in 0..d {
        table.push(1);
    }

    // iterate
    for j in 1..=k {
        // there are allways d cells
        let offset = d * j;
        // Set diag elements
        table.push(1);
        for i in 1..d {
            table.push(table[offset + i - 1] + table[offset - d + i]);
        }
    }

    table[size - 1]
}

///
/// Computes the stirling number of first type.
///
/// This function computes s(n, k),
/// the number of permutations with k cycles,
/// using a inlined recursiv startegy based on pascals triangle.
///
/// # Panics
///
/// This function panics, should k be greater than n.
///
/// # Examples
///
/// ```
/// use linalg::prelude::*;
///
/// let s1 = stirling1(4, 1);
/// // There are 6 permuts with one full cycle.
/// assert_eq!(s1, 6);
/// ```
///
pub fn stirling1(n: usize, k: usize) -> usize {
    assert!(k <= n);
    // There is a pascal-triag of size NxN
    // The only relevant parts of this triag are above and diagonal-left-up
    // from the anchor point (n, k)
    // that means there are n - k element per col
    // => k*d relevant cells

    if n == k {
        return 1;
    }
    if k.is_zero() {
        return 0;
    }

    let d: usize = (n - k) + 1;
    let k: usize = k;
    let size: usize = d * (k + 1);
    let mut table = Vec::with_capacity(size);

    // fill first col
    table.push(1);
    for _ in 1..d {
        table.push(0);
    }

    // iterate
    for j in 1..=k {
        // there are allways d cells
        let offset = d * j;
        // Set diag elements
        table.push(1);
        for i in 1..d {
            let n = j + i;
            table.push(table[offset + i - 1] * (n - 1) + table[offset - d + i]);
        }
    }

    table[size - 1]
}

///
/// Computes the stirling number of second type.
///
/// This function computes S(n, k),
/// the number of k-partitions of a set of size n,
/// using an inlined recursiv method, based on pascals triangle.
///
/// # Panics
///
/// This functions panics should k be greater than n.
///
/// # Example
///
/// ```
/// use linalg::prelude::*;
///
/// let s2 = stirling2(3, 2);
/// // There are 3 possible partitions:
/// // {(1, 2), (3)}, {(1, 3), (2)}, {(2, 3), (1)}
/// assert_eq!(s2, 3);
/// ```
///  
pub fn stirling2(n: usize, k: usize) -> usize {
    assert!(k <= n);
    // There is a pascal-triag of size NxN
    // The only relevant parts of this triag are above and diagonal-left-up
    // from the anchor point (n, k)
    // that means there are n - k element per col
    // => k*d relevant cells

    if n == k {
        return 1;
    }
    if k.is_zero() {
        return 0;
    }

    let d: usize = (n - k) + 1;
    let k: usize = k;
    let size: usize = d * (k + 1);
    let mut table = Vec::with_capacity(size);

    // fill first col
    table.push(1);
    for _ in 1..d {
        table.push(0);
    }

    // iterate
    for j in 1..=k {
        // there are allways d cells
        let offset = d * j;
        // Set diag elements
        table.push(1);
        for i in 1..d {
            table.push(table[offset + i - 1] * j + table[offset - d + i]);
        }
    }

    table[size - 1]
}
