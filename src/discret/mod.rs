//! A collection of function usfull for discrete mathematics.

use num_traits::{Num, Zero};

mod tests;

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
    if n.is_zero() {
        return T::one();
    }
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

    while fc < n + k - T::one() {
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
    let mut numerator = 1;
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
/// A toolset to caculate mutiple bionomial coefficeients efficently
///
/// This module provides the functionallity to caculate binomial coefficents,
/// and memorize the recursive steps, to reuse those steps (saved as a table)
/// to caculate the next coefficeient without duplicated steps.
///
/// ```
/// use linalg::prelude::binom;
/// use linalg::prelude::binom_tbl::*;
///
/// // Generates a table that at least holds the values up to binom(4, 4);
/// let mut table = gen(4, 4);
///
/// // Use constant time operations to get allready caculated values
/// let binom_3_3 = table[idx(3, 3)];
/// assert_eq!(binom_3_3, binom(3, 3));
///
/// // Use expand to increase the table to hold new values (at least (n, n))
/// let (n, k) = (5, 5);
/// expand(&mut table, n);
///
/// let binom_5_5 = table[idx(5, 5)];
/// assert_eq!(binom_5_5, binom(5, 5));
///
/// // Use getter to automaticly expand if nessecary
/// let binom_6_6 = get(&mut table, 6, 6);
/// assert_eq!(binom_6_6, binom(6, 6));
/// ```
pub mod binom_tbl {
    ///
    /// Generates a recursiv table to calculate binomial coefficents.
    ///
    /// Note thath this function will compute all values up to (n_cap, k_cap)
    /// meaning all values with n < n_cap and all values with n == n_cap and k <= k_cap.
    ///
    pub fn gen(n_cap: usize, k_cap: usize) -> Vec<usize> {
        // This function will generate a vec to be interprerted row-major
        // with truncated rows where row i has i+1 elements

        assert!(k_cap <= n_cap);

        let tb_size: usize = (n_cap + 1) * (n_cap + 1);
        let upper_triag: usize = (1..=n_cap).sum();
        let size = tb_size - upper_triag - (n_cap - k_cap);

        let mut result = Vec::with_capacity(size);

        for i in 0..=n_cap {
            result.push(1);
            for j in 1..=i {
                if j == i {
                    result.push(1)
                } else {
                    // (i, j) = (i-1, j-1) + (i-1, j)
                    let idx = result.len() - i;
                    result.push(result[idx - 1] + result[idx]);
                }
            }
        }

        result
    }

    ///
    /// Calculates a bionomial coefficent based on the values in the given table.
    ///
    /// This function either extracts the allready cacluated value from the table
    /// or expands the table until it holds the searched value.
    ///
    pub fn get(table: &mut Vec<usize>, n: usize, k: usize) -> usize {
        assert!(k <= n);
        let idx = idx(n, k);

        if idx >= table.len() {
            expand(table, n)
        }
        table[idx]
    }

    ///
    /// Calculates the index of a entry in the table.
    ///
    #[inline(always)]
    pub fn idx(n: usize, k: usize) -> usize {
        let n_offset: usize = (1..=n).sum();
        n_offset + k
    }

    ///
    /// Caculates the current capacity (thus all calculated values) of
    /// a given table.
    ///
    pub fn cap(table: &Vec<usize>) -> (usize, usize) {
        let mut k_cap = table.len();
        let mut n_cap = 1;
        while k_cap > n_cap {
            k_cap -= n_cap;
            n_cap += 1;
        }

        n_cap -= 1;
        k_cap -= 1;

        (n_cap, k_cap)
    }

    ///
    /// Expands the table to contain values up to (n, n).
    ///
    pub fn expand(table: &mut Vec<usize>, n: usize) {
        // extract n_cap, j_cap to append
        let (n_cap, k_cap) = cap(table);

        // fill up current row
        for j in k_cap..n_cap {
            if n_cap == j + 1 {
                table.push(1)
            } else {
                let idx = table.len() - n_cap;
                table.push(table[idx - 1] + table[idx]);
            }
        }

        for i in (n_cap + 1)..=n {
            table.push(1);
            for j in 1..=i {
                if i == j {
                    table.push(1)
                } else {
                    let idx = table.len() - i;
                    table.push(table[idx - 1] + table[idx]);
                }
            }
        }
    }
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
/// A toolset to caculate mutiple stirling(1) numbers efficently
///
/// This module provides the functionallity to caculate stirling(1) coefficents,
/// and memorize the recursive steps, to reuse those steps (saved as a table)
/// to caculate the next instances without duplicated steps.
///
/// ```
/// use linalg::prelude::stirling1;
/// use linalg::prelude::stirling1_tbl::*;
///
/// // Generates a table that at least holds the values up to stirling1(4, 4);
/// let mut table = gen(4, 4);
///
/// // Use constant time operations to get allready caculated values
/// let stir1_3_3 = table[idx(3, 3)];
/// assert_eq!(stir1_3_3, stirling1(3, 3));
///
/// // Use expand to increase the table to hold new values (at least (n, n))
/// let (n, k) = (5, 5);
/// expand(&mut table, n);
///
/// let stir1_5_5 = table[idx(n, k)];
/// assert_eq!(stir1_5_5, stirling1(n, k));
///
/// // Use getter to automaticly expand if nessecary
/// let stir1_6_6 = get(&mut table, 6, 6);
/// assert_eq!(stir1_6_6, stirling1(6, 6));
/// ```
pub mod stirling1_tbl {
    ///
    /// Generates a new stirling(1) table containing all values up to (n_cap, k_cap).
    ///
    pub fn gen(n_cap: usize, k_cap: usize) -> Vec<usize> {
        // This function will generate a vec to be interprerted row-major
        // with truncated rows where row i has i+1 elements

        assert!(k_cap <= n_cap);

        let tb_size: usize = (n_cap + 1) * (n_cap + 1);
        let upper_triag: usize = (1..=n_cap).sum();
        let size = tb_size - upper_triag - (n_cap - k_cap);

        let mut result = Vec::with_capacity(size);

        result.push(1);
        for i in 1..=n_cap {
            result.push(0);
            for j in 1..=i {
                if j == i {
                    result.push(1)
                } else {
                    // (i, j) = (i-1, j-1) + (i-1, j)
                    let idx = result.len() - i;
                    result.push(result[idx - 1] + result[idx] * (i - 1));
                }
            }
        }

        result
    }

    ///
    /// Calculates a stirling(1) number based on the values in the given table.
    ///
    /// This function either extracts the allready cacluated value from the table
    /// or expands the table until it holds the searched value.
    ///
    pub fn get(table: &mut Vec<usize>, n: usize, k: usize) -> usize {
        assert!(k <= n);
        let idx = idx(n, k);

        if idx >= table.len() {
            expand(table, n)
        }
        table[idx]
    }

    ///
    /// Calculates the index of a entry in the binom table.
    ///
    #[inline(always)]
    pub fn idx(n: usize, k: usize) -> usize {
        let n_offset: usize = (1..=n).sum();
        n_offset + k
    }

    ///
    /// Caculates the current capacity (thus all calculated values) of
    /// a given table.
    ///
    pub fn cap(table: &Vec<usize>) -> (usize, usize) {
        let mut k_cap = table.len();
        let mut n_cap = 1;
        while k_cap > n_cap {
            k_cap -= n_cap;
            n_cap += 1;
        }

        n_cap -= 1;
        k_cap -= 1;

        (n_cap, k_cap)
    }

    ///
    /// Expands the table to contain values up to (n, n).
    ///
    pub fn expand(table: &mut Vec<usize>, n: usize) {
        // extract n_cap, j_cap to append
        let (n_cap, k_cap) = cap(table);

        // fill up current row
        for j in k_cap..n_cap {
            if n_cap == j + 1 {
                table.push(1)
            } else {
                let idx = table.len() - n_cap;
                table.push(table[idx - 1] + table[idx] * (n_cap - 1));
            }
        }

        for i in (n_cap + 1)..=n {
            table.push(0);
            for j in 1..=i {
                if i == j {
                    table.push(1)
                } else {
                    let idx = table.len() - i;
                    table.push(table[idx - 1] + table[idx] * (i - 1));
                }
            }
        }
    }
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

///
/// A toolset to caculate mutiple stirling(2) numbers efficently
///
/// This module provides the functionallity to caculate stirling(2) coefficents,
/// and memorize the recursive steps, to reuse those steps (saved as a table)
/// to caculate the next instances without duplicated steps.
///
/// ```
/// use linalg::prelude::stirling2;
/// use linalg::prelude::stirling2_tbl::*;
///
/// // Generates a table that at least holds the values up to stirling1(4, 4);
/// let mut table = gen(4, 4);
///
/// // Use constant time operations to get allready caculated values
/// let stir2_3_3 = table[idx(3, 3)];
/// assert_eq!(stir2_3_3, stirling2(3, 3));
///
/// // Use expand to increase the table to hold new values (at least (n, n))
/// let (n, k) = (5, 5);
/// expand(&mut table, n);
///
/// let stir2_5_5 = table[idx(n, k)];
/// assert_eq!(stir2_5_5, stirling2(n, k));
///
/// // Use getter to automaticly expand if nessecary
/// let stir2_6_6 = get(&mut table, 6, 6);
/// assert_eq!(stir2_6_6, stirling2(6, 6));
/// ```
pub mod stirling2_tbl {
    ///
    /// Generates a new table containing stirling(2) values up to (n_cap, k_cap).
    ///
    pub fn gen(n_cap: usize, k_cap: usize) -> Vec<usize> {
        // This function will generate a vec to be interprerted row-major
        // with truncated rows where row i has i+1 elements

        assert!(k_cap <= n_cap);

        let tb_size: usize = (n_cap + 1) * (n_cap + 1);
        let upper_triag: usize = (1..=n_cap).sum();
        let size = tb_size - upper_triag - (n_cap - k_cap);

        let mut result = Vec::with_capacity(size);

        result.push(1);
        for i in 1..=n_cap {
            result.push(0);
            for j in 1..=i {
                if j == i {
                    result.push(1)
                } else {
                    // (i, j) = (i-1, j-1) + (i-1, j)
                    let idx = result.len() - i;
                    result.push(result[idx - 1] + result[idx] * j);
                }
            }
        }

        result
    }

    ///
    /// Calculates a stirling(2) numbers based on the values in the given table.
    ///
    /// This function either extracts the allready cacluated value from the table
    /// or expands the table until it holds the searched value.
    ///
    pub fn get(table: &mut Vec<usize>, n: usize, k: usize) -> usize {
        assert!(k <= n);
        let idx = idx(n, k);

        if idx >= table.len() {
            expand(table, n)
        }
        table[idx]
    }

    ///
    /// Calculates the index of a entry in the binom table.
    ///
    #[inline(always)]
    pub fn idx(n: usize, k: usize) -> usize {
        let n_offset: usize = (1..=n).sum();
        n_offset + k
    }

    ///
    /// Caculates the current capacity (thus all calculated values) of
    /// a given table.
    ///
    pub fn cap(table: &Vec<usize>) -> (usize, usize) {
        let mut k_cap = table.len();
        let mut n_cap = 1;
        while k_cap > n_cap {
            k_cap -= n_cap;
            n_cap += 1;
        }

        n_cap -= 1;
        k_cap -= 1;

        (n_cap, k_cap)
    }

    ///
    /// Expands the table to contain values up to (n, n)
    ///
    pub fn expand(table: &mut Vec<usize>, n: usize) {
        // extract n_cap, j_cap to append
        let (n_cap, k_cap) = cap(table);

        // fill up current row
        for j in k_cap..n_cap {
            if n_cap == j + 1 {
                table.push(1)
            } else {
                let idx = table.len() - n_cap;
                table.push(table[idx - 1] + table[idx] * (j + 1));
            }
        }

        for i in (n_cap + 1)..=n {
            table.push(0);
            for j in 1..=i {
                if i == j {
                    table.push(1)
                } else {
                    let idx = table.len() - i;
                    table.push(table[idx - 1] + table[idx] * j);
                }
            }
        }
    }
}

///
/// Computes the bell number of the given operand n.
///
pub fn bell(n: usize) -> usize {
    let mut result = 0;
    for k in 0..=n {
        result = result + stirling2(n, k)
    }
    result
}
