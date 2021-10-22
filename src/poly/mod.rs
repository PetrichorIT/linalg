//! A set of functions for dealing with polynoms.

use std::{
    fmt::Display,
    hash::Hash,
    ops::{Deref, DerefMut, Index, IndexMut},
};

use num_traits::{Float, Num};

use crate::prelude::{lr, Matrix};

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct Polynom<T> {
    coefficients: Vec<T>,
}

impl<T> Polynom<T> {
    /// The coefficients that define the polynom, using index as power.
    #[inline(always)]
    pub fn coefficients(&self) -> &Vec<T> {
        &self.coefficients
    }

    pub fn new(coefficients: Vec<T>) -> Self {
        assert!(!coefficients.is_empty());
        Self { coefficients }
    }
}

impl<T> Polynom<T>
where
    T: Num + Copy,
{
    ///
    /// The degree of the polynom, meaning the highes index of a relevant coefficient.
    ///
    /// # Examples
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// // y = f(x) = 3x^2 + 2x + 1
    /// let poly = Polynom::new(vec![1, 2, 3, 0, 0]);
    /// assert_eq!(poly.degree(), 2);
    /// ```
    pub fn degree(&self) -> usize {
        for i in 1..=self.coefficients.len() {
            let idx = self.coefficients.len() - i;
            if !self.coefficients[idx].is_zero() {
                return idx;
            }
        }
        return 0;
    }

    ///
    /// Evaluate the polynomial at a given point x and get the resulting y value.
    ///
    /// This functions evaluates all relevant coefficients (up to the degree),
    /// and iterativly computes the x^i to compute the final sum of coefficents*x^i = y.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// // y = f(x) = 3x^2 + 2x + 1
    /// let poly = Polynom::new(vec![1, 2, 3, 0, 0]);
    /// let y = poly.evaluate_at(3);
    ///
    /// assert_eq!(y, 1 + 2*(3) + 3*(3*3));
    /// ```
    pub fn evaluate_at(&self, x: T) -> T {
        let mut sum = T::zero();
        let mut pow = T::one();
        for i in 0..=self.degree() {
            sum = sum + (self.coefficients[i] * pow);
            pow = pow * x;
        }
        sum
    }
}

impl<T: Display> Display for Polynom<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut str = String::new();
        for i in (0..self.coefficients.len()).rev() {
            str.push_str(&format!("{}x^{}", self.coefficients[i], i))
        }
        write!(f, "y = {}", str)
    }
}

impl<T> Default for Polynom<T> {
    fn default() -> Self {
        Polynom {
            coefficients: vec![],
        }
    }
}

impl<T> Deref for Polynom<T> {
    type Target = [T];
    fn deref(&self) -> &Self::Target {
        self.coefficients.deref()
    }
}

impl<T> DerefMut for Polynom<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.coefficients.deref_mut()
    }
}

impl<T: Hash> Hash for Polynom<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.coefficients.hash(state)
    }
}

impl<T> Index<usize> for Polynom<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        self.coefficients.index(index)
    }
}

impl<T> IndexMut<usize> for Polynom<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.coefficients.index_mut(index)
    }
}

//
// ===========
//  Functions
// ===========
//

///
/// Interpolate the points using a lse construction to perfectly intersect all points.
///
/// Note that this function may fail to deliver results if the points contain duplicates (x-values)
/// or a numeric error happens while solving the lse (e.g. to similar x values).
///
pub fn pinterpol<T: Float + Copy>(points: &[(T, T)]) -> Result<Polynom<T>, &'static str> {
    let n = points.len();
    let mut a = Matrix::zeroed((n, n));
    let mut b = Matrix::zeroed((n, 1));

    for i in 0..n {
        let mut x = points[i].0;
        a[(i, 0)] = T::one();
        for j in 1..n {
            a[(i, j)] = x;
            x = x * points[i].0;
        }

        b[i] = points[i].1;
    }

    let lr = lr(a);
    let res = match lr {
        Some(lr) => {
            let x: Vec<T> = lr.solve(&b).into();
            Ok(Polynom::new(x))
        }
        None => Err("Cannot solve lse, due to either numeric error or duplicate points"),
    };

    res
}
