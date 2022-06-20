//!
//! Algebraic functions.
//!

use crate::{
    matrix::{Matrix, MatrixLayout},
    num::{NumConstants, NumFractionalConstants},
};
use num_traits::{Float, Num};
use std::{
    fmt::Display,
    ops::{Add, Mul},
};

#[allow(clippy::all)]
pub(crate) mod eigen;
pub use eigen::*;

/*** LR DECOMPOSITION ***/

///
/// A mathematical LR-Decomposition of a square matrix A into
/// L and R using pivot method encoded in P.
///
/// The struture is encoded in the formular:
/// L*R = P*A
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LrDecomposition<T: Num> {
    /// A lower triangular matrix enocding the gaussian steps.
    pub l: Matrix<T>,

    /// An upper triangular matrix, the result of the gaussian elimination.
    pub r: Matrix<T>,

    /// A encoder for the pivot changes.
    pub p: Matrix<T>,
}

impl<T: Num> LrDecomposition<T>
where
    T: Copy, // where
             // T: Float,
{
    pub fn solve(&self, b: &Matrix<T>) -> Matrix<T> {
        assert!(b.layout().is_colvec());
        assert!(b.layout().rows() == self.l.layout().rows());

        // Apply rowswap
        let mut pb = Matrix::zeroed(b.layout());
        for row in 0..self.p.layout().rows() {
            for j in 0..self.p.layout().size() {
                if self.p[(row, j)] == T::one() {
                    pb[row] = b[j];
                    break;
                }
            }
        }

        // Ly = b
        // Rx = y

        // Apply L
        let mut y = Matrix::zeroed(b.layout());
        for i in 0..self.l.layout().rows() {
            let mut s = T::zero();
            for j in 0..i {
                s = s + self.l[(i, j)] * y[j]
            }
            y[i] = (pb[i] - s) / self.l[(i, i)];
        }

        let mut x = Matrix::zeroed(b.layout());
        for i in (0..self.l.layout().rows()).rev() {
            let mut s = T::zero();
            for j in (i + 1)..self.l.layout().cols() {
                s = s + self.r[(i, j)] * x[j]
            }

            x[i] = (y[i] - s) / self.r[(i, i)];
        }

        x
    }
}

impl<T: Num> Display for LrDecomposition<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "L:\n {}\nR:\n{}\nP:\n{}\n", self.l, self.r, self.p)
    }
}

///
/// Creates a LR-Decomposition with a given matrix A.
///
/// Note that thus function will return None if the
/// LSE is not expilictly solvable, meaning
/// no valid pivot was found in one of its iterations.
///
/// # Panics
///
/// This function will panic when applied to a non-square matrix.
///
/// # Examples
///
/// ```
/// use linalg::prelude::*;
///
/// let matrix = matrix![
///     1.0, 2.0, 3.0;
///     2.0, 3.0, 4.0;
///     3.0, 4.0, 5.0;
/// ];
/// let lrdc = lr(matrix).unwrap();
///
/// assert_eq!(lrdc.r, matrix![
///     1.0,  2.0,  3.0;
///     0.0, -1.0, -2.0;
///     0.0,  0.0,  0.0;
/// ]);
///
/// assert!(lrdc.r.is_upper_triag());
/// assert!(lrdc.l.is_lower_triag());
///
/// // This example did not need to change another pivot.
/// assert_eq!(lrdc.p, Matrix::eye(3));
/// ```
///
pub fn lr<T>(mut matrix: Matrix<T>) -> Option<LrDecomposition<T>>
where
    T: Num + Copy,
{
    let (l, p) = _lr_impl(&mut matrix)?;
    Some(LrDecomposition { l, r: matrix, p })
}

///
/// Uses the gaussian elminiation process (lr process) to
/// generte a lower triangular matrix in place.
///
pub fn lr_dry<T>(a: &mut Matrix<T>) -> Option<()>
where
    T: Num + Copy,
{
    _lr_impl(a)?;
    Some(())
}

/// Internal function to perform a inplace reduction for ref bound values.
fn _lr_impl<T>(r: &mut Matrix<T>) -> Option<(Matrix<T>, Matrix<T>)>
where
    T: Num + Copy,
{
    assert!(r.layout().is_square());

    let mut p = Matrix::eye(r.layout().rows());
    let mut l = Matrix::eye(r.layout().rows());

    for i in 0..(r.layout().cols() - 1) {
        // Pivot check
        if r[(i, i)] == T::zero() {
            // search for valid pivot
            let pv_row = (i + 1..r.layout().cols()).find(|&k| r[(k, k)] != T::zero())?;

            p[(i, i)] = T::zero();
            p[(i, pv_row)] = T::one();
            p[(pv_row, i)] = T::one();
            p[(pv_row, pv_row)] = T::zero();

            // swap rows in r
            for j in 0..r.layout().cols() {
                let temp = r[(i, j)];
                r[(i, j)] = r[(pv_row, j)];
                r[(pv_row, j)] = temp;
            }

            // swap rows in l
            for j in 0..i {
                let temp = l[(i, j)];
                l[(i, j)] = l[(pv_row, j)];
                l[(pv_row, j)] = temp;
            }
        }

        // assert!(r[(i, i)] != 0.0)

        // Reduction
        for k in (i + 1)..r.layout().rows() {
            let a = r[(k, i)] / r[(i, i)];
            l[(k, i)] = a;

            for j in 0..r.layout().cols() {
                let d = a * r[(i, j)];
                r[(k, j)] = r[(k, j)] - d;
            }
        }
    }

    Some((l, p))
}

/*** QR DECOMPOSITION ***/

///
/// An algebraic QR-Decomposition of an abitraty matrix A.
///
/// ```math
/// QR = A
/// ```
#[derive(Debug, Clone)]
pub struct QrDecomposition<T: Num> {
    /// An othogonal matrix that represents the decomposition steps.
    pub q: Matrix<T>,
    /// A upper triangualr matrix, the result of the decomposition using Q.
    pub r: Matrix<T>,
}

impl<T> QrDecomposition<T>
where
    T: Float + Mul<Output = T> + Add<Output = T>,
{
    ///
    /// Uses the given decomposition and the argument to solve the matrix equation:
    /// ```math
    /// Ax = b <-> QRx = b <-> Rx = y and Qy = b
    /// ```
    ///
    /// # Panics
    ///
    /// This function panics should b not be a collum vector with
    /// the same amount of elements as q has rows (matching the linear sys. of eq.).
    /// This amount is equal to the number of collum in the inital matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let qrdc = qr(matrix![
    ///     1.0, 2.0, 3.0;
    ///     2.0, 3.0, 4.0;
    ///     3.0, 4.0, 5.0;
    /// ]);
    ///
    /// let x = qrdc.solve(&matrix![9.0; 8.0; 7.0;]);
    /// ```
    pub fn solve(&self, b: &Matrix<T>) -> Matrix<T> {
        assert!(b.layout().is_colvec());
        assert!(self.q.layout().rows() == b.size());

        // QR = A
        // ==> Qy = b ==> y = Q^T b
        // ==> Rx = y

        // preinit to prevent borrowing issues
        let mut x = Matrix::zeroed(b.layout());

        let y = Matrix::mmul(&self.q, b);

        for i in (0..self.r.layout().rows()).rev() {
            let mut s = T::zero();
            for j in (i + 1)..self.r.layout().cols() {
                s = s + self.r[(i, j)] * x[j]
            }

            x[i] = (y[i] - s) / self.r[(i, i)];
        }

        x
    }
}

impl<T: Num> Display for QrDecomposition<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "=== Q ===\n{}\n=== R ===\n{}\n", self.q, self.r)
    }
}

///
/// A QR-decomposition from a general matrix A into Q and R.
///
/// This function performs a qr-decomposition using the householder transformation
/// thereby consuming the matrix given to the function.
/// What results will be two matrices Q and R where:
/// - Q is othogonal (Q^T = Q^-1) and
/// - R a upper triagonal matrix.
///
/// # Guarntees
///
/// Even though the householder-qr algorithm is numericly prone to rounding errors
/// the upper triagonal form is explicitly enforced by overwriting any non-standart values
/// zero values of the given type T.
///
/// # Examples
///
/// ```
/// use linalg::prelude::*;
///
/// let matrix = matrix![
///     2.0,  4.0, -4.0;
///     1.0,  1.0,  2.0;
///     2.0, -3.0,  0.0;
/// ];
///
/// let decomp = qr(matrix.clone());
/// assert!(decomp.r.is_upper_triag());
///
/// let b = matrix![ 1.0; 2.0; 3.0; ];
///
/// // Solving Ax = b as QRx = b
/// let x = decomp.solve(&b);
///
/// // Reconstructing b from Ax
/// let b2 = Matrix::mmul(&matrix, &x);
///
/// // To account for numeric margins not equality but relativ equalit is needed
/// // A - B = C and C almost 0
/// let diff = b - b2;
/// assert!(diff.iter().fold(0.0, |x, &y| x + y) < 0.1);
/// ```
///
pub fn qr<T>(mut matrix: Matrix<T>) -> QrDecomposition<T>
where
    T: Float + Copy + NumConstants,
{
    let q = _qr_impl(&mut matrix);
    QrDecomposition { q, r: matrix }
}

/// Internal function for qr decomposition.
///
/// This functions performs an in-place qr decomposition
/// on the given mutable ref matrix R.
/// Additionaly it returns the Q matrix.
///
fn _qr_impl<T>(r: &mut Matrix<T>) -> Matrix<T>
where
    T: Float + Copy + NumConstants,
{
    let mut q = Matrix::eye(r.layout().rows());

    for c in 0..(r.layout().cols() - 1) {
        let mut v = Matrix::zeroed(MatrixLayout::new(r.layout().rows(), 1));
        for i in c..v.size() {
            v[i] = r[(i, c)];
        }

        let sign = v[c].signum();
        let norm = v.raw().iter().fold(T::zero(), |acc, &s| acc + s * s).sqrt();

        v[c] = v[c] + sign * norm;

        let norm = v.raw().iter().fold(T::zero(), |acc, &s| acc + s * s).sqrt();
        for i in c..v.size() {
            v[i] = v[i] / norm;
        }

        let mut h = v.clone() * v.transposed();
        h.scale(T::two());

        let qi = Matrix::eye(r.layout().rows()) - h;

        *r = Matrix::mmul(&qi, r);
        q = Matrix::mmul(&qi.transposed(), &q);
    }

    // Zero-Fragment cleanup

    for i in 1..r.layout().rows() {
        for j in 0..i {
            r[(i, j)] = T::zero();
        }
    }

    q
}

/*** DET  ***/

///
/// Calculates the determinat of a matrix using gaussian elminiation.
///
/// # Panics
///
/// This function will panic if appied to a non-square matrix.
///
/// # Examples
///
/// ```
/// use linalg::prelude::*;
///
/// let matrix = matrix![
///     1, 2, 1;
///     3, 4, 5;
///     5, 6, 7;
/// ];
///
/// // Following Leibniz rule the 2x2 matrix should have the following determinat.
/// let d = det(&matrix);
/// assert_eq!(d, 4);
/// ```
pub fn det<T>(matrix: &Matrix<T>) -> T
where
    T: Num + Copy,
{
    // Implement quick fixes for 1x1 2x2 matrices.
    assert!(matrix.layout().is_square());
    match matrix.layout().rows() {
        0 => T::zero(),
        1 => matrix[0],
        2 => matrix[0] * matrix[3] - matrix[1] * matrix[2],
        _ => {
            let mut matrix = matrix.clone();
            lr_dry(&mut matrix).unwrap();

            let mut product = T::one();
            for i in 0..matrix.layout().cols() {
                product = product * matrix[(i, i)]
            }

            product
        }
    }
}

///
/// Checks the given vectors for linear identpendence.
///
pub fn lua<T: Num + Copy>(set: &[Matrix<T>]) -> bool {
    // Assume its col_vecs
    let matrix = Matrix::cloned_from_parts_vertical(set);
    !self::det(&matrix).is_zero()
}

///
/// Returns the trace of the given square matrix.
///
/// # Panics
///
/// This function panics if applied to a non-square matrix.
///
/// # Examples
///
/// ```
/// use linalg::prelude::*;
///
/// let matrix = Matrix::diag(vec![ 1, 2, 3 ]);
/// let tr = trace(&matrix);
/// assert_eq!(tr, 1 + 2 + 3);
/// ```
pub fn trace<T>(matrix: &Matrix<T>) -> T
where
    T: Num + Copy,
{
    assert!(matrix.layout().is_square());
    let mut sum = T::zero();
    for i in 0..matrix.layout().rows() {
        sum = sum + matrix[(i, i)]
    }
    sum
}

/*** Gaus-Jordan Invert ***/

///
/// Returns a inverted matrix of the given matrix if existent or [None] otherwise.
///
/// This function uses the Gauss-Jordan algorithm to solve a system of linear equations,
/// to generate the inverted matrix. If the lse is not solvable, no inverted matrix exists.
///
/// # Panics
///
/// This function panics if applied to a non-square matrix.
///
/// # Example
///
/// ```
/// use linalg::prelude::*;
///
/// let matrix = matrix![
///     2.0, 1.0;
///     6.0, 4.0;
/// ];
///
/// let inverted = inv(matrix.clone()).unwrap();
///
/// assert_eq!(Matrix::mmul(&matrix, &inverted), Matrix::eye(2));
/// ```
pub fn inv<T>(mut matrix: Matrix<T>) -> Option<Matrix<T>>
where
    T: Num + Copy,
{
    assert!(matrix.layout().is_square());

    // TODO:
    // Make memory efficeient method using collum.major matrices

    // SAFTY:
    //  all cells will be filled up

    let n = matrix.layout().rows();
    matrix.resize(MatrixLayout::new(n, 2 * n));

    for i in 0..n {
        for j in 0..n {
            matrix[(i, n + j)] = if i == j { T::one() } else { T::zero() };
        }
    }

    _gauss_jordan_impl(&mut matrix)?;

    // SAFTY:
    // will be filled up
    let mut result = unsafe { Matrix::uninitalized(MatrixLayout::new(n, n)) };

    for i in 0..n {
        for j in 0..n {
            result[(i, j)] = matrix[(i, n + j)];
        }
    }

    Some(result)
}

/// Internal fn to apply the gaus jordan reduction to a given matrix.
/// Panics if cols < rows
fn _gauss_jordan_impl<T>(matrix: &mut Matrix<T>) -> Option<()>
where
    T: Num + Copy,
{
    assert!(matrix.layout().cols() >= matrix.layout().rows());

    for i in 0..matrix.layout().rows() {
        // Find pivot
        if matrix[(i, i)] == T::zero() {
            let pv_row = (i + 1..matrix.layout().rows()).find(|&k| matrix[(k, k)] != T::zero())?;

            // Swap lines
            for j in 0..matrix.layout().cols() {
                let temp = matrix[(i, j)];
                matrix[(i, j)] = matrix[(pv_row, j)];
                matrix[(pv_row, j)] = temp;
            }
        }

        // Normalize line (enforce pv = 1)
        let pv = matrix[(i, i)];
        for j in 0..matrix.layout().cols() {
            matrix[(i, j)] = matrix[(i, j)] / pv;
        }

        // Gauss el
        if i + 1 != matrix.layout().rows() {
            for k in (i + 1)..matrix.layout().rows() {
                let a = matrix[(k, i)];

                for j in 0..matrix.layout().cols() {
                    let d = a * matrix[(i, j)];
                    matrix[(k, j)] = matrix[(k, j)] - d;
                }
            }
        }

        // Jordan elim
        if i != 0 {
            for k in 0..i {
                let a = matrix[(k, i)];

                for j in 0..matrix.layout().cols() {
                    let d = a * matrix[(i, j)];
                    matrix[(k, j)] = matrix[(k, j)] - d;
                }
            }
        }
    }

    Some(())
}

// BROKEN
#[deprecated = "since its broken"]
pub fn tridiag<T>(matrix: Matrix<T>) -> Matrix<T>
where
    T: Float + Copy + NumFractionalConstants,
{
    assert!(matrix.layout().is_square());
    let n = matrix.layout().rows();

    let mut a = T::one();
    for j in 2..n {
        a = a + matrix[(j, 0)] * matrix[(j, 0)];
    }
    a = a.sqrt() * matrix[(1, 0)].signum().neg();

    let r = (T::HALF * (a * a - matrix[(1, 0)] * a)).sqrt();

    let mut v = Matrix::zeroed((n, 1));
    v[0] = T::zero();
    v[1] = (matrix[(1, 0)] - a) / T::two() * r;
    for k in 2..n {
        v[k] = matrix[(k, 0)] / T::two() * r;
    }

    let p1 = Matrix::eye(n) - Matrix::mmul(&v, &v.transposed()).scalar(T::two());
    let mut mt = Matrix::mmul(&p1, &Matrix::mmul(&matrix, &p1));

    for k in 1..(n - 2) {
        let mut a = T::one();
        for j in (k + 1)..n {
            a = a + matrix[(j, k)] * matrix[(j, k)];
        }
        a = a.sqrt() * matrix[(k + 1, k)].signum().neg();

        let r = (T::HALF * (a * a - matrix[(k + 1, k)] * a)).sqrt();
        let mut v = Matrix::zeroed((n, 1));
        // v[k+1]
        v[k] = (matrix[(k + 1, k)] - a) / T::two() * r;

        for j in (k + 2)..n {
            v[j] = matrix[(j, k)] / T::two() * r;
        }

        let p = Matrix::eye(n) - Matrix::mmul(&v, &v.transposed()).scalar(T::two());
        mt = Matrix::mmul(&p, &Matrix::mmul(&mt, &p));
    }

    mt
}

pub fn vandermonde<K: Num + Copy>(elements: &[K]) -> Matrix<K> {
    let n = elements.len();
    assert!(elements.len() >= 2);

    let mut matrix = Vec::with_capacity(n * n);
    // go by memory order
    for e in elements {
        matrix.push(K::one());
        let mut element = *e;
        matrix.push(element);
        for _ in 2..n {
            element = element * (*e);
            matrix.push(element);
        }
    }

    Matrix::new((n, n), matrix)
}
