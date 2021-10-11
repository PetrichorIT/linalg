use std::{
    fmt::Display,
    ops::{Add, AddAssign, DivAssign, Mul, SubAssign},
};

use num_traits::{Float, Num};

use crate::core::{Matrix, MatrixLayout};

///
/// A mathematical LR-Decomposition of a square matrix A into
/// L and R using pivot method encoded in P.
///
/// The struture is encoded in the formular:
/// L*R = P*A
///
#[derive(Debug, Clone)]
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
    T: Float + AddAssign + SubAssign,
{
    ///
    /// Creates a LR-Decomposition with a given matrix A.
    /// Will return None if the LSE is not expilictly solvable, meaning
    /// no valid pivot was found.
    ///
    /// # Panics
    ///
    /// Panics of applied to a non-square matrix.
    ///
    pub fn create(a: Matrix<T>) -> Option<Self> {
        assert!(a.layout().is_square());

        let mut p = Matrix::eye(a.layout().rows());
        let mut l = Matrix::eye(a.layout().rows());

        let mut r = a;

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
                    r[(k, j)] -= d;
                }
            }
        }

        Some(LrDecomposition { l, r, p })
    }

    pub fn solve(&self, b: Matrix<T>) -> Matrix<T> {
        assert!(b.layout().is_colvec());
        assert!(b.layout().rows() == self.l.layout().rows());

        // Apply rowswap
        let mut pb = Matrix::zeroed(b.layout().clone());
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
        let mut y = Matrix::zeroed(b.layout().clone());
        for i in 0..self.l.layout().rows() {
            let mut s = T::zero();
            for j in 0..i {
                s += self.l[(i, j)] * y[j]
            }
            y[i] = (pb[i] - s) / self.l[(i, i)];
        }

        let mut x = Matrix::zeroed(b.layout().clone());
        for i in (0..self.l.layout().rows()).rev() {
            let mut s = T::zero();
            for j in (i + 1)..self.l.layout().cols() {
                s += self.r[(i, j)] * x[j]
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
/// An algebraic QR-Decomposition of an abitraty matrix A.
///
#[derive(Debug, Clone)]
pub struct QrDecomposition<T: Num> {
    pub q: Matrix<T>,
    pub r: Matrix<T>,
}

impl<T> QrDecomposition<T>
where
    T: Float + AddAssign + DivAssign + Mul<Output = T> + Add<Output = T>,
{
    pub fn create(a: Matrix<T>) -> Self {
        let mut r = a;
        let mut q = Matrix::eye(r.layout().rows());

        for c in 0..(r.layout().cols() - 1) {
            let mut v = Matrix::zeroed(MatrixLayout::new(r.layout().rows(), 1));
            for i in c..v.size() {
                v[i] = r[(i, c)];
            }

            let sign = v[c].signum();
            let norm = v.raw().iter().fold(T::zero(), |acc, &s| acc + s * s).sqrt();

            v[c] += sign * norm;

            let norm = v.raw().iter().fold(T::zero(), |acc, &s| acc + s * s).sqrt();
            for i in c..v.size() {
                v[i] /= norm;
            }

            let mut h = v.clone() * v.transposed();
            h.scale(T::one() + T::one());

            let qi = Matrix::eye(r.layout().rows()) - h;

            r = Matrix::mmul(&qi, r);
            q = Matrix::mmul(&qi.transposed(), q);
        }

        // Zero-Fragment cleanup

        for i in 1..r.layout().rows() {
            for j in 0..i {
                r[(i, j)] = T::zero();
            }
        }

        // for k in 0..r.size() {
        //     if r[k].abs() < 0.0001 {
        //         r[k] = T::zero();
        //     }
        // }

        QrDecomposition { q, r }
    }

    pub fn solve(&self, b: Matrix<T>) -> Matrix<T> {
        assert!(b.layout().is_colvec());
        assert!(self.q.layout().rows() == b.size());

        // QR = A
        // ==> Qy = b ==> y = Q^T b
        // ==> Rx = y

        // preinit to prevent borrowing issues
        let mut x = Matrix::zeroed(b.layout().clone());

        let y = Matrix::mmul(&self.q, b);

        for i in (0..self.r.layout().rows()).rev() {
            let mut s = T::zero();
            for j in (i + 1)..self.r.layout().cols() {
                s += self.r[(i, j)] * x[j]
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
