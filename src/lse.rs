use std::fmt::Display;

use num_traits::Num;

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

impl LrDecomposition<f64> {
    ///
    /// Creates a LR-Decomposition with a given matrix A.
    /// Will return None if the LSE is not expilictly solvable, meaning
    /// no valid pivot was found.
    ///
    /// # Panics
    ///
    /// Panics of applied to a non-square matrix.
    ///
    pub fn create(a: Matrix<f64>) -> Option<Self> {
        assert!(a.layout().is_square());

        let mut p = Matrix::eye(a.layout().rows());
        let mut l = Matrix::eye(a.layout().rows());

        let mut r = a;

        for i in 0..(r.layout().cols() - 1) {
            // Pivot check
            if r[(i, i)] == 0.0 {
                // search for valid pivot
                let pv_row = (i + 1..r.layout().cols()).find(|&k| r[(k, k)] != 0.0)?;

                p[(i, i)] = 0.0;
                p[(i, pv_row)] = 1.0;
                p[(pv_row, i)] = 1.0;
                p[(pv_row, pv_row)] = 0.0;

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
                    r[(k, j)] -= a * r[(i, j)];
                }
            }
        }

        Some(LrDecomposition { l, r, p })
    }

    pub fn solve(&self, b: Matrix<f64>) -> Matrix<f64> {
        assert!(b.layout().is_colvec());
        assert!(b.layout().rows() == self.l.layout().rows());

        // Apply rowswap
        let mut pb = Matrix::fill(b.layout().clone(), 0.0);
        for row in 0..self.p.layout().rows() {
            for j in 0..self.p.layout().size() {
                if self.p[(row, j)] == 1.0 {
                    pb[row] = b[j];
                    break;
                }
            }
        }

        // Ly = b
        // Rx = y

        // Apply L
        let mut y = Matrix::fill(b.layout().clone(), 0.0);
        for i in 0..self.l.layout().rows() {
            let mut s = 0.0;
            for j in 0..i {
                s += self.l[(i, j)] * y[j]
            }
            y[i] = (pb[i] - s) / self.l[(i, i)];
        }

        let mut x = Matrix::fill(b.layout().clone(), 0.0);
        for i in (0..self.l.layout().rows()).rev() {
            let mut s = 0.0;
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

impl QrDecomposition<f64> {
    pub fn create(a: Matrix<f64>) -> Self {
        let mut r = a;
        let mut q = Matrix::eye(r.layout().rows());

        for c in 0..(r.layout().cols() - 1) {
            let mut v = Matrix::fill(MatrixLayout::new(r.layout().rows(), 1), 0.0);
            for i in c..v.size() {
                v[i] = r[(i, c)];
            }

            let sign = v[c].signum();
            let norm = v.raw().iter().fold(0.0, |acc, &s| acc + s * s).sqrt();

            v[c] += sign * norm;

            let norm = v.raw().iter().fold(0.0, |acc, &s| acc + s * s).sqrt();
            for i in c..v.size() {
                v[i] /= norm;
            }

            let mut h = v.clone() * v.transposed();
            h.scale(2.0);

            let qi = Matrix::eye(r.layout().rows()) - h;

            r = Matrix::mmul(&qi, r);
            q = Matrix::mmul(&qi.transposed(), q);
        }

        // Zero-Fragment cleanup

        for k in 0..r.size() {
            if r[k].abs() < 0.0001 {
                r[k] = 0.0;
            }
        }

        QrDecomposition { q, r }
    }

    pub fn solve(&self, b: Matrix<f64>) -> Matrix<f64> {
        assert!(b.layout().is_colvec());
        assert!(self.q.layout().rows() == b.size());

        // QR = A
        // ==> Qy = b ==> y = Q^T b
        // ==> Rx = y

        // preinit to prevent borrowing issues
        let mut x = Matrix::fill(b.layout().clone(), 0.0);

        let y = Matrix::mmul(&self.q, b);

        for i in (0..self.r.layout().rows()).rev() {
            let mut s = 0.0;
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
