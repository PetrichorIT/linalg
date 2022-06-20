use crate::{matrix::Matrix, num::*};
use num_traits::Float;
use std::{
    cmp,
    ops::{AddAssign, SubAssign},
};

///
/// Computes the eigenvalues of a given matrix using
/// the qr algorithm.
///
/// This function applies the QR-Algorithm to the given square matrix,
/// assumming all eigenvalues have no complex component (this is guaranteed with symmeric matrices).
/// note that the number of iterations is closely corrolated to the expected error of the computation.
/// As a rule of thumb, use the dimension of the matrix as iteration count.
///
/// # Panics
///
/// This functions panics should the number of iterations not be positiv,
/// or should the matrix be non-square.
///
#[deprecated]
pub fn eig_old<T>(mut matrix: Matrix<T>, itr: usize) -> Vec<T>
where
    T: Float + Copy + NumConstants,
{
    assert!(itr > 0);
    assert!(matrix.layout().is_square());

    for _i in 0..itr {
        let q = super::_qr_impl(&mut matrix);
        matrix = Matrix::mmul(&matrix, &q);
    }

    let mut values = Vec::with_capacity(matrix.layout().rows());
    for i in 0..matrix.layout().rows() {
        values.push(matrix[(i, i)])
    }

    values
}

///
/// Performs a limited power iteration using the second operand as start vector.
///
/// # Panics
///
/// This function panics should the first operand be a non-square matrix,
/// or should the second operand not be a colvec, or if the third operand
/// is not real-positiv.
///  
pub fn eigv_powitr<T>(matrix: &Matrix<T>, x: &Matrix<T>, itr: usize) -> (T, Matrix<T>)
where
    T: Float + Copy,
{
    assert!(itr > 0);
    assert!(matrix.layout().is_square());
    assert!(matrix.layout().is_colvec());

    let mut x = Matrix::mmul(matrix, x);
    for _ in 1..itr {
        x = Matrix::mmul(matrix, &x)
    }

    let norm = x.iter().fold(T::zero(), |acc, &c| acc + c * c).sqrt();
    let y = x / norm;

    let ay = Matrix::mmul(matrix, &y);
    let mut l = T::zero();
    for i in 0..y.size() {
        l = l + y[i] * ay[i]
    }

    (l, y)
}

///
/// Performs a limited inverse power iteration using the third operand as start vector.
///
/// Note that this function may not return any value if the inverse matrix cannot
/// be constructed using the given l (this generally only happens if the matrix itself is non-invertable).
///
/// # Panics
///
/// This function panics should the second operand be a non-square matrix,
/// or should the third operand not be a colvec, or if the fourth operand
/// is not real-positiv, or the first operand be zero.
///  
pub fn eigv_powitr_inv<T>(
    l: T,
    matrix: &Matrix<T>,
    x: &Matrix<T>,
    itr: usize,
) -> Option<(T, Matrix<T>)>
where
    T: Float + Copy,
{
    assert!(l != T::zero());
    assert!(itr > 0);
    assert!(matrix.layout().is_square());
    assert!(matrix.layout().is_colvec());

    let a = matrix - &Matrix::eye(matrix.layout().rows()).scalar(l);
    let a = super::inv(a)?;

    Some(eigv_powitr(&a, x, itr))
}

/*
This code is a modification of the implementation from Stepan Yakovenko (https://github.com/stiv-yakovenko).
See https://github.com/stiv-yakovenko/rust-eigen/blob/master/eigen.rs.
*/

fn cdiv<F: Float>(xr: F, xi: F, yr: F, yi: F) -> (F, F) {
    let r: F;
    let d: F;
    if yr.abs() > yi.abs() {
        r = yi / yr;
        d = yr + r * yi;
        ((xr + r * xi) / d, (xi - r * xr) / d)
    } else {
        r = yr / yi;
        d = yi + r * yr;
        ((r * xr + xi) / d, (r * xi - xr) / d)
    }
}

///  This is derived from the Algol procedure hqr2,
///  by Martin and Wilkinson, Handbook for Auto. Comp.,
///  Vol.ii-Linear Algebra, and the corresponding
///  Fortran subroutine in EISPACK.
pub fn hqr2<F: HQR2Ready>(
    n_in: usize,
    h: &mut Matrix<F>,
    v: &mut Matrix<F>,
    d: &mut Vec<F>,
    e: &mut Vec<F>,
) {
    //  This is derived from the Algol procedure hqr2,
    //  by Martin and Wilkinson, Handbook for Auto. Comp.,
    //  Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
    // Initialize
    let nn = n_in;
    let mut n = nn as i16 - 1;
    let low = 0;
    let high = nn - 1;
    let eps = F::epsilon();
    let mut exshift = F::zero();
    let mut p = F::zero();
    let mut q = F::zero();
    let mut r = F::zero();
    let mut s = F::zero();
    let mut z = F::zero();
    let mut t;
    let mut w;
    let mut x;
    let mut y;
    // Store roots isolated by balanc and compute matrix norm
    let mut norm = F::zero();
    let mut i = 0 as usize;
    while i < nn {
        if i < low || i > high {
            d[i] = h[[i, i]];
            e[i] = F::zero();
        }
        let mut j = (i as i16 - 1).max(0) as usize;
        while j < nn {
            norm = norm + (h[[i, j]]).abs();
            j = j + 1;
        }
        i = i + 1;
    }
    // Outer loop over eigenvalue index
    let mut iter = 0;
    while n >= low as i16 {
        // Look for single small sub-diagonal element
        let mut l = n;
        while l > low as i16 {
            s = (h[[l as usize - 1, l as usize - 1]]).abs() + (h[[l as usize, l as usize]]).abs();
            if s == F::zero() {
                s = norm;
            }
            if (h[[l as usize, l as usize - 1]]).abs() < eps * s {
                break;
            }
            l = l - 1;
        }
        // Check for convergence
        // One root found
        if l == n {
            h[[n as usize, n as usize]] = h[[n as usize, n as usize]] + exshift;
            d[n as usize] = h[[n as usize, n as usize]];
            e[n as usize] = F::zero();
            n = n - 1;
            iter = 0;
            // Two roots found
        } else if l == n - 1 {
            w = h[[n as usize, n as usize - 1]] * h[[n as usize - 1, n as usize]];
            p = (h[[n as usize - 1, n as usize - 1]] - h[[n as usize, n as usize]]) / F::two();
            q = p * p + w;
            z = (q).abs().sqrt();
            h[[n as usize, n as usize]] = h[[n as usize, n as usize]] + exshift;
            h[[n as usize - 1, n as usize - 1]] = h[[n as usize - 1, n as usize - 1]] + exshift;
            x = h[[n as usize, n as usize]];
            // Real pair
            if q >= F::zero() {
                if p >= F::zero() {
                    z = p + z;
                } else {
                    z = p - z;
                }
                d[n as usize - 1] = x + z;
                d[n as usize] = d[n as usize - 1];
                if z != F::zero() {
                    d[n as usize] = x - w / z;
                }
                e[n as usize - 1] = F::zero();
                e[n as usize] = F::zero();
                x = h[[n as usize, n as usize - 1]];
                s = (x).abs() + (z).abs();
                p = x / s;
                q = z / s;
                r = (p * p + q * q).sqrt();
                p = p / r;
                q = q / r;
                // Row modification
                let mut j = n - 1;
                while j < nn as i16 {
                    z = h[[n as usize - 1, j as usize]];
                    h[[n as usize - 1, j as usize]] = q * z + p * h[[n as usize, j as usize]];
                    h[[n as usize, j as usize]] = q * h[[n as usize, j as usize]] - p * z;
                    j = j + 1;
                }
                // Column modification
                let mut i = 0;
                while i <= n {
                    z = h[[i as usize, n as usize - 1]];
                    h[[i as usize, n as usize - 1]] = q * z + p * h[[i as usize, n as usize]];
                    h[[i as usize, n as usize]] = q * h[[i as usize, n as usize]] - p * z;
                    i = i + 1;
                }
                // Accumulate transformations
                let mut i = low;
                while i <= high {
                    z = v[[i as usize, n as usize - 1]];
                    v[[i as usize, n as usize - 1]] = q * z + p * v[[i as usize, n as usize]];
                    v[[i as usize, n as usize]] = q * v[[i as usize, n as usize]] - p * z;
                    i = i + 1;
                }
                // Complex pair
            } else {
                d[n as usize - 1] = x + p;
                d[n as usize] = x + p;
                e[n as usize - 1] = z;
                e[n as usize] = -z;
            }
            n = n - 2;
            iter = 0;
            // No convergence yet
        } else {
            // Form shift
            x = h[[n as usize, n as usize]];
            y = F::zero();
            w = F::zero();
            if l < n {
                y = h[[n as usize - 1, n as usize - 1]];
                w = h[[n as usize, n as usize - 1]] * h[[n as usize - 1, n as usize]];
            }
            // Wilkinson's original ad hoc shift
            if iter == 10 {
                exshift += x;
                let mut i = low;
                while i <= n as usize {
                    h[[i, i]] -= x;
                    i = i + 1;
                }
                s = (h[[n as usize, n as usize - 1]]).abs()
                    + (h[[n as usize - 1, n as usize - 2]]).abs();
                y = (F::three_quater()) * s;
                x = y;
                w = F::HQR2_CONST1 * s * s;
            }
            // MATLAB's new ad hoc shift
            if iter == 30 {
                s = (y - x) / F::two();
                s = s * s + w;
                if s > F::zero() {
                    s = s.sqrt();
                    if y < x {
                        s = -s;
                    }
                    s = x - w / ((y - x) / F::two() + s);
                    let mut i = low;
                    while i <= n as usize {
                        h[[i, i]] -= s;
                        i = i + 1;
                    }
                    exshift += s;
                    x = F::HQR2_CONST2;
                    y = x;
                    w = y;
                }
            }
            iter = iter + 1; // (Could check iteration count here.)
                             // Look for two consecutive small sub-diagonal elements
            let mut m = n - 2;
            while m >= l {
                z = h[[m as usize, m as usize]];
                r = x - z;
                s = y - z;
                p = (r * s - w) / h[[m as usize + 1, m as usize]] + h[[m as usize, m as usize + 1]];
                q = h[[m as usize + 1, m as usize + 1]] - z - r - s;
                r = h[[m as usize + 2, m as usize + 1]];
                s = (p).abs() + (q).abs() + (r).abs();
                p = p / s;
                q = q / s;
                r = r / s;
                if m == l {
                    break;
                }
                if h[[m as usize, m as usize - 1]].abs() * (q).abs() + (r).abs()
                    < eps
                        * ((p).abs()
                            * ((h[[m as usize - 1, m as usize - 1]]).abs()
                                + (z).abs()
                                + (h[[m as usize + 1, m as usize + 1]]).abs()))
                {
                    break;
                }
                m = m - 1;
            }
            let mut i = m + 2;
            while i <= n {
                h[[i as usize, i as usize - 2]] = F::zero();
                if i > m + 2 {
                    h[[i as usize, i as usize - 3]] = F::zero();
                }
                i = i + 1;
            }
            // Double QR step involving rows l:n and columns m:n
            let mut k = m;
            while k <= n - 1 {
                let notlast = if k != n - 1 { true } else { false };
                if k != m {
                    p = h[[k as usize, k as usize - 1]];
                    q = h[[k as usize + 1, k as usize - 1]];
                    r = if notlast {
                        h[[k as usize + 2, k as usize - 1]]
                    } else {
                        F::zero()
                    };
                    x = (p).abs() + (q).abs() + (r).abs();
                    if x == F::zero() {
                        k = k + 1;
                        continue;
                    }
                    p = p / x;
                    q = q / x;
                    r = r / x;
                }
                s = (p * p + q * q + r * r).sqrt();
                if p < F::zero() {
                    s = -s;
                }
                if s != F::zero() {
                    if k != m {
                        h[[k as usize, k as usize - 1]] = -s * x;
                    } else if l != m {
                        h[[k as usize, k as usize - 1]] = -h[[k as usize, k as usize - 1]];
                    }
                    p = p + s;
                    x = p / s;
                    y = q / s;
                    z = r / s;
                    q = q / p;
                    r = r / p;
                    // Row modification
                    let mut j = k;
                    while j < nn as i16 {
                        p = h[[k as usize, j as usize]] + q * h[[k as usize + 1, j as usize]];
                        if notlast {
                            p = p + r * h[[k as usize + 2, j as usize]];
                            h[[k as usize + 2, j as usize]] =
                                h[[k as usize + 2, j as usize]] - p * z;
                        }
                        h[[k as usize, j as usize]] = h[[k as usize, j as usize]] - p * x;
                        h[[k as usize + 1, j as usize]] = h[[k as usize + 1, j as usize]] - p * y;
                        j = j + 1;
                    }
                    // Column modification
                    let mut i = 0;
                    while i <= cmp::min(n as usize, k as usize + 3) {
                        p = x * h[[i, k as usize]] + y * h[[i as usize, k as usize + 1]];
                        if notlast {
                            p = p + z * h[[i, k as usize + 2]];
                            h[[i, k as usize + 2]] = h[[i, k as usize + 2]] - p * r;
                        }
                        h[[i, k as usize]] = h[[i, k as usize]] - p;
                        h[[i, k as usize + 1]] = h[[i, k as usize + 1]] - p * q;
                        i = i + 1;
                    }
                    // Accumulate transformations
                    let mut i = low;
                    while i <= high {
                        p = x * v[[i, k as usize]] + y * v[[i, k as usize + 1]];
                        if notlast {
                            p = p + z * v[[i as usize, k as usize + 2]];
                            v[[i as usize, k as usize + 2]] =
                                v[[i as usize, k as usize + 2]] - p * r;
                        }
                        v[[i, k as usize]] = v[[i, k as usize]] - p;
                        v[[i, k as usize + 1]] = v[[i, k as usize + 1]] - p * q;
                        i = i + 1;
                    }
                } // (s != 0)
                k = k + 1;
            } // k loop
        } // check convergence
    } // while n >= low
      // Backsubstitute to find vectors of upper triangular form
    if norm == F::zero() {
        return;
    }
    n = nn as i16 - 1;
    while n >= 0 {
        p = d[n as usize];
        q = e[n as usize];
        // Real vector
        if q == F::zero() {
            let mut l = n;
            h[[n as usize, n as usize]] = F::one();
            let mut i = n as i16 - 1;
            while i >= 0 {
                w = h[[i as usize, i as usize]] - p;
                r = F::zero();
                let mut j = l;
                while j <= n {
                    r = r + h[[i as usize, j as usize]] * h[[j as usize, n as usize]];
                    j = j + 1;
                }
                if e[i as usize] < F::zero() {
                    z = w;
                    s = r;
                } else {
                    l = i;
                    if e[i as usize] == F::zero() {
                        if w != F::zero() {
                            h[[i as usize, n as usize]] = -r / w;
                        } else {
                            h[[i as usize, n as usize]] = -r / (eps * norm);
                        }
                        // Solve real equations
                    } else {
                        x = h[[i as usize, i as usize + 1]];
                        y = h[[i as usize + 1, i as usize]];
                        q = (d[i as usize] - p) * (d[i as usize] - p)
                            + e[i as usize] * e[i as usize];
                        t = (x * s - z * r) / q;
                        h[[i as usize, n as usize]] = t;
                        if (x).abs() > (z).abs() {
                            h[[i as usize + 1, n as usize]] = (-r - w * t) / x;
                        } else {
                            h[[i as usize + 1, n as usize]] = (-s - y * t) / z;
                        }
                    }
                    // Overflow control
                    t = h[[i as usize, n as usize]];
                    if (eps * t).abs() * t > F::one() {
                        let mut j = i;
                        while j <= n as i16 {
                            h[[j as usize, n as usize]] = h[[j as usize, n as usize]] / t;
                            j = j + 1;
                        }
                    }
                }
                i = i - 1;
            }
            // Complex vector
        } else if q < F::zero() {
            let mut l = n - 1;
            // Last vector component imaginary so matrix is triangular
            if (h[[n as usize, n as usize - 1]]).abs() > (h[[n as usize - 1, n as usize]]).abs() {
                h[[n as usize - 1, n as usize - 1]] = q / h[[n as usize, n as usize - 1]];
                h[[n as usize - 1, n as usize]] =
                    -(h[[n as usize, n as usize]] - p) / h[[n as usize, n as usize - 1]];
            } else {
                let (cdivr, cdivi) = cdiv(
                    F::zero(),
                    -h[[n as usize - 1, n as usize]],
                    h[[n as usize - 1, n as usize - 1]] - p,
                    q,
                );
                h[[n as usize - 1, n as usize - 1]] = cdivr;
                h[[n as usize - 1, n as usize]] = cdivi;
            }
            h[[n as usize, n as usize - 1]] = F::zero();
            h[[n as usize, n as usize]] = F::one();
            let mut i = n - 2;
            while i >= 0 {
                let mut ra = F::zero();
                let mut sa = F::zero();
                let mut vr;
                let vi;
                let mut j = l;
                while j <= n {
                    ra += h[[i as usize, j as usize]] * h[[j as usize, n as usize - 1]];
                    sa += h[[i as usize, j as usize]] * h[[j as usize, n as usize]];
                    j = j + 1;
                }
                w = h[[i as usize, i as usize]] - p;
                if e[i as usize] < F::zero() {
                    z = w;
                    r = ra;
                    s = sa;
                } else {
                    l = i;
                    if e[i as usize] == F::zero() {
                        let (cdivr, cdivi) = cdiv(-ra, -sa, w, q);
                        h[[i as usize, n as usize - 1]] = cdivr;
                        h[[i as usize, n as usize]] = cdivi;
                    } else {
                        // Solve complex equations
                        x = h[[i as usize, i as usize + 1]];
                        y = h[[i as usize + 1, i as usize]];
                        vr = (d[i as usize] - p) * (d[i as usize] - p)
                            + e[i as usize] * e[i as usize]
                            - q * q;
                        vi = (d[i as usize] - p) * F::two() * q;
                        if vr == F::zero() && vi == F::zero() {
                            vr = eps
                                * norm
                                * ((w).abs() + (q).abs() + (x).abs() + (y).abs() + (z)).abs();
                        }
                        let (cdivr, cdivi) =
                            cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                        h[[i as usize, n as usize - 1]] = cdivr;
                        h[[i as usize, n as usize]] = cdivi;
                        if (x).abs() > ((z).abs() + (q).abs()) {
                            h[[i as usize + 1, n as usize - 1]] = (-ra
                                - w * h[[i as usize, n as usize - 1]]
                                + q * h[[i as usize, n as usize]])
                                / x;
                            h[[i as usize + 1, n as usize]] = (-sa
                                - w * h[[i as usize, n as usize]]
                                - q * h[[i as usize, n as usize - 1]])
                                / x;
                        } else {
                            let (cdivr, cdivi) = cdiv(
                                -r - y * h[[i as usize, n as usize - 1]],
                                -s - y * h[[i as usize, n as usize]],
                                z,
                                q,
                            );
                            h[[i as usize + 1, n as usize - 1]] = cdivr;
                            h[[i as usize + 1, n as usize]] = cdivi;
                        }
                    }
                    // Overflow control
                    t = (h[[i as usize, n as usize - 1]])
                        .abs()
                        .max(h[[i as usize, n as usize]].abs());
                    if (eps * t) * t > F::one() {
                        let mut j = i;
                        while j <= n {
                            j = j + 1;
                            h[[j as usize, n as usize - 1]] = h[[j as usize, n as usize - 1]] / t;
                            h[[j as usize, n as usize]] = h[[j as usize, n as usize]] / t;
                        }
                    }
                }
                i = i - 1;
            }
        }
        n = n - 1;
    }
    // Vectors of isolated roots
    let mut i = 0;
    while i < nn {
        if i < low || i > high {
            let mut j = i;
            while j < nn {
                v[[i, j]] = h[[i, j]];
                j = j + 1;
            }
        }
        i = i + 1;
    }
    // Back transformation to get eigenvectors of original matrix
    let mut j = nn as i16 - 1;
    while j >= low as i16 {
        let mut i = low;
        while i <= high {
            z = F::zero();
            let mut k = low;
            while k <= (j as usize).min(high) {
                z = z + v[[i, k]] * h[[k, j as usize]];
                k = k + 1;
            }
            v[[i, j as usize]] = z;
            i = i + 1;
        }
        j = j - 1;
    }
}

///  This is derived from the Algol procedures orthes and ortran,
///  by Martin and Wilkinson, Handbook for Auto. Comp.,
///  Vol.ii-Linear Algebra, and the corresponding
///  Fortran subroutines in EISPACK.
pub fn orthes<F: Float + AddAssign + SubAssign>(
    m: &mut Matrix<F>,
    h_mat: &mut Matrix<F>,
    v_mat: &mut Matrix<F>,
) {
    let low = 0;
    let n = m.layout().cols();
    let high = n - 1;
    let mut m = low + 1;
    let mut ort = vec![F::zero(); n];
    while m < high - 1 {
        // Scale column.
        let mut scale = F::zero();
        let mut i = m;
        //for (int        i = m;        i < = high;        i + +)
        while i <= high {
            scale += (h_mat[[i, m - 1]]).abs();
            i = i + 1;
        }
        if scale != F::zero() {
            // Compute Householder transformation.
            let mut h = F::zero();
            let mut i = high;
            while i >= m {
                ort[i] = h_mat[[i, m - 1]] / scale;
                h += ort[i] * ort[i];
                i = i - 1;
            }
            let mut g = h.sqrt();
            if ort[m] > F::zero() {
                g = -g;
            }
            h = h - ort[m] * g;
            ort[m] = ort[m] - g;
            // Apply Householder similarity transformation
            // H = (I-u*u'/h)*H*(I-u*u')/h)
            let mut j = m;
            while j < n {
                let mut f = F::zero();
                let mut i = high;
                while i >= m {
                    f += ort[i] * h_mat[[i, j]];
                    i = i - 1;
                }
                f = f / h;
                let mut i = m;
                while i <= high {
                    h_mat[[i, j]] -= f * ort[i];
                    i = i + 1;
                }
                j = j + 1;
            }
            let mut i = 0;
            while i <= high {
                let mut f = F::zero();
                let mut j = high;
                while j >= m {
                    f += ort[j] * h_mat[[i, j]];
                    j = j - 1;
                }
                f = f / h;
                let mut j = m;
                while j <= high {
                    h_mat[[i, j]] -= f * ort[j];
                    j = j + 1;
                }
                i = i + 1;
            }
            ort[m] = scale * ort[m];
            h_mat[[m, m - 1]] = scale * g;
        }
        m = m + 1;
    }
    // Accumulate transformations (Algol's ortran).
    for i in 0..n {
        for j in 0..n {
            v_mat[[i, j]] = if i == j { F::one() } else { F::zero() };
        }
    }
    let mut m = high - 1;
    while m >= low + 1 {
        if h_mat[[m, m - 1]] != F::zero() {
            let mut i = m + 1;
            while i <= high {
                ort[i] = h_mat[[i, m - 1]];
                i = i + 1;
            }
            let mut j = m;
            while j <= high {
                let mut g = F::zero();
                let mut i = m;
                while i <= high {
                    g += ort[i] * v_mat[[i, j]];
                    i = i + 1;
                }
                // Double division avoids possible underflow
                g = (g / ort[m]) / h_mat[[m, m - 1]];
                let mut i = m;
                while i <= high {
                    v_mat[[i, j]] += g * ort[i];
                    i = i + 1;
                }
                j = j + 1;
            }
        }
        m = m - 1;
    }
}

///
/// A collection of traits and constanst used for the
/// subroutine HQR2.
///
pub trait HQR2Ready: Float + AddAssign + SubAssign + NumConstants + NumFractionalConstants {
    /// -0.4375
    const HQR2_CONST1: Self;
    /// 0.964
    const HQR2_CONST2: Self;
}

impl HQR2Ready for f32 {
    const HQR2_CONST1: Self = -0.4375_f32;
    const HQR2_CONST2: Self = 0.964_f32;
}

impl HQR2Ready for f64 {
    const HQR2_CONST1: Self = -0.4375_f64;
    const HQR2_CONST2: Self = 0.964_f64;
}

///
/// Computes the eigen vaues of a matrix using HQR2 and orthes
/// by Martin and Wilkinson, Handbook for Auto. Comp.,
//  Vol.ii-Linear Algebra, and the corresponding
//  Fortran subroutines in EISPACK.
///
pub fn eig<F: HQR2Ready>(matrix: Matrix<F>) -> Vec<(F, F)> {
    let mut matrix = matrix;
    _eig_impl(&mut matrix)
}

pub(crate) fn _eig_impl<F: HQR2Ready>(m: &mut Matrix<F>) -> Vec<(F, F)> {
    let n = m.layout().cols();
    let mut h_mat = Matrix::new(
        (m.layout().cols(), m.layout().cols()),
        vec![F::zero(); m.layout().cols() * m.layout().cols()],
    );
    let mut v_mat = Matrix::new(
        (m.layout().cols(), m.layout().cols()),
        vec![F::zero(); m.layout().cols() * m.layout().cols()],
    );
    //let mut ort = vec!(0.;n);
    let mut d = vec![F::zero(); n];
    let mut e = vec![F::zero(); n];
    for i in 0..n {
        for j in 0..n {
            h_mat[[i, j]] = m[[i, j]];
        }
    }
    orthes(m, &mut h_mat, &mut v_mat);
    hqr2(n, &mut h_mat, &mut v_mat, &mut d, &mut e);
    let mut r = vec![(F::zero(), F::zero()); n];
    for i in 0..n {
        r[i] = (d[i], e[i])
    }
    r
}

///
/// Solves a polynome using eigenvalues computed by [eig].
///
pub fn solve_poly(c: Vec<f64>) -> Vec<(f64, f64)> {
    let n = c.len();
    let mut m = Matrix::new((n, n), vec![0.0; n * n]);
    for i in 0..(n - 1) {
        m[[i + 1, i]] = 1.;
    }
    for i in 0..(n) {
        m[[i, n - 1]] = -c[i];
    }
    println!("{:?}", m);
    let r = _eig_impl(&mut m);
    r
}

#[cfg(test)]
mod tests {
    #[test]
    #[allow(deprecated)]
    fn test_eig() {
        use super::*;
        use crate::prelude::*;

        let m = matrix![
             3.0,  -1.0,  0.0;
             2.0,   0.0,  0.0;
             -2.0,  2.0, -1.0;
        ];

        let a = eig_old(m.clone(), 4);
        let b = eig(m);

        dbg!(a, b);
    }
}
