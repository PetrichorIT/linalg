//!
//! A module for solving linear optimization problems.
//!
//! A linear optimization problem in (semi-)normal form consists of a
//! linear cost function f defined by the vector c and a constant offset l,
//! and two sets of constraints.
//! The border constraints are defined by the "a" martix and the "b" vector
//! following the inqeuality Ax <= b.
//! The equality constrains are defined by the "a_eq" matrix and the "b_eq" vector
//! following the equality Ax = b
//!
//! # Solving a LOP
//!
//! When solving a given LOP a simplified simplex algorithm will be used.
//! This may result in a solution of type [LOPSolution] consisting of a
//! solution vector x, its allocated cost value, and processing increments if
//! the solver was used in verbose mode.
//! If no solution is found then the problem is either unbound, or iteration limits
//! are exceeded.
//!

use std::{
    fmt::Display,
    mem::swap,
    ops::{AddAssign, Neg},
};

use num_traits::{Float, Num};

use crate::matrix::{Matrix, MatrixLayout};

// min c^T*x + l with a_eq*x=b_eq and a*x <= b

/// A linear optimization problem type alias for [LinearOpimizationProblem].
pub type LOP<T> = LinearOpimizationProblem<T>;

///
/// A linear optimization problem in normal form.
///
/// This structure describes a minimation problem of type
/// ```math
/// min c^T*x + l
///     with Ax <= b and A{_eq}x <= b{_eq} and x >= 0
/// ```
///
/// This structure also holds the function [LinearOpimizationProblem::solve] to
/// solve the optimization problem an return the caclulated result.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LinearOpimizationProblem<T: Float> {
    /// The non-constant coefficiant vector of the cost function.
    pub c: Matrix<T>,
    /// The constant offset coefficiant of the cost function.
    pub l: T,
    /// The coefficant matrix of the border constraints.
    pub a: Matrix<T>,
    /// The offset vector for the border constraints.
    pub b: Matrix<T>,
    /// The coefficant matrix of the equality constrains.
    pub a_eq: Matrix<T>,
    /// The offset vector of the equality constrains.
    pub b_eq: Matrix<T>,
}

impl<T: Float> LinearOpimizationProblem<T> {
    ///
    /// Creates a new LOP using the given partial values.
    ///
    /// # Panics
    ///
    /// This function will panic should the size of the c vector not equal the
    /// collum dimension of the coefficient matrices, or should the row dimension of each
    /// coefficent matrix not equal the size of its b vector.
    /// Note that matrics of size 0, are exempted from rules defined by the c vector.
    ///
    #[inline]
    pub fn new(
        c: Matrix<T>,
        l: T,
        a: Matrix<T>,
        b: Matrix<T>,
        a_eq: Matrix<T>,
        b_eq: Matrix<T>,
    ) -> Self {
        assert!(c.size() == a.layout().cols() || a.size() == 0);
        assert!(c.size() == a_eq.layout().cols() || a_eq.size() == 0);
        assert_eq!(a.layout().rows(), b.size());
        assert_eq!(a_eq.layout().rows(), b_eq.size());

        Self {
            c,
            l,
            a,
            b,
            a_eq,
            b_eq,
        }
    }

    ///
    /// Creates a new LOP using the given partial values converting
    /// them to matrices or vectors (subcase of matrices) using the
    /// best implementation of [Into].
    ///
    /// # Panics
    ///
    /// This function will panic should either a conversion fails
    /// or any of the constrains layed out in [LinearOpimizationProblem::new()] is not given.
    ///
    #[inline]
    pub fn build<VecType, MatType>(
        c: VecType,
        l: T,
        a: MatType,
        b: VecType,
        a_eq: MatType,
        b_eq: VecType,
    ) -> Self
    where
        VecType: Into<Matrix<T>>,
        MatType: Into<Matrix<T>>,
    {
        Self::new(c.into(), l, a.into(), b.into(), a_eq.into(), b_eq.into())
    }
}

impl<T> LinearOpimizationProblem<T>
where
    T: Float + AddAssign + Neg + Display,
{
    ///
    /// Tries to solve the linear optimization problem using the default
    /// solving configuration.
    ///
    pub fn solve(&self) -> Result<LOPSolution<T>, &'static str> {
        self.solve_with(LOPOptions::default())
    }

    ///
    /// Tries to solve the linear optimization problem using the provided
    /// options and observers.
    ///
    /// This solving algorithm will terminated after at
    /// most [LOPOptions.max_p1_iterations] + [LOPOptions.max_p2_iterations].
    /// Possible errors include max. iteration count exeeded, or and unbound problem (type 1 / 2).
    ///
    /// # Panics
    ///
    /// This function will panic should any of the constrains layed out in
    /// [LinearOpimizationProblem::new()] be violated, or should the c and b vector
    /// sieze to be collumvectors.
    ///
    pub fn solve_with(&self, options: LOPOptions) -> Result<LOPSolution<T>, &'static str> {
        assert!(self.c.layout().is_colvec());
        assert!(self.b.layout().is_colvec());
        assert!(self.b_eq.layout().is_colvec());
        assert!(self.a.layout().rows() == self.b.layout().rows());
        assert!(self.a.layout().cols() == self.c.size() || self.a.layout().cols() == 0);
        assert!(self.a_eq.layout().rows() == self.b_eq.layout().rows());
        assert!(self.a_eq.layout().cols() == self.c.size() || self.a_eq.layout().cols() == 0);

        let n = self.c.layout().rows();
        let number_of_vars = n + self.a.layout().rows() + self.a_eq.layout().rows();
        let number_of_eq = self.a.layout().rows() + self.a_eq.layout().rows();

        let mut increments = Vec::<LOPIncrement<T>>::new();

        // SAFTY:
        // All positions in the given matrix will be filled in the following
        // steps without using the matrix values as getter
        let layout = MatrixLayout::new(number_of_eq + 2, n + 1);
        let mut mat = Matrix::<T>::zeroed(layout);

        // Fill the Axy=b matrix
        if self.a.layout().rows() != 0 {
            mat.insert(0.., n.., &self.b);
            mat.insert(0.., 0.., &self.a);
        }

        // Fill A_eqx=b_eq
        if self.a_eq.layout().rows() != 0 {
            mat.insert(self.b.layout().rows().., n.., &self.b_eq);
            mat.insert(self.b.layout().rows().., 0.., &self.a_eq);
        }

        mat.insert((number_of_eq + 1).., 0.., &self.c.transposed().neg());

        // Sum up all rows of a_eq
        let mut ch = Matrix::<T>::zeroed(MatrixLayout::new(n + 1, 1));
        for o in self.a.layout().rows()..number_of_eq {
            ch += mat.extract(o..=o, 0..).transposed()
        }

        // Fill in ch (transposed) into matrix at anchor (number_of_eq, 0)
        mat.insert(number_of_eq.., 0.., &ch.transposed());

        let mut base_vars = Vec::with_capacity(n);
        for i in 1..=n {
            base_vars.push(i);
        }

        let mut non_base_vars = Vec::with_capacity(number_of_vars - n);
        for i in (n + 1)..=number_of_vars {
            non_base_vars.push(i);
        }

        // Phase 1
        // let mut p1_hash_list = Vec::new();
        let mut p1_itr = 0;

        loop {
            // Get Pivot

            let mut col_select = T::infinity();
            let mut pv_col = None;

            for i in 0..n {
                let raw_index = mat.layout().index((number_of_eq, i));
                if mat[raw_index] <= T::zero() {
                    continue;
                }

                if mat[raw_index] < col_select {
                    col_select = mat[raw_index];
                    pv_col = Some(i);
                }
            }

            if pv_col.is_none() {
                break;
            }
            let pv_col = pv_col.unwrap();

            let mut row_select = T::infinity();
            let mut pv_row = None;

            for j in 0..number_of_eq {
                if mat[(j, pv_col)] == T::zero() {
                    continue;
                }

                let r_value = mat[(j, n)] / mat[(j, pv_col)];
                if r_value < T::zero() {
                    continue;
                }

                if r_value < row_select {
                    row_select = r_value;
                    pv_row = Some(j);
                }
            }

            if pv_row.is_none() {
                break;
            }
            let pv_row = pv_row.unwrap();

            // Got Pivot, Do Iteration

            if p1_itr >= options.max_p1_iterations {
                return Err("Exceeded maximum phase one iteration count");
            }
            p1_itr += 1;

            // TODO: Hasher to prevent endless loops

            // Log current iteration
            if options.verbose {
                increments.push(LOPIncrement {
                    phase: 1,
                    iteration: p1_itr,
                    pv_row,
                    pv_col,
                    base_vars: base_vars.clone(),
                    non_base_vars: non_base_vars.clone(),
                    matrix: mat.clone(),
                })
            }

            let pv = mat[(pv_row, pv_col)];

            // Swap indices
            swap(&mut base_vars[pv_col], &mut non_base_vars[pv_row]);

            let old = mat.clone();
            for row in 0..mat.layout().rows() {
                for col in 0..mat.layout().cols() {
                    match (row == pv_row, col == pv_col) {
                        (true, true) => mat[(row, col)] = T::one() / pv,
                        (true, false) => mat[(row, col)] = old[(row, col)] / pv,
                        (false, true) => mat[(row, col)] = (old[(row, col)] / pv).neg(),
                        (false, false) => {
                            mat[(row, col)] =
                                old[(row, col)] - (old[(pv_row, col)] * old[(row, pv_col)]) / pv
                        }
                    }
                }
            }
        }

        if options.verbose {
            increments.push(LOPIncrement {
                phase: 1,
                iteration: p1_itr,
                pv_row: usize::MAX,
                pv_col: usize::MAX,
                base_vars: base_vars.clone(),
                non_base_vars: non_base_vars.clone(),
                matrix: mat.clone(),
            })
        }

        if mat[(number_of_eq, n)] != T::zero() {
            return Err("Problem is unbound (1)");
        }

        let help_var = number_of_vars - self.a_eq.layout().rows() + 1;
        if non_base_vars.iter().fold(0usize, |x, &y| usize::max(x, y)) >= help_var {
            return Err("Problem is unbound (2)");
        }

        // Phase 2

        // SAFTY:
        // Will be filled up in the next steps
        let layout = MatrixLayout::new(number_of_eq + 1, n + 1 - self.a_eq.layout().rows());
        let mut mx = Matrix::<T>::zeroed(layout);

        // Copy releveant Colums
        let mut col_idx = 0;
        for col in 0..mat.layout().cols() {
            if (col != mat.layout().cols() - 1)
                && (base_vars[col] >= (number_of_vars - self.a_eq.layout().rows()))
            {
                continue;
            }

            // Copy colum mat(0..<number_eq, col) to (0..., col_idx)
            for i in 0..number_of_eq {
                mx[(i, col_idx)] = mat[(i, col)];
            }
            // Copy mat(number_eq, col_idx)
            mx[(number_of_eq, col_idx)] = mat[(number_of_eq + 1, col)];

            col_idx += 1;
        }

        base_vars = base_vars
            .iter()
            .filter(|&&x| x < (number_of_vars - self.a_eq.layout().rows()))
            .copied()
            .collect();

        // let p2_hash_list = Vec::new()
        let mut p2_itr = 0;

        loop {
            // Get pivot

            let mut col_select = T::infinity();
            let mut pv_col = None;

            for i in 0..(mx.layout().cols() - 1) {
                if mx[(number_of_eq, i)] <= T::zero() {
                    continue;
                }
                if mx[(number_of_eq, i)] < col_select {
                    col_select = mx[(number_of_eq, i)];
                    pv_col = Some(i);
                }
            }

            if pv_col.is_none() {
                break;
            }
            let pv_col = pv_col.unwrap();

            let mut row_select = T::infinity();
            let mut pv_row = None;

            for j in 0..number_of_eq {
                if mx[(j, pv_col)] == T::zero() {
                    continue;
                }
                let r_value = mx[(j, mx.layout().cols() - 1)] / mx[(j, pv_col)];
                if r_value < T::zero() {
                    continue;
                }
                if r_value < row_select {
                    row_select = r_value;
                    pv_row = Some(j);
                }
            }

            if pv_row.is_none() {
                break;
            }
            let pv_row = pv_row.unwrap();

            if p2_itr >= options.max_p2_iterations {
                return Err("Exceeded maximum phase two iteration count");
            }
            p2_itr += 1;

            // TODO: Hash loop detection

            // Log increment
            if options.verbose {
                increments.push(LOPIncrement {
                    phase: 2,
                    iteration: p2_itr,
                    pv_row,
                    pv_col,
                    base_vars: base_vars.clone(),
                    non_base_vars: non_base_vars.clone(),
                    matrix: mx.clone(),
                })
            }

            let pv = mx[(pv_row, pv_col)];

            // Index swap
            swap(&mut base_vars[pv_col], &mut non_base_vars[pv_row]);

            let old = mx.clone();

            for row in 0..mx.layout().rows() {
                for col in 0..mx.layout().cols() {
                    match (row == pv_row, col == pv_col) {
                        (true, true) => mx[(row, col)] = T::one() / pv,
                        (true, false) => mx[(row, col)] = old[(row, col)] / pv,
                        (false, true) => mx[(row, col)] = (old[(row, col)] / pv).neg(),
                        (false, false) => {
                            mx[(row, col)] =
                                old[(row, col)] - (old[(pv_row, col)] * old[(row, pv_col)]) / pv;
                        }
                    }
                }
            }
        }

        if options.verbose {
            increments.push(LOPIncrement {
                phase: 2,
                iteration: p2_itr,
                pv_row: usize::MAX,
                pv_col: usize::MAX,
                base_vars: base_vars.clone(),
                non_base_vars: non_base_vars.clone(),
                matrix: mx.clone(),
            })
        }

        let mut solution = Matrix::zeroed(MatrixLayout::new(n, 1));
        for i in 0..non_base_vars.len() {
            let k = non_base_vars[i];
            if k > self.c.layout().rows() {
                continue;
            }
            solution[k - 1] = mx[(i, mx.layout().cols() - 1)];
        }

        Ok(LOPSolution {
            x: solution,
            fval: mx[(mx.layout().rows() - 1, mx.layout().cols() - 1)] + self.l,
            increments,
        })
    }
}

///
/// A set of options to define the solving process of a linear optimization problem.
///
/// # Example
///
/// ```no_run
/// use linalg::prelude::*;
///
/// # fn get_some_lop() -> LOP<f64> { todo!() }
/// let lop = get_some_lop();
/// let opts = LOPOptions {
///     max_p1_iterations: 8,
///     max_p2_iterations: 8,
///     verbose: true,
/// };
///
/// let sol = lop.solve_with(opts).unwrap();
/// ```
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LOPOptions {
    /// A limit to the number of iterations for the first normalization iteration (reducing Ax=b).
    pub max_p1_iterations: usize,
    /// A limit to the number of iterations for the final solving iteration (solving Aeqx=beq).
    pub max_p2_iterations: usize,
    /// A flag to enable logging of increments.
    pub verbose: bool,
}

impl Default for LOPOptions {
    fn default() -> Self {
        Self {
            max_p1_iterations: 8,
            max_p2_iterations: 8,
            verbose: false,
        }
    }
}

///
/// A step in the simplex algorithm saved.
///
/// # Example Output
///
/// ```txt
/// N:2 |       x7       x8       x3       x4       x5 |       
/// ----------------------------------------------------------------
///  x6 |      -1       -1        1        1        1  |      2
///  x1 |       1        1       -2       -1       -1  |      2
///  x2 |       0        1       -2        0       -1  |      1
/// (H) |      -1       -1        0        0        0  |      0
/// (P) |       1       -2        2       -1        2  |     -1
/// ```
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LOPIncrement<T: Num> {
    /// The phase the increment was recorded (phase 1 - normilization / phase 2 - solving).
    pub phase: usize,
    /// The exlicit number of the iteration (should be implicitly known from vector position).
    pub iteration: usize,
    /// The row index of the current pivot.
    pub pv_row: usize,
    /// The collum index of the current pivot.
    pub pv_col: usize,
    /// Non-Zero vars in simplex
    pub base_vars: Vec<usize>,
    /// Zero vars in simplex
    pub non_base_vars: Vec<usize>,
    /// The current simplex tablau.
    pub matrix: Matrix<T>,
}

impl<T: Num> Display for LOPIncrement<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // MAIN

        // A Table cell will have 7 characters with 1 character padding
        // meaning 8 characters per cell
        // exepte the last line that cannot be pivoted so - 2
        // -> self.matrix.size() * 8 - 2
        let line_size = self.matrix.layout().cols() * 9;

        // addionaly ther will be tow vertical splitter before and after the main part a 2 chars
        // -> self.matrix.layout().rows() * 2 * 2
        let line_size = line_size + 2 * 2;

        // there will be a sign prefix (5 char) and a linebreak (1 char)
        // -> self.matrix.layout().rows() * 6
        let line_size = line_size + 6;

        // Line 1

        let phase_label = format!(
            "{:>4} ",
            if self.phase == 1 {
                format!("N:{}", self.iteration)
            } else {
                format!("S:{}", self.iteration)
            }
        );
        write!(f, "{}", phase_label)?;
        write!(f, "| ")?;
        for bv in &self.base_vars {
            let s = format!("{:>8} ", format!("x{}", bv));
            write!(f, "{}", s)?;
        }
        writeln!(f, "|       ")?;

        // Line 2

        for _ in 0..(line_size) {
            write!(f, "-")?
        }
        writeln!(f)?;

        // N lines
        for i in 0..self.matrix.layout().rows() {
            let label = match i >= self.non_base_vars.len() {
                true => match i - self.non_base_vars.len() {
                    0 => String::from("(H)"),
                    _ => String::from("(P)"),
                },
                false => format!("x{}", self.non_base_vars[i]),
            };

            let label = format!("{:>4} ", label);
            write!(f, "{}| ", label)?;

            // for m - 1 cells
            for j in 0..(self.matrix.layout().cols() - 1) {
                let element = if (i, j) == (self.pv_row, self.pv_col) {
                    format!("[{:>6}] ", self.matrix[(i, j)])
                } else {
                    format!(" {:>6}  ", self.matrix[(i, j)])
                };
                write!(f, "{}", element)?;
            }

            writeln!(
                f,
                "| {:>6}",
                self.matrix[(i, self.matrix.layout().cols() - 1)]
            )?;
        }

        Ok(())
    }
}

///
/// A solution to a linear optimization problem.
///
/// TODO
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LOPSolution<T: Num> {
    /// The optimal vector to minimize the linear function.
    pub x: Matrix<T>,
    /// The resulting global minimum under the given side conditions.
    pub fval: T,
    /// A vector of incremental steps of the simplex algorihtm (empty if non-verbose run).
    pub increments: Vec<LOPIncrement<T>>,
}

impl<T: Num> Display for LOPSolution<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "=== LOP Solution ===")?;

        if self.increments.is_empty() {
            "<> No logged incr <>".fmt(f)?
        } else {
            for x in &self.increments {
                x.fmt(f)?
            }
        }

        write!(f, "\n {} ==> {}", self.x, self.fval)
    }
}

#[cfg(test)]
mod tests {

    use crate::{lop::*, matrix::*};

    #[test]
    fn test_builder() {
        let explicit: LOP<f64> = LOP::new(
            Matrix::from(vec![1.0, -3.0, 2.0, 0.0, 0.0]),
            0.0,
            Matrix::new((1, 5), vec![1.0, 0.0, -1.0, 0.0, 0.0]),
            Matrix::from(vec![4.0]),
            Matrix::new(
                (2, 5),
                vec![1.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, -2.0, 0.0, -1.0],
            ),
            Matrix::from(vec![1.0, 1.0]),
        );

        let implict: LOP<f64> = LOP::build(
            vec![1.0, -3.0, 2.0, 0.0, 0.0],
            0.0,
            Matrix::new((1, 5), vec![1.0, 0.0, -1.0, 0.0, 0.0]),
            vec![4.0],
            Matrix::new(
                (2, 5),
                vec![1.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, -2.0, 0.0, -1.0],
            ),
            vec![1.0, 1.0],
        );

        assert_eq!(explicit, implict);
    }

    #[test]
    fn check_short_two_phase_success() {
        let lop = LOP::new(
            Matrix::colvec(vec![1.0, -3.0, 2.0, 0.0, 0.0]),
            0.0,
            Matrix::rowvec(vec![1.0, 0.0, -1.0, 0.0, 0.0]),
            Matrix::colvec(vec![4.0]),
            Matrix::new(
                (2, 5),
                vec![1.0, -1.0, 0.0, -1.0, 0.0, 0.0, 1.0, -2.0, 0.0, -1.0],
            ),
            Matrix::colvec(vec![1.0, 1.0]),
        );

        let x: Result<LOPSolution<f64>, &'static str> = lop.solve_with(LOPOptions {
            max_p1_iterations: 2,
            max_p2_iterations: 2,
            verbose: false,
        });

        assert!(x.is_ok());
        let x = x.unwrap();

        assert_eq!(x.fval, -5.0);
        assert_eq!(*x.x.layout(), MatrixLayout::new(5, 1));
        assert_eq!(*x.x.raw(), vec![6.0, 5.0, 2.0, 0.0, 0.0]);
        assert_eq!(x.increments, vec![]);
    }

    #[test]
    fn check_failure_unbound() {
        let lop = LOP::new(
            Matrix::colvec(vec![8.0, 8.0, -9.0, 0.0, 0.0]),
            0.0f64,
            Matrix::rowvec(vec![1.0, 1.0, 1.0, 0.0, 0.0]),
            Matrix::colvec(vec![1.0]),
            Matrix::new(
                MatrixLayout::new(2, 5),
                vec![2.0, 4.0, 1.0, -1.0, 0.0, 1.0, -1.0, -1.0, 0.0, -1.0],
            ),
            Matrix::colvec(vec![8.0, 2.0]),
        );

        let x = lop.solve();
        assert!(x.is_err());

        let err = x.unwrap_err();
        assert_eq!(err, "Problem is unbound (1)")
    }

    #[test]
    fn check_verbose_succes_example() {
        let lop = LOP::new(
            Matrix::colvec(vec![4.0, -2.0, -5.0, 0.0]),
            0.0f64,
            Matrix::new(
                MatrixLayout::new(2, 4),
                vec![-5.0, 2.0, 9.0, 0.0, -2.0, 1.0, 4.0, 0.0],
            ),
            Matrix::colvec(vec![2.0, 1.0]),
            Matrix::rowvec(vec![-13.0, 7.0, 27.0, -1.0]),
            Matrix::colvec(vec![3.0]),
        );

        let mut options = LOPOptions::default();
        options.verbose = true;

        let x = lop.solve_with(options);
        assert!(x.is_ok());

        let x = x.unwrap();
        assert_eq!(x.fval, -2.0);
        assert_eq!(*x.x.layout(), MatrixLayout::new(4, 1));
        assert_eq!(*x.x.raw(), vec![0.0, 1.0, 0.0, 4.0]);

        assert!(x.increments.len() == 5);
    }
}
