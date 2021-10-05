use std::{mem::swap, ops::Neg};

use crate::matrix::{Matrix, MatrixLayout};

// min c^T*x + l with a_eq*x=b_eq and a*x <= b

/// A linear optimization problem type alias for [LinearOpimizationProblem].
pub type LOP<T> = LinearOpimizationProblem<T>;

///
/// A linear optimization problem in normal form.
///
/// This structure describes a minimation problem of type
/// $ \text{min} c^Tx + l with A_{eq}x=b_{eq} and Ax <= b and x >= 0 $
///
/// This structure also holds the function [LinearOpimizationProblem::solve] to
/// solve the optimization problem an return the caclulated result.
///
#[derive(Debug, Clone)]
pub struct LinearOpimizationProblem<T> {
    pub c: Matrix<T>,
    pub l: T,
    pub a: Matrix<T>,
    pub b: Matrix<T>,
    pub a_eq: Matrix<T>,
    pub b_eq: Matrix<T>,
}

impl<T> LinearOpimizationProblem<T> {
    /// Creates a new linear optimization problem.
    pub fn new(
        c: Matrix<T>,
        l: T,
        a: Matrix<T>,
        b: Matrix<T>,
        a_eq: Matrix<T>,
        b_eq: Matrix<T>,
    ) -> Self {
        Self {
            c,
            l,
            a,
            b,
            a_eq,
            b_eq,
        }
    }
}

impl LinearOpimizationProblem<f64> {
    ///
    /// Tries to solve the linear optimization problem using the provided
    /// options and observers.
    ///
    pub fn solve(&self) -> Result<LOPSolution<f64>, &'static str> {
        self.solve_with(LOPOptions::default())
    }

    ///
    /// Tries to solve the linear optimization problem using the provided
    /// options and observers.
    ///
    pub fn solve_with(&self, options: LOPOptions) -> Result<LOPSolution<f64>, &'static str> {
        assert!(self.c.layout().is_colvec());
        assert!(self.a.layout().rows() == self.b.layout().rows());
        assert!(self.a.layout().cols() == self.c.size() || self.a.layout().cols() == 0);
        assert!(self.a_eq.layout().rows() == self.b_eq.layout().rows());
        assert!(self.a_eq.layout().cols() == self.c.size() || self.a_eq.layout().cols() == 0);

        let n = self.c.layout().rows();
        let number_of_vars = n + self.a.layout().rows() + self.a_eq.layout().rows();
        let number_of_eq = self.a.layout().rows() + self.a_eq.layout().rows();

        let mut increments = Vec::new();

        // SAFTY:
        // All positions in the given matrix will be filled in the following
        // steps without using the matrix values as getter
        let layout = MatrixLayout::new(number_of_eq + 2, n + 1);
        let mut mat = unsafe { Matrix::<f64>::uninitalized(layout) };

        // Fill the Axy=b matrix
        if self.a.layout().rows() != 0 {
            // Place b at anchor (0, n)
            for i in 0..self.b.size() {
                mat[(i, n)] = self.b[i];
            }
            // Place a at anchor (0, 0)
            for i in 0..self.a.layout().rows() {
                for j in 0..self.a.layout().cols() {
                    mat[(i, j)] = self.a[(i, j)];
                }
            }
        }

        // Fill A_eqx=b_eq
        if self.a_eq.layout().rows() != 0 {
            // place b_eq at anchor (b::rows, n)
            for i in 0..self.b_eq.size() {
                mat[(self.b.size() + i, n)] = self.b_eq[i];
            }
            // place a_eq at anchor (a.rows(), 0)
            for i in 0..self.a_eq.layout().rows() {
                for j in 0..self.a_eq.layout().cols() {
                    mat[(self.a.layout().rows() + i, j)] = self.a_eq[(i, j)];
                }
            }
        }

        // Fill c matrix
        // place c transposed * -1 at anchor (number_of_eq + 1,0)
        for j in 0..self.c.size() {
            mat[(number_of_eq + 1, j)] = self.c[j].neg();
        }

        // SAFTY:
        // Concerning matrix only the sum fields are still uninitialized

        // Sum up all rows of a_eq
        let mut ch = Matrix::<f64>::fill(MatrixLayout::new(n + 1, 1), 0.0);
        for o in self.a.layout().rows()..number_of_eq {
            for j in 0..mat.layout().cols() {
                ch[j] += mat[(o, j)];
            }
        }

        // Fill in ch (transposed) into matrix at anchor (number_of_eq, 0)
        for j in 0..ch.size() {
            mat[(number_of_eq, j)] = ch[j];
        }

        // Prepare base va counters

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

            let mut col_select = f64::INFINITY;
            let mut pv_col = None;

            for i in 0..n {
                let raw_index = mat.layout().index((number_of_eq, i));
                if mat[raw_index] <= 0.0 {
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

            let mut row_select = f64::INFINITY;
            let mut pv_row = None;

            for j in 0..number_of_eq {
                if mat[(j, pv_col)] == 0.0 {
                    continue;
                }

                let r_value = mat[(j, n)] / mat[(j, pv_col)];
                if r_value < 0.0 {
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

            if p1_itr > options.max_p1_iterations {
                return Err(&"Exceeded maximum phase one iteration count");
            }
            p1_itr += 1;

            // TODO: Hasher to prevent endless loops

            // Log current iteration
            if options.verbose {
                increments.push(LOPIncrement {
                    phase: 1,
                    iteration: p1_itr,
                    pv_row: pv_row,
                    pv_col: pv_col,
                    matrix: mat.clone(),
                })
            }

            let pv = mat[(pv_row, pv_col)];

            // Swap indices
            swap(&mut base_vars[pv_col], &mut non_base_vars[pv_col]);

            let old = mat.clone();
            for row in 0..mat.layout().rows() {
                for col in 0..mat.layout().cols() {
                    match (row == pv_row, col == pv_col) {
                        (true, true) => mat[(row, col)] = 1.0 / pv,
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
                matrix: mat.clone(),
            })
        }

        if mat[(number_of_eq, n)] != 0.0 {
            return Err(&"Problem is unbound (1)");
        }

        let help_var = number_of_vars - self.a_eq.layout().rows() + 1;
        if non_base_vars.iter().fold(0usize, |x, &y| usize::max(x, y)) >= help_var {
            return Err(&"Problem is unbound (2)");
        }

        // Phase 2

        // SAFTY:
        // Will be filled up in the next steps
        let layout = MatrixLayout::new(number_of_eq + 1, n + 1 - self.a_eq.layout().rows());
        let mut mx = unsafe { Matrix::<f64>::uninitalized(layout) };

        // Copy releveant Colums
        let mut col_idx = 0;
        for col in 0..mat.layout().cols() {
            if (col != mat.layout().cols() - 1)
                && (base_vars[col] < (number_of_vars - self.a_eq.layout().rows()))
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

            let mut col_select = f64::INFINITY;
            let mut pv_col = None;

            for i in 0..(mx.layout().cols() - 1) {
                if mx[(number_of_eq, i)] <= 0.0 {
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

            let mut row_select = f64::INFINITY;
            let mut pv_row = None;

            for j in 0..number_of_eq {
                if mx[(j, pv_col)] == 0.0 {
                    continue;
                }
                let r_value = mx[(j, mx.layout().cols() - 1)] / mx[(j, pv_col)];
                if r_value < 0.0 {
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

            if p2_itr > options.max_p2_iterations {
                return Err(&"Exceeded maximum phase two iteration count");
            }
            p2_itr += 1;

            // TODO: Hash loop detection

            // Log increment
            if options.verbose {
                increments.push(LOPIncrement {
                    phase: 2,
                    iteration: p2_itr,
                    pv_row: pv_row,
                    pv_col: pv_col,
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
                        (true, true) => mx[(row, col)] = 1.0 / pv,
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
                matrix: mx.clone(),
            })
        }

        let mut solution = Matrix::fill(MatrixLayout::new(n, 1), 0.0);
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
            increments: increments,
        })
    }
}

///
/// A set of options to define the solving process of a linear optimization problem.
///
/// # Example
///
/// ```no_run
/// use linalg::{lop::*, matrix::*};
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
/// TODO
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LOPIncrement<T> {
    /// The phase the increment was recorded (phase 1 - normilization / phase 2 - solving).
    pub phase: usize,
    /// The exlicit number of the iteration (should be implicitly known from vector position).
    pub iteration: usize,
    /// The row index of the current pivot.
    pub pv_row: usize,
    /// The collum index of the current pivot.
    pub pv_col: usize,
    /// The current simplex tablau.
    pub matrix: Matrix<T>,
}

///
/// A solution to a linear optimization problem.
///
/// TODO
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct LOPSolution<T> {
    /// The optimal vector to minimize the linear function.
    pub x: Matrix<T>,
    /// The resulting global minimum under the given side conditions.
    pub fval: T,
    /// A vector of incremental steps of the simplex algorihtm (empty if non-verbose run).
    pub increments: Vec<LOPIncrement<T>>,
}
