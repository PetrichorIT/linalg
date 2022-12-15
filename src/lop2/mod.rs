mod defs;
use std::fmt::Display;
use std::ops::Neg;

pub use defs::{LOPParameterBound, LOP};

use crate::matrix::Matrix;
use num_traits::Num;

pub struct Simplex<T: Num + Copy> {
    tableau: Matrix<T>,
}

impl<T: Num + Copy + PartialOrd + Neg<Output = T>> Simplex<T> {
    pub fn new(lop: LOP<T>) -> Simplex<T> {
        let lop = lop.into_standard_form();
        let LOP::StandardForm { c, a, b } = lop else {
            panic!("This should not have happened")
        };

        let mut tableau = Matrix::fill((a.layout().rows() + 1, a.layout().cols() + 1), T::zero());
        tableau.insert(0.., 0.., &a);
        tableau.insert(0.., a.layout().cols().., &b);
        tableau.insert(a.layout().rows().., 0.., &c.transposed());

        let x = Matrix::fill((a.layout().cols(), 1), T::zero());
        let z = Matrix::mmul(&c.transposed(), &x);
        assert_eq!(z.len(), 1);
        tableau[(a.layout().rows(), a.layout().cols())] = z[0];

        Self { tableau }
    }

    pub fn solve(mut self, max_steps: usize) -> Result<(), &'static str>
    where
        T: Display,
    {
        // Asumes step 2

        let n = self.tableau.layout().cols() - 1;
        let m = self.tableau.layout().rows() - 1;
        for _ in 0..max_steps {
            println!("{}", self.tableau);
            // Chosse j such that dj < 0
            // Since last rows is -dj search for a positive one
            let mut j = usize::MAX;
            for i in 0..n {
                if self.tableau[(m, i)] > T::zero() {
                    j = i;
                    break;
                }
            }

            if j == usize::MAX {
                return Ok(());
            }

            // Choose an l such that A[l, j] > 0 and
            let mut can = Vec::new();
            for i in 0..m {
                if self.tableau[(i, j)] > T::zero() {
                    can.push((i, self.tableau[(i, n)] / self.tableau[(i, j)]))
                }
            }

            let min = can
                .into_iter()
                .min_by(|lhs, rhs| lhs.1.partial_cmp(&rhs.1).unwrap());

            let Some((l, _)) = min else {
                return Err("Problem is unbounded");
            };

            // Pivot A[l, j]
            // (1) Make this row basr worth by making the pivot 1
            let factor = T::one() / self.tableau[(l, j)];
            for i in 0..=n {
                self.tableau[(l, i)] = self.tableau[(l, i)] * factor;
            }

            // (2) For all other rows, make zero
            for k in 0..m {
                if k == l {
                    continue;
                }
                let factor = -self.tableau[(k, j)];
                for i in 0..=n {
                    self.tableau[(k, i)] = self.tableau[(k, i)] + factor * self.tableau[(l, i)];
                }
            }

            let factor = -self.tableau[(m, j)];
            for i in 0..=n {
                self.tableau[(m, i)] = self.tableau[(m, i)] + factor * self.tableau[(l, i)];
            }
        }

        Ok(())
    }
}

impl<T: Num + Copy + Display> Display for Simplex<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.tableau.fmt(f)
    }
}

#[cfg(test)]
mod tests;
