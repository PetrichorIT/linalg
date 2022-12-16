mod defs;
use std::fmt::Display;
use std::ops::Neg;

pub use defs::LOP;

use crate::matrix::Matrix;
use num_traits::Num;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SimplexOptions {
    pub max_steps: usize,
}

pub fn simplex<T>(lop: LOP<T>, options: SimplexOptions) -> Result<Matrix<T>, &'static str>
where
    T: Num + Copy + PartialOrd + Neg<Output = T> + Display + std::fmt::Debug,
{
    let SimplexOptions { mut max_steps } = options;

    let lop = lop.into_standard_form();
    let LOP::StandardForm { c, a, b, .. } = lop else {
        panic!("LOP::into_standard_form has not returned a LOP in standard form")
    };

    let mut tableau = Matrix::fill((a.layout().rows() + 1, a.layout().cols() + 1), T::zero());
    tableau.insert(0.., 0.., &a);
    tableau.insert(0.., a.layout().cols().., &b);
    tableau.insert(a.layout().rows().., 0.., &c.transposed());

    let x = Matrix::fill((a.layout().cols(), 1), T::zero());

    let z = Matrix::mmul(&c.transposed(), &x);
    assert_eq!(z.len(), 1);
    tableau[(a.layout().rows(), a.layout().cols())] = z[0];

    let n = tableau.layout().cols() - 1;
    let m = tableau.layout().rows() - 1;

    // Phase 1
    let mut based = Vec::new();
    let mut non_based = Vec::new();
    for row in 0..m {
        let base_col = (0..n).find(|&col| {
            let col = tableau.extract(0..n, col..(col + 1));
            let mut r = true;
            for i in 0..m {
                if i == row {
                    r &= col[i] == T::one()
                } else {
                    r &= col[i] == T::zero()
                }
            }
            r
        });

        match base_col {
            Some(col) => based.push(col),
            None => non_based.push(row),
        };
    }

    println!("{:?} {:?}", based, non_based);

    let mut base_state = vec![true; n];
    for i in 0..(n - m) {
        base_state[i] = false
    }

    if !non_based.is_empty() {
        // Phase 1 must be caclulated

        let p1n = n + non_based.len();

        let mut p1 = Matrix::zeroed((tableau.layout().rows() + 1, p1n + 1));
        p1.insert(0.., 0.., &tableau);
        p1.insert(
            0..,
            (tableau.layout().cols() - 1)..,
            &Matrix::colvec(vec![T::zero(); m + 1]),
        );
        p1.insert(
            0..,
            p1n..,
            &tableau.extract(0.., (tableau.layout().cols() - 1)..),
        );

        // Insert h1...hn
        let offset = tableau.layout().cols() - 1;
        for i in 0..non_based.len() {
            p1[(non_based[i], offset + i)] = T::one();
        }

        // Insert z2
        for i in 0..non_based.len() {
            p1[(m + 1, offset + i)] = T::one();
        }
        let mut base_state_p1 = vec![false; p1n];
        for i in offset..p1n {
            base_state_p1[i] = true;
        }
        for col in based {
            base_state_p1[col] = true;
        }

        let z2 = p1.layout().rows() - 1;

        // Correct z2 so that base elements are zero
        for i in 0..non_based.len() {
            let row = non_based[i];
            // all elements are one
            for i in 0..=p1n {
                p1[(z2, i)] = p1[(z2, i)] - p1[(row, i)]
            }
        }

        println!(">{}", p1);

        match simplex_iteration(&mut p1, &mut base_state_p1, p1n, m, m + 1, &mut max_steps) {
            Ok(_) => {
                let z2 = p1[(m + 1, p1n)];
                if z2 != T::zero() {
                    return Err("soloutions space is empty");
                }
            }
            Err(str) => eprintln!("{str}"),
        }

        // Copy old tabelau data
        tableau.insert(0.., 0.., &p1);
        tableau.insert(0.., n.., &p1.extract(0.., p1n..));
        base_state_p1.truncate(n);
        base_state = base_state_p1;

        // todo!()
    }

    // Asumes step 2

    simplex_iteration(&mut tableau, &mut base_state, n, m, m, &mut max_steps)
}

fn simplex_iteration<T>(
    tableau: &mut Matrix<T>,
    base_state: &mut Vec<bool>,
    n: usize,
    m: usize,
    c: usize,
    itr: &mut usize,
) -> Result<Matrix<T>, &'static str>
where
    T: Num + Copy + PartialOrd + Neg<Output = T> + Display + std::fmt::Debug,
{
    loop {
        if *itr == 0 {
            return Err("iteration counter exceeded");
        }
        *itr -= 1;

        println!("{}", tableau);
        // Chosse j such that dj < 0
        // Since last rows is -dj search for a positive one
        let mut j = usize::MAX;
        for i in 0..n {
            if tableau[(c, i)] < T::zero() {
                j = i;
                break;
            }
        }

        if j == usize::MAX {
            let mut x = Matrix::zeroed((n, 1));
            for i in 0..n {
                if base_state[i] {
                    let row = (0..m)
                        .into_iter()
                        .find(|j| tableau[(*j, i)] == T::one())
                        .unwrap();
                    x[(i, 0)] = tableau[(row, n)];
                }
            }

            return Ok(x);
        }

        // Choose an l such that A[l, j] > 0 and
        let mut can = Vec::new();
        for i in 0..m {
            if tableau[(i, j)] > T::zero() {
                can.push((i, tableau[(i, n)] / tableau[(i, j)]))
            }
        }

        let min = can
            .into_iter()
            .min_by(|lhs, rhs| lhs.1.partial_cmp(&rhs.1).unwrap());

        let Some((l, _)) = min else {
            return Err("problem is unbounded");
        };

        // Find replace base col
        let (idx, _) = base_state
            .iter()
            .enumerate()
            .find(|(i, c)| **c && tableau[(l, *i)] == T::one())
            .unwrap();

        base_state[idx] = false;
        base_state[j] = true;

        println!("pivot = {l},{j}");

        // Pivot A[l, j]
        // (1) Make this row basr worth by making the pivot 1
        let factor = T::one() / tableau[(l, j)];
        for i in 0..=n {
            tableau[(l, i)] = tableau[(l, i)] * factor;
        }
        // Ensure that pivot is 1 (it should be either way but this is for Eq)
        tableau[(l, j)] = T::one();

        // (2) For all other rows, make zero
        for k in 0..tableau.layout().rows() {
            if k == l {
                continue;
            }
            let factor = -tableau[(k, j)];
            for i in 0..=n {
                tableau[(k, i)] = tableau[(k, i)] + factor * tableau[(l, i)];
            }
            tableau[(k, j)] = T::zero();
        }

        // let factor = -tableau[(m, j)];
        // for i in 0..=n {
        //     tableau[(m, i)] = tableau[(m, i)] + factor * tableau[(l, i)];
        // }
    }
}
#[cfg(test)]
mod tests;
