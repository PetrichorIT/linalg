use std::{
    fmt::{Debug, Display},
    ops::{Deref, Neg},
};

use crate::matrix::Matrix;
use num_traits::Num;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LOP<T: Num + Copy> {
    /// An optimization problem defined by:
    ///  min c^T * x subject to Ax <= b and x >= 0
    ///  dim(c) = (n, 1)
    ///  dim(A) = (n, m)
    ///  dim(b) = (m, 1)
    ///
    /// Note that b has no explicit bounds
    CanonicalForm {
        c: Matrix<T>,
        a: Matrix<T>,
        b: Matrix<T>,
        labels: Vec<String>,
    },
    /// An optimization problem defined by:
    ///  min c^T * x subject to Ax = b and x >= 0 and b >= 0
    ///  dim(c) = (n, 1)
    ///  dim(A) = (n, m)
    ///  dim(b) = (m, 1)
    ///
    StandardForm {
        c: Matrix<T>,
        a: Matrix<T>,
        b: Matrix<T>,
        labels: Vec<String>,
    },
    /// An optimization problem defined by:
    ///  min c^T * x subject to A_lt x <= b_lt and A_gt x >= b_gt and A_eq x <= b_eq
    Mixed {
        c: Matrix<T>,
        a_lt: Matrix<T>,
        b_lt: Matrix<T>,
        a_gt: Matrix<T>,
        b_gt: Matrix<T>,
        a_eq: Matrix<T>,
        b_eq: Matrix<T>,
        labels: Vec<String>,
        // x_bounds: Vec<LOPParameterBound<T>>,
    },
}

impl<T: Num + Copy + PartialOrd + Neg<Output = T>> LOP<T> {
    fn check(&self) {
        match self {
            Self::CanonicalForm { c, a, b, labels } => {
                let n = a.layout().cols();
                let m = a.layout().rows();

                assert_eq!(c.layout().cols(), 1);
                assert_eq!(c.layout().rows(), n);

                // assert!(m <= n);

                assert_eq!(b.layout().rows(), m);
                assert_eq!(b.layout().cols(), 1);

                assert!([0, n].contains(&labels.len()));
            }
            Self::StandardForm { c, a, b, labels } => {
                let n = c.layout().rows();
                assert_eq!(c.layout().cols(), 1);
                let m = a.layout().rows();
                assert_eq!(a.layout().cols(), n);
                assert!(m <= n);

                assert_eq!(b.layout().rows(), m);
                assert_eq!(b.layout().cols(), 1);
                assert!(b.iter().all(|v| *v >= T::zero()));

                assert!([0, n].contains(&labels.len()))
            }
            Self::Mixed {
                c,
                a_lt,
                b_lt,
                a_gt,
                b_gt,
                a_eq,
                b_eq,
                labels,
            } => {
                let num_vars = a_eq
                    .layout()
                    .cols()
                    .max(a_lt.layout().cols())
                    .max(a_gt.layout().cols());

                assert!([num_vars, 0].contains(&a_eq.layout().cols()));
                assert!([num_vars, 0].contains(&a_lt.layout().cols()));
                assert!([num_vars, 0].contains(&a_gt.layout().cols()));

                assert_eq!(num_vars, c.layout().rows());
                assert_eq!(1, c.layout().cols());

                assert_eq!(a_eq.layout().rows(), b_eq.layout().rows());
                assert!([0, 1].contains(&b_eq.layout().cols()));

                assert_eq!(a_lt.layout().rows(), b_lt.layout().rows());
                assert!([0, 1].contains(&b_lt.layout().cols()));

                assert_eq!(a_gt.layout().rows(), b_gt.layout().rows());
                assert!([0, 1].contains(&b_gt.layout().cols()));

                assert!([0, num_vars].contains(&labels.len()))
            }
        }
    }

    pub fn into_standard_form(&self) -> Self {
        self.check();

        match self {
            Self::CanonicalForm { c, a, b, labels } => {
                // Firstly, add y
                let m = a.layout().rows();
                let n = a.layout().cols();
                let mut ra = Matrix::fill((m, n + m), T::zero());
                let mut rc = Matrix::fill((m + n, 1), T::zero());

                ra.insert(0.., 0.., &a);
                for i in 0..m {
                    ra[(i, n + i)] = T::one();
                }
                rc.insert(0.., 0.., c);

                let mut rb = b.clone();
                let mut labels = labels.clone();
                if labels.is_empty() {
                    labels = dummy_labels(n)
                }

                for i in 0..m {
                    if rb[i] < T::zero() {
                        rb[i] = rb[i] * -T::one();
                        for j in 0..(n + m) {
                            ra[(i, j)] = ra[(i, j)] * -T::one();
                        }
                    }
                    labels.push(format!("s{}", i));
                }

                Self::StandardForm {
                    c: rc,
                    a: ra,
                    b: rb,
                    labels,
                }
            }
            Self::StandardForm { .. } => self.clone(),
            Self::Mixed {
                c,
                a_lt,
                b_lt,
                a_gt,
                b_gt,
                a_eq,
                b_eq,
                labels,
            } => {
                let num_vars = a_eq
                    .layout()
                    .cols()
                    .max(a_lt.layout().cols())
                    .max(a_gt.layout().cols());
                let n = num_vars + a_lt.layout().rows() + a_gt.layout().rows();
                let m = a_eq.layout().rows() + a_lt.layout().rows() + a_gt.layout().rows();

                let mut a = Matrix::zeroed((m, n));
                let mut b = Matrix::zeroed((m, 1));

                let mut labels = labels.clone();
                if labels.is_empty() {
                    labels = dummy_labels(num_vars);
                }

                // Add a_eq to a
                for i in 0..a_eq.layout().rows() {
                    if b_eq[(i, 0)] < T::zero() {
                        let mut line = a_eq.extract(i.., 0..);
                        for k in 0..line.len() {
                            line[k] = -line[k];
                        }
                        a.insert(i.., 0.., &line);
                        b[(i, 0)] = -b_eq[(i, 0)];
                    } else {
                        a.insert(i.., 0.., &a_eq.extract(i.., 0..));
                        b[(i, 0)] = b_eq[(i, 0)];
                    }
                }

                let mut e_ptr = num_vars;
                let mut ptr = a_eq.layout().rows();

                // Fetch LT from lt
                for i in 0..a_lt.layout().rows() {
                    if b_lt[(i, 0)] < T::zero() {
                        continue;
                    }

                    a.insert(ptr.., 0.., &a_lt.extract(i.., 0..));
                    a[(ptr, e_ptr)] = T::one();
                    b[(ptr, 0)] = b_lt[(i, 0)];
                    labels.push(format!("s{}", e_ptr - num_vars));
                    ptr += 1;
                    e_ptr += 1;
                }

                // FETCH LT from gt
                for i in 0..a_gt.layout().rows() {
                    if b_gt[(i, 0)] >= T::zero() {
                        continue;
                    }
                    let mut line = a_gt.extract(i.., 0..);
                    for k in 0..line.len() {
                        line[k] = -line[k];
                    }

                    a.insert(ptr.., 0.., &line);
                    a[(ptr, e_ptr)] = T::one();
                    b[(ptr, 0)] = -b_gt[(i, 0)];
                    labels.push(format!("s{}", e_ptr - num_vars));
                    ptr += 1;
                    e_ptr += 1;
                }

                // FETCH GT from gt
                for i in 0..a_gt.layout().rows() {
                    if b_gt[(i, 0)] < T::zero() {
                        continue;
                    }

                    a.insert(ptr.., 0.., &a_gt.extract(i.., 0..));
                    a[(ptr, e_ptr)] = -T::one();
                    b[(ptr, 0)] = b_gt[(i, 0)];
                    labels.push(format!("r{}", e_ptr - num_vars));
                    ptr += 1;
                    e_ptr += 1;
                }

                // FETCH GT from lt
                for i in 0..a_lt.layout().rows() {
                    if b_lt[(i, 0)] >= T::zero() {
                        continue;
                    }
                    let mut line = a_lt.extract(i.., 0..);
                    for k in 0..line.len() {
                        line[k] = -line[k];
                    }

                    a.insert(ptr.., 0.., &line);
                    a[(ptr, e_ptr)] = -T::one();
                    b[(ptr, 0)] = -b_lt[(i, 0)];
                    labels.push(format!("r{}", e_ptr - num_vars));
                    ptr += 1;
                    e_ptr += 1;
                }

                assert_eq!(e_ptr, n);
                assert_eq!(ptr, m);

                let mut new_c = Matrix::zeroed((n, 1));
                new_c.insert(0.., 0.., &c);

                Self::StandardForm {
                    c: new_c,
                    a,
                    b,
                    labels,
                }
            } // Self::Mixed {
              //     c,
              //     a_lt,
              //     b_lt,
              //     a_gt,
              //     b_gt,
              //     a_eq,
              //     b_eq,
              // } => {
              //     let n = a.layout().cols();
              //     let m = a.layout().rows();

              //     // base size + bounds vars + y
              //     let mm = a.layout().rows() + a_eq.layout().rows();
              //     let nn = n + m + x_bounds.iter().filter(|v| v.need_s()).count();

              //     let mut ra = Matrix::fill((mm, nn), T::zero());
              //     ra.insert(0.., 0.., &a);
              //     ra.insert(a.layout().rows().., 0.., &a_eq);

              //     let mut rb = Matrix::fill((mm, 1), T::zero());
              //     rb.insert(0.., 0.., b);
              //     rb.insert(b.len().., 0.., b_eq);

              //     let mut rc = Matrix::fill((nn, 1), T::zero());
              //     rc.insert(0.., 0.., &c);

              //     let mut ptr = n;
              //     for i in 0..n {
              //         match x_bounds[i] {
              //             LOPParameterBound::Unbounded => {
              //                 // Since x_i is unbounded
              //                 // let x_i be x_i(+) >= 0
              //                 // and x_ptr be x_i(-) >= 0
              //                 for r in 0..mm {
              //                     ra[(r, ptr)] = ra[(r, i)]
              //                 }
              //                 rc[ptr] = -rc[i];
              //                 ptr += 1;
              //             }
              //             LOPParameterBound::GreaterOrEqualThan(v) if v != T::zero() => {
              //                 unimplemented!()
              //             }
              //             LOPParameterBound::GreaterOrEqualThan(v) if v == T::zero() => { /* NOP */ }
              //             LOPParameterBound::LesserOrEqualThan(v) if v != T::zero() => {
              //                 unimplemented!()
              //             }
              //             LOPParameterBound::LesserOrEqualThan(v) if v == T::zero() => {
              //                 for r in 0..mm {
              //                     ra[(r, i)] = -ra[(r, i)]
              //                 }
              //                 rc[i] = -rc[i]
              //             }
              //             _ => unreachable!(),
              //         }
              //     }

              //     // assert_eq!(nn - ptr, mm);

              //     for i in 0..a.layout().rows() {
              //         ra[(i, ptr + i)] = T::one();
              //     }

              //     for i in 0..mm {
              //         if rb[i] < T::zero() {
              //             rb[i] = rb[i] * -T::one();
              //             for j in 0..nn {
              //                 ra[(i, j)] = ra[(i, j)] * -T::one();
              //             }
              //         }
              //     }

              //     Self::StandardForm {
              //         c: rc,
              //         a: ra,
              //         b: rb,
              //     }
              // }
        }
    }
}

impl<T: Num + Copy + Debug + Display> Display for LOP<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CanonicalForm { c, a, b, labels } => {
                writeln!(f, "min {:?} * x subject to", c.deref())?;
                writeln!(f, " Ax <= b with x >= 0")?;
                writeln!(f, " {}", a)?;
                writeln!(f, " {}", b.transposed())?;
                writeln!(f, " labels: {:?}", labels)?;
            }
            Self::StandardForm { c, a, b, labels } => {
                writeln!(f, "min {:?} * x subject to", c.deref())?;
                writeln!(f, " Ax = b with x,b >= 0")?;
                writeln!(f, " {}", a)?;
                writeln!(f, " {}", b.transposed())?;
                writeln!(f, " labels: {:?}", labels)?;
            }
            Self::Mixed {
                c,
                a_lt,
                b_lt,
                a_gt,
                b_gt,
                a_eq,
                b_eq,
                labels,
            } => {
                writeln!(f, "min {:?} * x subject to", c.deref())?;
                writeln!(f, " A_lt*x <= b_lt")?;
                writeln!(f, " {}", a_lt)?;
                writeln!(f, " {}", b_lt.transposed())?;
                writeln!(f, " A_gt*x >= b_gt")?;
                writeln!(f, " {}", a_gt)?;
                writeln!(f, " {}", b_gt.transposed())?;
                writeln!(f, " A_eq*x = b_eq")?;
                writeln!(f, " {}", a_eq)?;
                writeln!(f, " {}", b_eq.transposed())?;
                writeln!(f, " labels: {:?}", labels)?;
            }
        }

        Ok(())
    }
}

// #[derive(Debug, Clone, PartialEq, Eq)]
// pub enum LOPParameterBound<T: Num> {
//     Unbounded,
//     GreaterOrEqualThan(T),
//     LesserOrEqualThan(T),
// }

// impl<T: Num> LOPParameterBound<T> {
//     fn need_s(&self) -> bool {
//         match self {
//             Self::Unbounded => true,
//             Self::GreaterOrEqualThan(v) if *v == T::zero() => false,
//             Self::LesserOrEqualThan(v) if *v == T::zero() => false,
//             _ => true,
//         }
//     }
// }

fn dummy_labels(n: usize) -> Vec<String> {
    (0..n).map(|i| format!("x{}", i)).collect()
}
