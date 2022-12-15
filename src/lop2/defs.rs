use std::{
    fmt::{Debug, Display},
    ops::{Deref, Neg},
};

use crate::matrix::Matrix;
use num_traits::Num;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LOP<T: Num + Copy> {
    /// An optimization problem defined by:
    ///  max c^T * x subject to Ax <= b and x >= 0
    ///  dim(c) = (n, 1)
    ///  dim(A) = (n, m)
    ///  dim(b) = (m, 1)
    ///
    /// Note that b has no explicit bounds
    CanonicalForm {
        c: Matrix<T>,
        a: Matrix<T>,
        b: Matrix<T>,
    },
    /// An optimization problem defined by:
    ///  max c^T * x subject to Ax = b and x >= 0 and b >= 0
    ///  dim(c) = (n, 1)
    ///  dim(A) = (n, m)
    ///  dim(b) = (m, 1)
    ///
    StandardForm {
        c: Matrix<T>,
        a: Matrix<T>,
        b: Matrix<T>,
    },
    /// An optimization problem defined by:
    ///  max c^T * x subject to Ax<=b and A_eq <= b
    Mixed {
        c: Matrix<T>,
        a: Matrix<T>,
        b: Matrix<T>,
        a_eq: Matrix<T>,
        b_eq: Matrix<T>,
        x_bounds: Vec<LOPParameterBound<T>>,
    },
}

impl<T: Num + Copy + PartialOrd + Neg<Output = T>> LOP<T> {
    fn check(&self) {
        match self {
            Self::CanonicalForm { c, a, b } => {
                let n = c.layout().rows();
                assert_eq!(c.layout().cols(), 1);
                let m = a.layout().rows();
                assert_eq!(a.layout().cols(), n);
                assert!(m <= n);

                assert_eq!(b.layout().rows(), m);
                assert_eq!(b.layout().cols(), 1);
            }
            Self::StandardForm { c, a, b } => {
                let n = c.layout().rows();
                assert_eq!(c.layout().cols(), 1);
                let m = a.layout().rows();
                assert_eq!(a.layout().cols(), n);
                assert!(m <= n);

                assert_eq!(b.layout().rows(), m);
                assert_eq!(b.layout().cols(), 1);
                assert!(b.iter().all(|v| *v >= T::zero()))
            }
            Self::Mixed { .. } => {}
        }
    }

    pub fn into_standard_form(&self) -> Self {
        self.check();

        match self {
            Self::CanonicalForm { c, a, b } => {
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
                for i in 0..m {
                    if rb[i] < T::zero() {
                        rb[i] = rb[i] * -T::one();
                        for j in 0..(n + m) {
                            ra[(i, j)] = ra[(i, j)] * -T::one();
                        }
                    }
                }

                Self::StandardForm {
                    c: rc,
                    a: ra,
                    b: rb,
                }
            }
            Self::StandardForm { .. } => self.clone(),
            Self::Mixed {
                c,
                a,
                b,
                a_eq,
                b_eq,
                x_bounds,
            } => {
                let n = a.layout().cols();
                let m = a.layout().rows();

                // base size + bounds vars + y
                let mm = a.layout().rows() + a_eq.layout().rows();
                let nn = n + m + x_bounds.iter().filter(|v| v.need_s()).count();

                let mut ra = Matrix::fill((mm, nn), T::zero());
                ra.insert(0.., 0.., &a);
                ra.insert(a.layout().rows().., 0.., &a_eq);

                let mut rb = Matrix::fill((mm, 1), T::zero());
                rb.insert(0.., 0.., b);
                rb.insert(b.len().., 0.., b_eq);

                let mut rc = Matrix::fill((nn, 1), T::zero());
                rc.insert(0.., 0.., &c);

                let mut ptr = n;
                for i in 0..n {
                    match x_bounds[i] {
                        LOPParameterBound::Unbounded => {
                            // Since x_i is unbounded
                            // let x_i be x_i(+) >= 0
                            // and x_ptr be x_i(-) >= 0
                            for r in 0..mm {
                                ra[(r, ptr)] = ra[(r, i)]
                            }
                            rc[ptr] = -rc[i];
                            ptr += 1;
                        }
                        LOPParameterBound::GreaterOrEqualThan(v) if v != T::zero() => {
                            unimplemented!()
                        }
                        LOPParameterBound::GreaterOrEqualThan(v) if v == T::zero() => { /* NOP */ }
                        LOPParameterBound::LesserOrEqualThan(v) if v != T::zero() => {
                            unimplemented!()
                        }
                        LOPParameterBound::LesserOrEqualThan(v) if v == T::zero() => {
                            for r in 0..mm {
                                ra[(r, i)] = -ra[(r, i)]
                            }
                            rc[i] = -rc[i]
                        }
                        _ => unreachable!(),
                    }
                }

                // assert_eq!(nn - ptr, mm);

                for i in 0..a.layout().rows() {
                    ra[(i, ptr + i)] = T::one();
                }

                for i in 0..mm {
                    if rb[i] < T::zero() {
                        rb[i] = rb[i] * -T::one();
                        for j in 0..nn {
                            ra[(i, j)] = ra[(i, j)] * -T::one();
                        }
                    }
                }

                Self::StandardForm {
                    c: rc,
                    a: ra,
                    b: rb,
                }
            }
        }
    }
}

impl<T: Num + Copy + Debug + Display> Display for LOP<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::CanonicalForm { c, a, b } => {
                writeln!(f, "max {:?} * x subject to", c.deref())?;
                writeln!(f, " Ax <= b with x >= 0")?;
                writeln!(f, " {}", a)?;
                writeln!(f, " {}", b)?;
            }
            Self::StandardForm { c, a, b } => {
                writeln!(f, "max {:?} * x subject to", c.deref())?;
                writeln!(f, " Ax = b with x,b >= 0")?;
                writeln!(f, " {}", a)?;
                writeln!(f, " {}", b)?;
            }
            Self::Mixed {
                c,
                a,
                b,
                a_eq,
                b_eq,
                x_bounds,
            } => {
                writeln!(f, "max {:?} * x subject to", c.deref())?;
                writeln!(f, " A*x <= b")?;
                writeln!(f, " {}", a)?;
                writeln!(f, " {}", b)?;
                writeln!(f, " A_eq*x = b_eq")?;
                writeln!(f, " {}", a_eq)?;
                writeln!(f, " {}", b_eq)?;
                writeln!(f, " with {:?}", x_bounds)?;
            }
        }

        Ok(())
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum LOPParameterBound<T: Num> {
    Unbounded,
    GreaterOrEqualThan(T),
    LesserOrEqualThan(T),
}

impl<T: Num> LOPParameterBound<T> {
    fn need_s(&self) -> bool {
        match self {
            Self::Unbounded => true,
            Self::GreaterOrEqualThan(v) if *v == T::zero() => false,
            Self::LesserOrEqualThan(v) if *v == T::zero() => false,
            _ => true,
        }
    }
}
