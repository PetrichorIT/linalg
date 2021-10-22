mod marco;
mod num;

pub mod core;
pub mod discret;
pub mod lop;
pub mod lse;
pub mod poly;
pub mod prelude;

#[allow(non_camel_case_types)]
pub type c32 = num_complex::Complex32;
#[allow(non_camel_case_types)]
pub type c64 = num_complex::Complex64;

#[cfg(test)]
mod tests {

    use crate::core::*;
    use crate::lop::LOPOptions;
    use crate::lop::LOP;
    use crate::lse::inv;
    use crate::lse::tridiag;
    use crate::mat;
    use crate::prelude::eigv;
    use crate::prelude::stirling2;

    #[test]
    fn it_works() {
        let mut m = Matrix::zeroed((5, 5));
        for i in 0..5 {
            for j in 0..=i {
                m[(i, j)] = stirling2(i, j)
            }
        }

        println!("{}", m)
    }

    #[test]
    fn lr_test() {
        let matrix = mat![
            2.0, 1.0;
            6.0, 4.0;
        ];

        let gj = inv(matrix);

        println!("{}", gj.unwrap());
    }

    #[test]
    fn eig_test() {
        let matrix = mat![
             2.0, -1.0;
             -4.0, 2.0;
        ];

        println!("{:?}", eigv(matrix, 100))
    }

    #[test]
    fn tdg_test() {
        let m = mat![
            1.0, 2.0, 6.0, -2.0;
            3.0, 1.0, 0.0, 65.0;
            (-2.0), -9.0, -0.4, 1.0;
            1.0, 1.0, 1.0, 1.0;
        ];

        let t = tridiag(m);
        println!("{}", t);
    }

    #[test]
    fn lop() {
        let lop = LOP::new(
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

        let x = lop.solve_with(LOPOptions {
            max_p1_iterations: 2,
            max_p2_iterations: 2,
            verbose: true,
        });

        println!("{}", x.unwrap());
    }
}
