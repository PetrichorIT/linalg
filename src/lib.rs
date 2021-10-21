mod marco;
mod num;

pub mod core;
pub mod lop;
pub mod lse;
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
    use crate::prelude::eig;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);

        let m = mat![ 1, 2, 0; 1, 2, 0;];

        println!("{}", m.is_lower_triag());
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

        println!("{:?}", eig(matrix, 100))
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
