mod marco;

pub mod core;
pub mod lop;
pub mod lse;

#[cfg(test)]
mod tests {

    use crate::core::*;
    use crate::lop::LOPOptions;
    use crate::lop::LOP;
    use crate::lse::*;
    use crate::mat;
    use std::convert::TryFrom;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);

        let v = mat![1, 2, 3; 1, 2, 3;];

        println!("{:?}", v)
    }

    #[test]
    fn lr() {
        let matrix: Matrix<f64> = Matrix::try_from(vec![
            vec![2.0, 4.0, -4.0],
            vec![1.0, 1.0, 2.0],
            vec![2.0, -3.0, 0.0],
        ])
        .expect("Failed");

        let qr = QrDecomposition::create(matrix.clone());

        println!("{}", qr);

        let b = Matrix::<f64>::from(vec![1.0, 2.0, 3.0]);

        let x = qr.solve(b);

        let bb = matrix * x;

        println!("{}", bb);
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
