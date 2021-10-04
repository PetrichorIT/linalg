pub mod lse;
pub mod matrix;

#[cfg(test)]
mod tests {

    use crate::lse::*;
    use crate::matrix::*;
    use std::convert::TryFrom;

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn lr() {
        let matrix: Matrix<f64> = Matrix::try_from(vec![
            vec![2.0, 4.0, -2.0],
            vec![4.0, 9.0, -3.0],
            vec![-2.0, -3.0, 7.0],
        ])
        .expect("Failed");

        let lr = LrDecomposition::create(matrix.clone()).unwrap();

        println!("{}", lr);

        let x = lr.solve(Matrix::from(vec![1.0, 2.0, 3.0]));

        println!("{}", x);

        let bb = matrix * x;

        println!("{}", bb);
    }
}
