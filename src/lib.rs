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
}
