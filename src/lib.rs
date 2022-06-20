mod marco;
mod num;

pub mod algebra;
pub mod codes;
pub mod discret;
pub mod lop;
pub mod matrix;
pub mod poly;
pub mod prelude;
pub mod sets;
pub mod vector;

#[cfg(test)]
mod tests {

    use crate::algebra::inv;
    use crate::algebra::vandermonde;
    use crate::lop::LOPOptions;
    use crate::lop::LOP;
    use crate::mat;
    use crate::matrix;
    use crate::matrix::*;
    use crate::prelude::binom;
    use crate::prelude::eig;
    use crate::prelude::GF;
    use crate::vector::VectorSpace;

    #[test]
    fn vandermonde_test() {
        let m = vandermonde(&[10, 20, 30, 40, 50]);
        println!("{}", m)
    }

    #[test]
    fn vector_space() {
        let _e = VectorSpace::<GF<2>, 3>::unit_vectors();
        let _a = VectorSpace::<GF<2>, 3>::all_vectors();

        let _eq = VectorSpace::<GF<3>, 2>::equivalence_classes();
        println!("{:?}", _eq)
    }

    #[test]
    fn it_works() {
        let mut matrix = Matrix::zeroed((7, 7));
        for i in 0..7 {
            for j in 0..=i {
                matrix[(i, j)] = binom(i, j);
            }
        }

        println!("{}", matrix);

        let mt = Matrix::fill((4, 4), 9);

        matrix.insert(2..=3, 2.., &mt);

        println!("{}", matrix);
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

        println!("{:?}", eig(matrix))
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

    #[test]
    fn lua() {
        assert_eq!(
            crate::algebra::lua(&[matrix![1;0;0;], matrix![0;1;0;], matrix![0;0;1;]]),
            true
        );

        assert_eq!(
            crate::algebra::lua(&[matrix![1;0;0;], matrix![0;1;0;], matrix![-8;2;0;]]),
            false
        );
    }
}
