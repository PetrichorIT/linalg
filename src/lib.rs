mod marco;
mod num;

pub mod algebra;
pub mod array;
pub mod codes;
pub mod discret;
pub mod ea;
pub mod lop;
pub mod lop2;
pub mod matrix;
pub mod poly;
pub mod prelude;
pub mod sets;
pub mod vector;

#[cfg(test)]
mod tests {

    use num_traits::Pow;

    use crate::algebra::gauss;
    use crate::algebra::inv;
    use crate::algebra::vandermonde;
    use crate::codes::generator_to_controll;
    use crate::codes::Code;
    use crate::codes::GenericSource;
    use crate::codes::LinearCode;
    use crate::codes::Word;
    use crate::lop::LOPOptions;
    use crate::lop::LOP;
    use crate::mat;
    use crate::matrix;
    use crate::matrix::*;
    use crate::prelude::binom;
    use crate::prelude::eig;
    use crate::prelude::GF;
    use crate::vector::VectorOrientation;
    use crate::vector::VectorSpace;

    #[test]
    fn p() {
        let mut s = 0;
        for k in 0..=4 {
            s += binom(7, k) as i32 * 2.pow(k)
        }
        println!("|B4(a)| = {}", s);
        println!("3^7 = {}", 3.pow(7u32));

        let code = LinearCode::<3, 7, 3>::with_control(mat![
            GF(1), GF(2), GF(0), GF(0);
            GF(0), GF(1), GF(2), GF(1);
            GF(1), GF(0), GF(1), GF(2);
            GF(2), GF(1), GF(0), GF(1);
            GF(0), GF(1), GF(2), GF(0);
            GF(0), GF(0), GF(1), GF(2);
            GF(2), GF(0), GF(0), GF(1);
        ]);

        println!("d(C) = {}", code.distance());

        let src = GenericSource::from_distr(vec![27, 9, 9, 9, 3, 3, 2, 2], |i| i);
        println!("H(Q) = {}", src.entropy(2.));

        let code = LinearCode::<2, 7, 3>::with_generator(mat![
            GF(1), GF(0), GF(0), GF(1), GF(0), GF(1), GF(1);
            GF(0), GF(1), GF(0), GF(1), GF(1), GF(0), GF(1);
            GF(0), GF(0), GF(1), GF(1), GF(1), GF(1), GF(0);
        ]);

        for w in code.all_codewords() {
            println!("{:?}", w)
        }

        for i in 1..=12 {
            println!("{}", 2.0f64.powi(i) as i32 % 11);
        }

        // for d in 1..=11 {
        //     // println!("2^{} = {}", i, 2.0.pow(i) as i32 % 11)
        //     for delta in 1..=11 {
        //         println!("{}, {}", d, delta);
        //         assert_ne!(0, (2.0f64.powi(2 * delta + 1) as i32 * d) % 11);
        //     }
        // }

        // for a in 3..=9 {
        //     for i in 1..=11 {
        //         assert_ne!(
        //             0,
        //             (2.0f64.powi(i) as i32 + a) % 11,
        //             "Failed at a = {}, i = {}",
        //             a,
        //             i
        //         )
        //     }
        // }

        for a in 0..11 {
            for c in 0..11 {
                if a == c {
                    continue;
                }

                for b in 0..11 {
                    let a: GF<11> = GF(a);
                    let b: GF<11> = GF(b);
                    let c: GF<11> = GF(c);

                    if a + GF(2) * b + GF(4) * c == c + GF(2) * b + GF(4) * a {
                        println!("a = {} b = {} c = {}", a, b, c)
                    }
                }
            }
        }
    }

    #[test]
    fn series_7_task_2() {
        println!("phi({}) => phi({})", GF::<11>(1), GF::<11>(1) / GF::<11>(1));
        println!("phi({}) => phi({})", GF::<11>(2), GF::<11>(1) / GF::<11>(2));
        println!("phi({}) => phi({})", GF::<11>(3), GF::<11>(1) / GF::<11>(3));
        println!("phi({}) => phi({})", GF::<11>(4), GF::<11>(1) / GF::<11>(4));
        println!("phi({}) => phi({})", GF::<11>(5), GF::<11>(1) / GF::<11>(5));
        println!("phi({}) => phi({})", GF::<11>(6), GF::<11>(1) / GF::<11>(6));
        println!("phi({}) => phi({})", GF::<11>(7), GF::<11>(1) / GF::<11>(7));
        println!("phi({}) => phi({})", GF::<11>(8), GF::<11>(1) / GF::<11>(8));
        println!("phi({}) => phi({})", GF::<11>(9), GF::<11>(1) / GF::<11>(9));
        println!(
            "phi({}) => phi({})",
            GF::<11>(10),
            GF::<11>(1) / GF::<11>(10)
        );

        println!("phi({}) => phi({})", GF::<10>(1), GF::<10>(1) / GF::<10>(1));
        println!("phi({}) => phi({})", GF::<10>(3), GF::<10>(1) / GF::<10>(3));
        println!("phi({}) => phi({})", GF::<10>(7), GF::<10>(1) / GF::<10>(7))
    }

    #[test]
    fn series_6_task_5() {
        let code = LinearCode::<3, 4, 2>::hamming().unwrap();
        let h = code.controll_matrix().unwrap();

        println!("{}", h);
        for c in code.all_codewords() {
            println!("- {:?}", c)
        }

        // 9 elements = 3^2
        // b1, b2 from

        let g = LinearCode::<3, 4, 2>::with_generator(mat![
            GF(1), GF(1), GF(1), GF(0);
            GF(1), GF(2), GF(0), GF(1);
        ]);

        println!("{}", code == g);

        'outer: for a in VectorSpace::<GF<3>>::new(4).all_vectors() {
            for b in VectorSpace::<GF<3>>::new(4).all_vectors() {
                let mut matrix = Matrix::fill((2, 4), GF(0));
                matrix.insert_vector(0.., 0.., a.clone().orientation(VectorOrientation::Row));
                matrix.insert_vector(1.., 0.., b.orientation(VectorOrientation::Row));

                let gencode = LinearCode::<3, 4, 2>::with_generator(matrix);
                if gencode == code {
                    println!("{}", gencode.generator_matrix().unwrap());
                    for c in gencode.all_codewords() {
                        println!("- {:?}", c)
                    }
                    break 'outer;
                }
            }
        }

        // println!("{}", code == gencode)
    }

    #[test]
    fn series_6_task_4() {
        let code = LinearCode::<2, 3, 1>::hamming().unwrap();
        println!("Ham_2(3, 1)");
        println!("{}", code.controll_matrix().unwrap());
        println!("|C| = {}", code.all_codewords().len());

        let code = LinearCode::<2, 15, 11>::hamming().unwrap();
        println!("Ham_2(15, 11)");
        println!("{}", code.controll_matrix().unwrap());
        println!("|C| = {}", code.all_codewords().len());

        let _code = LinearCode::<2, 30, 25>::hamming().is_none();

        let code = LinearCode::<3, 4, 2>::hamming().unwrap();
        println!("Ham_3(4, 2)");
        println!("{}", code.controll_matrix().unwrap());
        println!("|C| = {}", 3.0.pow(2.0));

        let code = LinearCode::<3, 13, 10>::hamming().unwrap();
        println!("Ham_3(13, 10)");
        println!("{}", code.controll_matrix().unwrap());
        println!("|C| = {}", 3.0.pow(3.0));
    }

    #[test]
    fn series_6_task_3() {
        let code = LinearCode::<2, 6, 0>::with_control(mat![
            GF(1), GF(1), GF(0);
            GF(1), GF(0), GF(1);
            GF(1), GF(0), GF(0);
            GF(1), GF(1), GF(1);
            GF(0), GF(1), GF(0);
            GF(1), GF(0), GF(0);
        ]);

        let all = code.all_codewords();
        let len = all.len();
        println!("|C| = {} ... dim(C) = {}", len, (len as f64).log2());
        println!("d(C) = {}", code.distance());

        let code = Code {
            words: all.into_iter().map(|v| Word(v.data)).collect::<Vec<_>>(),
        };

        println!(
            "{}",
            code.match_for(Word(Vec::from([0, 1, 1, 0, 0, 1].map(|v| GF(v)))))
                .0
        );
        println!(
            "{}",
            code.match_for(Word(Vec::from([0, 1, 1, 1, 1, 1].map(|v| GF(v)))))
                .0
        );
        println!(
            "{}",
            code.match_for(Word(Vec::from([1, 1, 1, 0, 1, 0].map(|v| GF(v)))))
                .0
        );
    }

    #[test]
    fn series_6_task_2() {
        let code = LinearCode::<2, 5, 3>::with_control(mat![
            GF(1), GF(1), GF(0);
            GF(1), GF(0), GF(1);
            GF(1), GF(1), GF(0);
            GF(1), GF(1), GF(0);
            GF(1), GF(1), GF(0);
        ]);

        for c in code.all_codewords() {
            println!("- {:?}", c)
        }

        println!();

        let code2 = LinearCode::<2, 5, 3>::with_generator(mat![
            GF(1), GF(0), GF(1), GF(0), GF(0);
            GF(1), GF(0), GF(0), GF(1), GF(0);
            GF(1), GF(0), GF(0), GF(0), GF(1);
        ]);

        for c in code2.all_codewords() {
            println!("- {:?}", c)
        }

        println!("Eq {}", code == code2);

        println!("{:?} = 2^3", code.all_codewords().len());
        println!("d(C) = {} = {}", code.distance(), code2.distance());
    }

    #[test]
    fn g_to_h() {
        let g = mat![
            GF(0), GF(1), GF(2);
            GF(2), GF(0), GF(0);
        ];

        let h = generator_to_controll(mat![
            GF::<3>(0), GF(1), GF(2);
            GF(2), GF(0), GF(0);
        ]);
        println!("{h}");

        println!(
            "{:?}",
            LinearCode::<3, 3, 2>::with_generator(g).all_codewords()
        );
        println!(
            "{:?}",
            LinearCode::<3, 3, 2>::with_control(h).all_codewords()
        );
    }

    #[test]
    fn gauss_test() {
        // let m = matrix![
        //     2.0, 1.0, -1.0, 8.0;
        //     -3.0, -1.0, 2.0, -11.0;
        //     -2.0, 1.0, 2.0, -3.0;
        // ];
        // println!("{}", gauss(m));

        let m = matrix![
            1.0, 3.0, 1.0, 9.0;
            1.0, 1.0, -1.0, 1.0;
            3.0, 11.0, 5.0, 35.0;
        ];
        println!("{}", gauss(m));
    }

    #[test]
    fn reed_solomon() {
        let _rs = LinearCode::<256, 257, 247>::reed_solomon(11).unwrap();
    }

    #[test]
    fn hamming() {
        let _lc = LinearCode::<3, 4, 2>::hamming().unwrap();
        // println!("lc")
    }

    #[test]
    fn vandermonde_test() {
        let m = vandermonde(&[10, 20, 30, 40, 50]);
        println!("{}", m)
    }

    #[test]
    fn vector_space() {
        let _e = VectorSpace::<GF<2>>::new(3).unit_vectors();
        let _a = VectorSpace::<GF<2>>::new(3).all_vectors();

        let _eq = VectorSpace::<GF<3>>::new(2).equivalence_classes();
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
