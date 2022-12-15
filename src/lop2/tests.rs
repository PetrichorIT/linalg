use crate::matrix;

use super::defs::LOP;

#[test]
fn tabelau() {
    //
    let lop = LOP::Mixed {
        c: matrix![6.0; -8.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; 5.0;],
        a_gt: matrix![
            1.0,  1.0, -1.0;
        ],
        b_gt: matrix![-3.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        labels: Vec::new(),
    };

    let lop = lop.into_standard_form();
    println!("{lop}");
    // panic!("A");
    // return;

    // let simplex = super::Simplex::new(lop);

    // println!("{:?}", simplex.solve(1000));
    // panic!("A");
}

#[test]
fn lop_example_a() {
    //
    let lop = LOP::Mixed {
        c: matrix![-6.0; 8.0; -1.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        a_lt: matrix![
            3.0,  1.0, 0.0;
            4.0, -1.0, 0.0;
        ],
        b_lt: matrix![10.0; 5.0;],
        a_gt: matrix![
            1.0,  1.0, -1.0;
        ],
        b_gt: matrix![-3.0;],
        labels: Vec::new(),
    };

    let simplex = super::Simplex::new(lop);
    let x = simplex.solve(1000).unwrap();

    assert_eq!(x, matrix![1.25; 0.0; 4.25; 6.25; 0.0; 0.0;]);
}

#[test]
fn lop_example_b() {
    //
    let lop = LOP::CanonicalForm {
        c: matrix![1.0; -3.0; 2.0;],
        a: matrix![
            3.0, -1.0, 2.0;
            -2.0, 4.0, 0.0;
            -4.0, 3.0, 8.0;
        ],
        b: matrix![7.0; 12.0; 10.0;],
        labels: Vec::new(),
    };

    let simplex = super::Simplex::new(lop);
    let x = simplex.solve(1000).unwrap();

    assert_eq!(x, matrix![4.0;5.0;0.0;0.0;0.0;11.0;]);
}

#[test]
fn lop_mixed_to_standard_fixed_groups() {
    let lop = LOP::Mixed {
        c: matrix![6.0; -8.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; 5.0;],
        a_gt: matrix![
            1.0,  1.0, -1.0;
            2.0,  2.0, -2.0;
        ],
        b_gt: matrix![3.0; 9.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        labels: Vec::new(),
    }
    .into_standard_form();

    let LOP::StandardForm { c, a, b,.. } = lop else {
        panic!("into_standard_form() did not return standard form")
    };
    assert_eq!(c, matrix![6.0; -8.0; 1.0;]);
    assert_eq!(b, matrix![10.0; 5.0; 3.0; 9.0;]);
    assert_eq!(
        a,
        matrix![
            3.0,  1.0,  0.0, 1.0, 0.0,  0.0,  0.0;
            4.0, -1.0,  0.0, 0.0, 1.0,  0.0,  0.0;
            1.0,  1.0, -1.0, 0.0, 0.0, -1.0,  0.0;
            2.0,  2.0, -2.0, 0.0, 0.0,  0.0, -1.0;
        ]
    );

    let lop = LOP::Mixed {
        c: matrix![6.0; -8.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
        ],
        b_lt: matrix![10.0;],
        a_gt: matrix![
            1.0,  1.0, -1.0;
        ],
        b_gt: matrix![3.0;],
        a_eq: matrix![
            4.0, -1.0,  0.0;
            2.0,  2.0, -2.0;
        ],
        b_eq: matrix![5.0; 9.0;],
        labels: Vec::new(),
    }
    .into_standard_form();

    let LOP::StandardForm { c, a, b,.. } = lop else {
        panic!("into_standard_form() did not return standard form")
    };
    assert_eq!(c, matrix![6.0; -8.0; 1.0;]);
    assert_eq!(b, matrix![5.0; 9.0; 10.0; 3.0;]);
    assert_eq!(
        a,
        matrix![
            4.0, -1.0,  0.0,  0.0,  0.0;
            2.0,  2.0, -2.0,  0.0,  0.0;
            3.0,  1.0,  0.0,  1.0,  0.0;
            1.0,  1.0, -1.0,  0.0, -1.0;
        ]
    );
}

#[test]
fn lop_mixed_to_standard_mixed_groups() {
    let lop = LOP::Mixed {
        c: matrix![6.0; -8.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; -5.0;],
        a_gt: matrix![
            1.0,  1.0, -1.0;
            2.0,  2.0, -2.0;
        ],
        b_gt: matrix![-3.0; 9.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        labels: Vec::new(),
    }
    .into_standard_form();

    let LOP::StandardForm { c, a, b, .. } = lop else {
        panic!("into_standard_form() did not return standard form")
    };
    assert_eq!(c, matrix![6.0; -8.0; 1.0;]);
    assert_eq!(b, matrix![10.0; 3.0; 9.0; 5.0;]);
    assert_eq!(
        a,
        matrix![
            3.0,  1.0,  0.0, 1.0, 0.0,  0.0,  0.0;
           -1.0, -1.0,  1.0, 0.0, 1.0,  0.0,  0.0;
            2.0,  2.0, -2.0, 0.0, 0.0, -1.0,  0.0;
           -4.0,  1.0,  0.0, 0.0, 0.0,  0.0, -1.0;
        ]
    );

    let lop = LOP::Mixed {
        c: matrix![6.0; -8.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; -5.0;],
        a_gt: matrix![
            1.0,  1.0, -1.0;
            2.0,  2.0, -2.0;
        ],
        b_gt: matrix![-3.0; 9.0;],
        a_eq: matrix![
            8.0, 8.0, 8.0;
        ],
        b_eq: matrix![-8.0;],
        labels: Vec::new(),
    }
    .into_standard_form();

    let LOP::StandardForm { c, a, b, .. } = lop else {
        panic!("into_standard_form() did not return standard form")
    };
    assert_eq!(c, matrix![6.0; -8.0; 1.0;]);
    assert_eq!(b, matrix![8.0; 10.0; 3.0; 9.0; 5.0;]);
    assert_eq!(
        a,
        matrix![
           -8.0, -8.0, -8.0, 0.0, 0.0,  0.0,  0.0;
            3.0,  1.0,  0.0, 1.0, 0.0,  0.0,  0.0;
           -1.0, -1.0,  1.0, 0.0, 1.0,  0.0,  0.0;
            2.0,  2.0, -2.0, 0.0, 0.0, -1.0,  0.0;
           -4.0,  1.0,  0.0, 0.0, 0.0,  0.0, -1.0;
        ]
    );
}

#[test]
#[should_panic]
fn lop_mixed_to_standard_dim_missmatch() {
    let _lop = LOP::Mixed {
        c: matrix![6.0; -8.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; -5.0;],
        a_gt: matrix![
            1.0,  1.0;
            2.0,  2.0;
        ],
        b_gt: matrix![-3.0; 9.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        labels: Vec::new(),
    }
    .into_standard_form();
}

#[test]
#[should_panic]
fn lop_mixed_to_standard_c_missmatch() {
    let _lop = LOP::Mixed {
        c: matrix![6.0; 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; -5.0;],
        a_gt: matrix![
            1.0,  1.0, 0.0;
            2.0,  2.0, 0.0;
        ],
        b_gt: matrix![-3.0; 9.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        labels: Vec::new(),
    }
    .into_standard_form();
}

#[test]
#[should_panic]
fn lop_mixed_to_standard_b_missmatch() {
    let _lop = LOP::Mixed {
        c: matrix![6.0; -8.0, 1.0;],
        a_lt: matrix![
            3.0,  1.0,  0.0;
            4.0, -1.0,  0.0;
        ],
        b_lt: matrix![10.0; -5.0; 1.0;],
        a_gt: matrix![
            1.0,  1.0, 0.0;
            2.0,  2.0, 0.0;
        ],
        b_gt: matrix![-3.0; 9.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        labels: Vec::new(),
    }
    .into_standard_form();
}

#[test]
fn lop_canonical_into_standard() {
    let lop = LOP::CanonicalForm {
        c: matrix![1;2;],
        a: matrix![1,0;2,3;],
        b: matrix![1;-2;],
        labels: Vec::new(),
    };

    println!("{:?}", lop);

    let s = lop.into_standard_form();
    println!("{:?}", s);

    // Simplex::new(s).solve(1000).unwrap();
}

// #[test]
// fn lop_mixed_into_standard() {
//     let lop = LOP::Mixed {
//         c: matrix![1.0; 2.0;],
//         a: matrix![-3.0, 2.0;9.0,8.0;],
//         b: matrix![3.0;100.0;],
//         a_eq: matrix![],
//         b_eq: matrix![],
//         x_bounds: vec![
//             LOPParameterBound::GreaterOrEqualThan(0.0),
//             LOPParameterBound::GreaterOrEqualThan(0.0),
//         ],
//     };

//     println!("{}", lop);
//     println!("{}", lop.into_standard_form())
// }

// #[test]
// fn lop_mixed_into_standard2() {
//     let lop = LOP::Mixed {
//         c: matrix![1.0; 2.0;],
//         a: matrix![-3.0, 2.0;9.0,8.0;],
//         b: matrix![3.0;100.0;],
//         a_eq: matrix![],
//         b_eq: matrix![],
//         x_bounds: vec![
//             LOPParameterBound::GreaterOrEqualThan(0.0),
//             LOPParameterBound::LesserOrEqualThan(0.0),
//         ],
//     };

//     println!("{}", lop);
//     println!("{}", lop.into_standard_form())
// }
