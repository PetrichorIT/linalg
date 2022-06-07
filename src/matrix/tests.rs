#[allow(unused_imports)]
use super::{Matrix, MatrixLayout};
#[allow(unused_imports)]
use crate::{mat, sets::GF};

// [MatrixLayout] Tests

#[test]
fn layout_test_inline() {
    let layout = MatrixLayout { rows: 3, cols: 4 };

    assert_eq!(layout.rows(), layout.rows);
    assert_eq!(layout.cols(), layout.cols);
    assert_eq!(layout.size(), layout.rows * layout.cols);
    assert_eq!(layout.is_square(), false);

    let layout = MatrixLayout { rows: 1, cols: 1 };

    assert_eq!(layout.is_square(), true);
    assert_eq!(layout.is_rowvec(), true);
    assert_eq!(layout.is_colvec(), true);
    assert_eq!(layout.is_vec(), true);

    let layout = MatrixLayout { rows: 5, cols: 1 };

    assert_eq!(layout.is_square(), false);
    assert_eq!(layout.is_rowvec(), false);
    assert_eq!(layout.is_colvec(), true);
    assert_eq!(layout.is_vec(), true);

    let mut layout = layout;
    layout.transpose();

    assert_eq!(layout.is_square(), false);
    assert_eq!(layout.is_rowvec(), true);
    assert_eq!(layout.is_colvec(), false);
    assert_eq!(layout.is_vec(), true);

    let layout = layout.transposed();

    assert_eq!(layout.is_square(), false);
    assert_eq!(layout.is_rowvec(), false);
    assert_eq!(layout.is_colvec(), true);
    assert_eq!(layout.is_vec(), true);
}

#[test]
fn layout_safe_constructors() {
    let layout = MatrixLayout::new(7, 2);
    assert_eq!(layout, MatrixLayout { rows: 7, cols: 2 });

    let layout = MatrixLayout::square(4);
    assert_eq!(layout, MatrixLayout { rows: 4, cols: 4 });
    assert_eq!(layout.is_square(), true);

    let layout: MatrixLayout = (3, 2).into();
    assert_eq!(layout, MatrixLayout { rows: 3, cols: 2 })
}

#[test]
fn layout_assumed_from_success() {
    let v = vec![vec![1, 2], vec![3, 4]];

    let layout = MatrixLayout::assumed_from(&v);
    assert!(layout.is_ok());
    assert_eq!(layout.unwrap(), MatrixLayout { rows: 2, cols: 2 });
}

#[test]
fn layout_assumed_from_failure() {
    let v = vec![vec![1, 2], vec![3, 4, 5]];

    let layout = MatrixLayout::assumed_from(&v);
    assert!(layout.is_err());
    assert_eq!(layout.unwrap_err(), "Structured data has irregulare format");
}

// [Matrix] Tests

#[test]
fn matrix_inline() {
    let buf = vec![1, 2, 3, 4];
    let layout = MatrixLayout::new(2, 2);

    let matrix = Matrix {
        raw: buf.clone(),
        layout: layout.clone(),
    };

    assert_eq!(*matrix.layout(), layout);
    assert_eq!(*matrix.raw(), buf);
    assert_eq!(matrix.size(), layout.size());
}

#[test]
fn matrix_new_succes() {
    let buf = vec![1, 2, 3, 4];

    let matrix = Matrix::new((2, 2), buf);
    assert_eq!(matrix.layout, MatrixLayout::new(2, 2));
    assert_eq!(matrix.raw, vec![1, 2, 3, 4]);

    // Should be possible with unsafe buffer

    let mut buf = Vec::<usize>::with_capacity(9);
    unsafe { buf.set_len(9) }

    let matrix = Matrix::new((3, 3), buf);
    assert_eq!(matrix.layout, MatrixLayout::new(3, 3));
}

#[test]
#[should_panic]
fn matrix_new_failure() {
    let buf = vec![1, 2, 3, 4];
    let _matrix = Matrix::new((2, 3), buf);
}

/*
Following the ::new() tests we will simply assume that the colvec/rowvec
are defined correctly.
Additionaly we will not test unsafe functions.
*/

#[test]
fn matrix_symmerty() {
    let matrix = mat![
        1, 2, 3;
        2, 3, 5;
        3, 5, 0;
    ];
    assert_eq!(matrix.is_symetric(), true);

    let matrix = mat![
        1, 3, 3;
        2, 3, 5;
        3, 5, 0;
    ];
    assert_eq!(matrix.is_symetric(), false);

    let matrix = mat![
        1, 3, 3, 1;
        2, 3, 5, 1;
        3, 5, 0, 1;
    ];
    assert_eq!(matrix.is_symetric(), false);
}

#[test]
fn matrix_upper_triag() {
    let matrix = mat![
        1, 2, 3;
        0, 1, 2;
        0, 0, 1;
    ];
    assert_eq!(matrix.is_upper_triag(), true);

    let matrix = mat![
        1, 2, 3;
        0, 1, 2;
        0, 1, 1;
    ];
    assert_eq!(matrix.is_upper_triag(), false);

    let matrix = mat![
        1, 2, 3, 1;
        0, 1, 2, 1;
        0, 0, 1, 1;
    ];
    assert_eq!(matrix.is_upper_triag(), true);

    let matrix = mat![
        1, 2, 3;
        0, 1, 2;
        0, 0, 1;
        0, 0, 0;
    ];
    assert_eq!(matrix.is_upper_triag(), true);
}

#[test]
fn matrix_lower_triag() {
    let matrix = mat![
        1, 0, 0;
        2, 1, 0;
        3, 5, 1;
    ];
    assert_eq!(matrix.is_lower_triag(), true);

    let matrix = mat![
        1, 0, 1;
        2, 1, 0;
        3, 5, 1;
    ];
    assert_eq!(matrix.is_lower_triag(), false);

    let matrix = mat![
        1, 0, 0;
        2, 1, 0;
        3, 5, 1;
        1, 1, 1;
    ];
    assert_eq!(matrix.is_lower_triag(), true);

    let matrix = mat![
        1, 0, 0, 0;
        2, 1, 0, 0;
        3, 5, 1, 0;
    ];
    assert_eq!(matrix.is_lower_triag(), true);
}

#[test]
fn matrix_is_diag() {
    let matrix = mat![
        1, 0, 0;
        0, 2, 0;
        0, 0, 3;
    ];
    assert_eq!(matrix.is_diag(), true);

    let matrix = mat![
        1, 0, 1;
        0, 2, 0;
        0, 0, 3;
    ];
    assert_eq!(matrix.is_diag(), false);

    let matrix = mat![
        1, 0, 0, 0;
        0, 2, 0, 0;
        0, 0, 3, 0;
    ];
    assert_eq!(matrix.is_diag(), true);

    let matrix = mat![
        1, 0, 0;
        0, 2, 0;
        0, 0, 3;
        0, 0, 0;
    ];
    assert_eq!(matrix.is_diag(), true);
}

/*
TODO: Tests and actual use cases for from_slice(es)

Allso we assume that fill, diag, eye and zeroed will in fact work
*/

#[test]
fn matrix_resize() {
    let mut matrix = mat![
        1, 2, 3;
        4, 5, 6;
        7, 8, 9;
    ];

    // memory efficient varian
    matrix.resize((4, 3));
    assert_eq!(*matrix.layout(), MatrixLayout::new(4, 3));

    assert_eq!(matrix[(1, 1)], 5);
    assert_eq!(matrix[(1, 0)], 4);
    assert_eq!(matrix[(0, 1)], 2);

    assert_eq!(matrix[(3, 0)], 0);

    let mut matrix = mat![
        1, 2, 3;
        4, 5, 6;
        7, 8, 9;
    ];

    // memory in-efficient varian
    matrix.resize((3, 4));
    assert_eq!(*matrix.layout(), MatrixLayout::new(3, 4));

    assert_eq!(matrix[(1, 1)], 5);
    assert_eq!(matrix[(1, 0)], 4);
    assert_eq!(matrix[(0, 1)], 2);

    assert_eq!(matrix[(0, 3)], 0);
}

#[test]
fn matrix_partial_ord() {
    // Restricted to vectors

    let v1 = mat![ 1, 2, 3;];
    let v2 = mat![1, 3, 4;];
    let v3 = mat![9, 9, 9;];

    assert!(v1 <= v2);
    assert!(!(v1 < v2));

    assert!(v1 < v3);
    assert!(v2 < v3);

    assert!(v3 <= v3);
    assert!(!(v3 < v3))
}

/*
TODO: Make basic vector opertaion test for add bit* sub ...
 */

#[test]
fn matrix_trace() {
    let matrix = Matrix::diag(vec![1, 2, 3, 5]);
    assert_eq!(matrix.trace(), 1 + 2 + 3 + 5);

    let mut matrix = matrix;
    matrix[(1, 2)] = 4;
    assert_eq!(matrix.trace(), 1 + 2 + 3 + 5);
}

#[test]
#[should_panic]
fn matrix_trace_panic() {
    let matrix = mat![
        1, 2, 3;
        1, 2, 3;
    ];

    let _a = matrix.trace();
}

#[test]
fn matrix_transpose() {
    let matrix = mat![
        1, 2, 3;
        4, 5, 6;
    ];

    let result = matrix.transposed();
    assert_eq!(*result.layout(), MatrixLayout::new(3, 2));
    assert_eq!(*result.raw(), vec![1, 4, 2, 5, 3, 6]);

    let mut matrix = matrix;
    matrix.transpose();

    assert_eq!(matrix, result);

    // check transpositional vector optiomizations;

    let v = mat![1, 2, 3;];
    let mut vt = v.clone();
    vt.transpose();

    assert_eq!(vt, v.transposed());

    // check into operator

    let mut buf = Matrix::zeroed((3, 1));
    v.transpose_into(&mut buf);

    assert_eq!(buf, vt);
}

#[test]
fn matrix_scalar_mul() {
    let matrix = mat![
        1, 2, 3;
        4, 5, 6;
    ];

    // Scale and Mul
    let double = matrix.clone() * 2usize;
    assert_eq!(*double.layout(), *matrix.layout());
    assert_eq!(double[(1, 1)], 2 * matrix[(1, 1)]);

    // scale with MulAssign
    let mut mutmat = matrix.clone();
    mutmat *= 8;
    assert_eq!(*mutmat.layout(), *matrix.layout());
    assert_eq!(mutmat[(1, 1)], 8 * matrix[(1, 1)]);

    // scalar (no clone needed ref based)
    let result = matrix.scalar(9);
    assert_eq!(*result.layout(), *matrix.layout());
    assert_eq!(result[(1, 1)], 9 * matrix[(1, 1)]);
}

#[test]
fn matrix_mmul() {
    let a = mat![
        1, 2, 3;
        4, 5, 6;
    ];
    let b = mat![
        1, 2;
        3,4;
        5, 6;
    ];

    let c = a.clone() * b.clone();
    assert_eq!(*c.layout(), MatrixLayout::new(2, 2));
    assert_eq!(c[(0, 0)], 1 + 6 + 15);

    let d = b * a;
    assert_eq!(*d.layout(), MatrixLayout::new(3, 3));
    assert_eq!(d[(0, 0)], 9);
}

#[test]
#[should_panic]
fn matrix_mmul_panic() {
    let a = mat![1, 2, 3;4, 5, 6;];
    let b = mat![1, 2; 3,4;];

    let _c = a * b;
}

#[allow(unused_imports)]
use std::ops::Neg;

#[test]
fn matrix_neg() {
    let m = mat![
        1, 2, 3isize;
        4, 5, 6;
    ];

    // checking direct impl
    let n = m.clone().neg();
    assert_eq!(*n.layout(), *m.layout());
    for k in 0..n.layout().size() {
        assert!(n[k] * -1 == m[k])
    }

    let new_m = (&n).neg();
    assert_eq!(m, new_m);
}

#[test]
fn matrix_extract() {
    let m = mat![
        1, 2, 3, 4, 5;
        2, 3, 4, 5, 6;
        3, 4, 5, 6, 7;
        4, 5, 6, 7, 8;
        5, 6, 7, 8, 9;
    ];

    let sb = m.extract(2.., 3..);
    assert_eq!(
        sb,
        mat![
            6, 7;
            7, 8;
            8, 9;
        ]
    );

    let sb = m.extract(3.., ..);
    assert_eq!(
        sb,
        mat![
            4, 5, 6, 7, 8;
            5, 6, 7, 8, 9;
        ]
    );

    let sb = m.extract(3.., ..=2);
    assert_eq!(
        sb,
        mat![
            4, 5, 6;
            5, 6, 7;
        ]
    );

    let sb = m.extract(..3, 0..=2);
    assert_eq!(
        sb,
        mat![
            1, 2, 3;
            2, 3, 4;
            3, 4, 5;
        ]
    );
}

#[test]
fn matrix_insert() {
    let mut matrix = Matrix::zeroed((5, 5));

    let mut insert = Matrix::fill((3, 3), 8);

    matrix.insert(2.., 2.., &insert);
    assert_eq!(
        matrix,
        mat![
            0, 0, 0, 0, 0;
            0, 0, 0, 0, 0;
            0, 0, 8, 8, 8;
            0, 0, 8, 8, 8;
            0, 0, 8, 8, 8;
        ]
    );

    insert *= 2;

    matrix.insert(0..=1, 0.., &insert);
    assert_eq!(
        matrix,
        mat![
            16, 16, 16, 0, 0;
            16, 16, 16, 0, 0;
            0, 0, 8, 8, 8;
            0, 0, 8, 8, 8;
            0, 0, 8, 8, 8;
        ]
    );

    matrix.insert(3.., 3.., &insert);
    assert_eq!(
        matrix,
        mat![
            16, 16, 16, 0, 0;
            16, 16, 16, 0, 0;
            0, 0, 8, 8, 8;
            0, 0, 8, 16, 16;
            0, 0, 8, 16, 16;
        ]
    );
}

#[test]
fn matrix_gf_2() {
    let matrix: Matrix<GF<2>> = mat![
        GF(1), GF(1), GF(0), GF(0);
        GF(1), GF(0), GF(0), GF(1);
    ];

    let v: Matrix<GF<2>> = mat![
        GF(1), GF(0);
    ];

    println!("{}", matrix);
    println!("{}", v);
    println!("{}", v * matrix)
}

#[test]
fn series_5_task_1() {
    let input: Vec<Matrix<GF<2>>> = vec![
        mat![ GF(0), GF(0), GF(0), GF(0); ],
        mat![ GF(0), GF(0), GF(0), GF(1); ],
        mat![ GF(0), GF(0), GF(1), GF(0); ],
        mat![ GF(0), GF(0), GF(1), GF(1); ],
        mat![ GF(0), GF(1), GF(0), GF(0); ],
        mat![ GF(0), GF(1), GF(0), GF(1); ],
        mat![ GF(0), GF(1), GF(1), GF(0); ],
        mat![ GF(0), GF(1), GF(1), GF(1); ],
        mat![ GF(1), GF(0), GF(0), GF(0); ],
        mat![ GF(1), GF(0), GF(0), GF(1); ],
        mat![ GF(1), GF(0), GF(1), GF(0); ],
        mat![ GF(1), GF(0), GF(1), GF(1); ],
        mat![ GF(1), GF(1), GF(0), GF(0); ],
        mat![ GF(1), GF(1), GF(0), GF(1); ],
        mat![ GF(1), GF(1), GF(1), GF(0); ],
        mat![ GF(1), GF(1), GF(1), GF(1); ],
    ];

    let generator: Matrix<GF<2>> = mat![
        GF(1), GF(0), GF(0), GF(0), GF(1), GF(1), GF(1);
        GF(0), GF(1), GF(0), GF(0), GF(1), GF(1), GF(0);
        GF(0), GF(0), GF(1), GF(0), GF(1), GF(0), GF(1);
        GF(0), GF(0), GF(0), GF(1), GF(0), GF(1), GF(1);
    ];

    for value in input {
        println!("{} ==> {}", &value, &value * &generator)
    }
}

#[test]
fn series_5_task_3() {
    let input: Vec<Matrix<GF<3>>> = vec![
        mat! [ GF(0), GF(0); ],
        mat! [ GF(0), GF(1); ],
        mat! [ GF(0), GF(2); ],
        mat! [ GF(1), GF(0); ],
        mat! [ GF(1), GF(1); ],
        mat! [ GF(1), GF(2); ],
        mat! [ GF(2), GF(0); ],
        mat! [ GF(2), GF(1); ],
        mat! [ GF(2), GF(2); ],
    ];

    let gen: Matrix<GF<3>> = mat![
        GF(2), GF(1), GF(2), GF(0);
        GF(1), GF(1), GF(0), GF(1);
    ];

    for value in input {
        println!("{} ==> {}", &value, &value * &gen);
    }
}
