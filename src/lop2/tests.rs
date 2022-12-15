use crate::matrix;

use super::defs::{LOPParameterBound, LOP};

#[test]
fn tabelau() {
    let lop = LOP::CanonicalForm {
        c: matrix![20.0; 15.0;],
        a: matrix![
            1.0, 1.0;
            9.0, 5.0;
            2.0, 1.0;
        ],
        b: matrix![7.0; 45.0; 8.0;],
    };

    let simplex = super::Simplex::new(lop);
    println!("{}", simplex);

    println!("{:?}", simplex.solve(1000));
    panic!("A");
}

#[test]
fn tabelau2() {
    let lop = LOP::CanonicalForm {
        c: matrix![-1.0; 3.0; -2.0;],
        a: matrix![
            3.0, -1.0, 2.0;
            -2.0, 4.0, 0.0;
            -4.0, 3.0, 8.0;
        ],
        b: matrix![7.0; 12.0; 10.0;],
    };

    let simplex = super::Simplex::new(lop);
    println!("{}", simplex);

    println!("{:?}", simplex.solve(1000));
    panic!("A");
}

#[test]
fn lop_canonical_into_standard() {
    let lop = LOP::CanonicalForm {
        c: matrix![1;2;],
        a: matrix![1,0;2,3;],
        b: matrix![1;-2;],
    };

    println!("{:?}", lop);

    let s = lop.into_standard_form();
    println!("{:?}", s);

    // Simplex::new(s).solve(1000).unwrap();
}

#[test]
fn lop_mixed_into_standard() {
    let lop = LOP::Mixed {
        c: matrix![1.0; 2.0;],
        a: matrix![-3.0, 2.0;9.0,8.0;],
        b: matrix![3.0;100.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        x_bounds: vec![
            LOPParameterBound::GreaterOrEqualThan(0.0),
            LOPParameterBound::GreaterOrEqualThan(0.0),
        ],
    };

    println!("{}", lop);
    println!("{}", lop.into_standard_form())
}

#[test]
fn lop_mixed_into_standard2() {
    let lop = LOP::Mixed {
        c: matrix![1.0; 2.0;],
        a: matrix![-3.0, 2.0;9.0,8.0;],
        b: matrix![3.0;100.0;],
        a_eq: matrix![],
        b_eq: matrix![],
        x_bounds: vec![
            LOPParameterBound::GreaterOrEqualThan(0.0),
            LOPParameterBound::LesserOrEqualThan(0.0),
        ],
    };

    println!("{}", lop);
    println!("{}", lop.into_standard_form())
}
