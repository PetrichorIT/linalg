#[allow(unused_imports)]
use super::*;

#[test]
fn test_fac() {
    // Base cases
    assert_eq!(1, fac(0));
    assert_eq!(1, fac(1));

    // Sample cases
    assert_eq!(2, fac(2));
    assert_eq!(24, fac(4));
}

#[test]
#[should_panic]
fn test_fac_panic() {
    let _ = fac(-1);
}

#[test]
fn test_prt_under() {
    // Edge cases
    assert_eq!(0, prt_under(0, 0));
    assert_eq!(1, prt_under(1, 0));
    assert_eq!(1, prt_under(1, 1));

    // Sample cases
    assert_eq!(5 * 4, prt_under(5, 2));
    assert_eq!(10 * 9 * 8 * 7, prt_under(10, 4));
}

#[test]
#[should_panic]
fn test_ptr_under_panic1() {
    prt_under(3, 8);
}

#[test]
#[should_panic]
fn test_ptr_under_panic2() {
    prt_under(3, -2);
}

#[test]
#[should_panic]
fn test_ptr_under_panic3() {
    prt_under(-3, 2);
}

#[test]
fn test_ptr_over() {
    // Edge cases
    assert_eq!(0, prt_over(0, 0));
    assert_eq!(1, prt_over(1, 0));
    assert_eq!(1, prt_over(1, 1));

    // Sample cases
    assert_eq!(5 * 6, prt_over(5, 2));
    assert_eq!(10 * 11 * 12 * 13, prt_over(10, 4));
}

#[test]
#[should_panic]
fn test_ptr_over_panic1() {
    prt_over(3, -2);
}

#[test]
#[should_panic]
fn test_ptr_over_panic2() {
    prt_over(-3, 2);
}
