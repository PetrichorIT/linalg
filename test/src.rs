use linalg::prelude::*;

fn main() {
    let mut p1 = Polynom::new(vec![0; 85]);
    *p1.coefficent(1) = 1;
    *p1.coefficent(9) = 1;
    *p1.coefficent(25) = 1;
    *p1.coefficent(49) = 1;
    *p1.coefficent(81) = 1;

    let mut p2 = Polynom::new(vec![0; 5]);
    *p2.coefficent(1) = -1;
    *p2.coefficent(3) = 1;

    let (val, r) = pdiv(p1, p2);
    println!("{} {}", val, r)
}
