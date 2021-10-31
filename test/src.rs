use linalg::prelude::{stirling1, Stirling1Tbl};
use random::Source;
use std::time::Instant;

fn main() {
    let n = 100_000;
    let upper_bound = 20;

    type T = usize;

    let t0 = Instant::now();
    let mut nums = Vec::<(T, T)>::with_capacity(n);

    let mut source = random::default().seed([42, 69]);
    for _ in 0..n {
        let r0 = source.read::<T>() % upper_bound;
        let r1 = source.read::<T>() % upper_bound;
        nums.push((r0.max(r1), r1.min(r0)));
    }

    let t1 = Instant::now();

    let mut tbl = Stirling1Tbl::new(1);

    for &(n, k) in &nums {
        let _ = stirling1(n, k);
    }

    let t2 = Instant::now();

    for &(n, k) in &nums {
        let _ = tbl.get(n, k);
    }

    let t3 = Instant::now();

    // Eval

    let gen_time = t1.duration_since(t0);
    let binom_time = t2.duration_since(t1);
    let binom2_time = t3.duration_since(t2);

    println!("=== Results ===");
    println!("(t) gentime := {}ns", gen_time.as_nanos());
    println!("(t) binom_1 := {}ns", binom_time.as_nanos());
    println!("(t) binom_2 := {}ns", binom2_time.as_nanos());

    let prct = (binom_time.as_nanos() as f64) / (binom2_time.as_nanos() as f64);

    println!("(p) {}x", prct - 1.0);
}
