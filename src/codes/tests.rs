#[allow(unused)]
use super::{permut::Permutable, *};
#[allow(unused)]
use crate::{mat, matrix::Matrix, sets::GF};

#[test]
fn series_3_task_1() {
    let code = Code {
        words: vec![
            Word(vec![0; 5]),
            Word(vec![1, 1, 1, 1, 0]),
            Word(vec![0, 1, 1, 1, 1]),
            Word(vec![1, 0, 0, 0, 1]),
        ],
    };

    println!("d(C) = {}", code.distance());

    println!("w1 => {:?}", code.match_for(Word(vec![1, 1, 1, 1, 1])));
    println!("w1 => {:?}", code.match_for(Word(vec![1, 0, 1, 0, 1])));
    println!("w1 => {:?}", code.match_for(Word(vec![0, 1, 1, 1, 0])));
}

#[test]
fn series_3_task_3() {
    let code = Code::new(&[0, 1], 5).repeat(3);

    println!(
        "010101101010110 => {:?}",
        code.match_for(Word(vec![0usize, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0]))
    );

    println!(
        "010100111101110 => {:?}",
        code.match_for(Word(vec![0usize, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0]))
    );
}

#[test]
fn series_3_task_6() {
    let code = Code {
        words: vec![
            Word(vec![0, 0, 0, 0]),
            Word(vec![0, 1, 1, 2]),
            Word(vec![0, 2, 2, 1]),
            Word(vec![1, 1, 0, 1]),
            Word(vec![1, 2, 1, 0]),
            Word(vec![1, 0, 2, 2]),
            Word(vec![2, 2, 0, 2]),
            Word(vec![2, 0, 1, 1]),
            Word(vec![2, 1, 2, 0]),
        ],
    };

    println!("d(C) = {}", code.distance())
}

#[test]
fn permut() {
    for p in [1, 2, 3, 4].permut() {
        println!("{:?}", p)
    }
}

#[test]
fn series_4_a_2() {
    let code = Code::new(&[0, 1], 6);

    'main: for a in 0..code.len() {
        println!("{}", a);
        for b in a..code.len() {
            if a == b {
                continue;
            }
            for c in b..code.len() {
                if b == c {
                    continue;
                }
                for d in c..code.len() {
                    if c == d {
                        continue;
                    }
                    for e in d..code.len() {
                        if d == e {
                            continue;
                        }

                        //
                        // generate subset
                        let new_words = Code {
                            words: vec![
                                code.words[a].clone(),
                                code.words[b].clone(),
                                code.words[c].clone(),
                                code.words[d].clone(),
                                code.words[e].clone(),
                            ],
                        };

                        if new_words.distance() >= 3 {
                            dbg!(new_words);
                            break 'main;
                        }
                    }
                }
            }
        }
    }
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

    for value in &input {
        println!("{} ==> {}", value, value * &generator)
    }

    let input: Vec<Word<GF<2>>> = input
        .into_iter()
        .map(|value| value * &generator)
        .map(|matrix| Word(matrix.to_vec()))
        .collect();

    let code = Code { words: input };
    println!(
        "{:?}",
        code.match_for(Word(vec![GF(0), GF(0), GF(0), GF(1), GF(1), GF(0), GF(0)]))
    );
    println!(
        "{:?}",
        code.match_for(Word(vec![GF(0), GF(0), GF(0), GF(1), GF(0), GF(1), GF(1)]))
    );
    println!(
        "{:?}",
        code.match_for(Word(vec![GF(1), GF(1), GF(1), GF(0), GF(1), GF(0), GF(0)]))
    );
    println!(
        "{:?}",
        code.match_for(Word(vec![GF(1), GF(1), GF(1), GF(1), GF(1), GF(1), GF(0)]))
    );
}

#[test]
fn series_5_task_2() {
    let input: Vec<Matrix<GF<2>>> = vec![
        mat! [ GF(0), GF(0), GF(0); ],
        mat! [ GF(0), GF(0), GF(1); ],
        mat! [ GF(0), GF(1), GF(0); ],
        mat! [ GF(0), GF(1), GF(1); ],
        mat! [ GF(1), GF(0), GF(0); ],
        mat! [ GF(1), GF(0), GF(1); ],
        mat! [ GF(1), GF(1), GF(0); ],
        mat! [ GF(1), GF(1), GF(1); ],
    ];

    let gen: Matrix<GF<2>> = mat![
        GF(1), GF(1), GF(1), GF(0), GF(1);
        GF(1), GF(0), GF(1), GF(1), GF(0);
        GF(1), GF(1), GF(0), GF(1), GF(0);
    ];

    let input: Vec<Word<GF<2>>> = input
        .into_iter()
        .map(|value| value * &gen)
        .map(|matrix| Word(matrix.to_vec()))
        .collect();

    let g: Vec<usize> = input
        .iter()
        .map(|w| {
            w.0.iter()
                .fold(0, |acc, c| acc + if *c == GF(0) { 0 } else { 1 })
        })
        .collect();

    println!("{:?}", g);

    let code = Code { words: input };
    println!("{:?}", code.min_distance());
}

#[test]
fn series_5_task_3() {
    let input: Vec<Matrix<GF<3>>> = Code::new(&[GF(0), GF(1), GF(2)], 2)
        .into_iter()
        .map(|w| Matrix::rowvec(w.0))
        .collect();

    // let input: Vec<Matrix<GF<3>>> = vec![
    //     mat! [ GF(0), GF(0); ],
    //     mat! [ GF(0), GF(1); ],
    //     mat! [ GF(0), GF(2); ],
    //     mat! [ GF(1), GF(0); ],
    //     mat! [ GF(1), GF(1); ],
    //     mat! [ GF(1), GF(2); ],
    //     mat! [ GF(2), GF(0); ],
    //     mat! [ GF(2), GF(1); ],
    //     mat! [ GF(2), GF(2); ],
    // ];

    let gen: Matrix<GF<3>> = mat![
        GF(2), GF(1), GF(2), GF(0);
        GF(1), GF(1), GF(0), GF(1);
    ];

    let input: Vec<Word<GF<3>>> = input
        .into_iter()
        .map(|value| value * &gen)
        .map(|matrix| dbg!(Word(matrix.to_vec())))
        .collect();

    let code = Code { words: input };
    println!("{}", code.distance());
}

#[test]
fn series_5_task_4() {
    let alphabeth: Code<GF<3>> = Code::new(&[GF(0), GF(1), GF(2)], 2);
    let input: Vec<Matrix<GF<3>>> = alphabeth.into_iter().map(|v| Matrix::rowvec(v.0)).collect();

    let gen: Matrix<GF<3>> = mat![
        GF(0), GF(1), GF(1), GF(2);
        GF(1), GF(0), GF(1), GF(2);
    ];

    let input: Vec<Word<GF<3>>> = input
        .into_iter()
        .map(|value| value * &gen)
        .map(|matrix| Word(matrix.to_vec()))
        .collect();

    let code = Code { words: input };
    println!("{:?}", code.min_distance());
}

#[test]
fn series_5_task_5() {
    let alphabeth: Code<GF<2>> = Code::new(&[GF(0), GF(1)], 4);
    let Code { words } = alphabeth;
    let input: Vec<Matrix<GF<2>>> = words.into_iter().map(|v| Matrix::rowvec(v.0)).collect();

    let gen: Matrix<GF<2>> = mat![
        GF(1), GF(0), GF(0), GF(0), GF(1), GF(1), GF(1);
        GF(0), GF(1), GF(0), GF(0), GF(1), GF(0), GF(1);
        GF(0), GF(0), GF(1), GF(0), GF(0), GF(1), GF(1);
        GF(0), GF(0), GF(0), GF(1), GF(1), GF(1), GF(0);
    ];

    let input: Vec<Word<GF<2>>> = input
        .into_iter()
        .map(|value| value * &gen)
        .map(|matrix| Word(matrix.to_vec()))
        .collect();

    let code = Code { words: input };
    println!("{}", code.distance());
    println!(
        "{:?}",
        code.match_for(Word(vec![GF(1), GF(0), GF(1), GF(0), GF(0), GF(1), GF(0)]))
    );
    println!(
        "{:?}",
        code.match_for(Word(vec![GF(1), GF(1), GF(1), GF(1), GF(1), GF(0), GF(0)]))
    );
}
