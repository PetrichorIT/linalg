use std::{
    borrow::Borrow,
    fmt::{Debug, Display},
};

use crate::prelude::binom;

///
/// A Word defined as a finite sequence of characters from an alphabeth A.
///
#[derive(Clone, PartialEq, Eq)]
pub struct Word<A>(pub Vec<A>);

impl<A: PartialEq> Word<A> {
    pub fn distance(&self, other: &Self) -> usize {
        let mut result = self.0.len().abs_diff(other.0.len());
        for (lhs, rhs) in self.0.iter().zip(other.0.iter()) {
            if lhs != rhs {
                result += 1;
            }
        }

        result
    }
}

impl<A: Debug> Debug for Word<A> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Word(")?;
        for c in self.0.iter() {
            write!(f, "{:?}", c)?
        }

        write!(f, ")")
    }
}

impl<A: Display> Display for Word<A> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Word(")?;
        for c in self.0.iter() {
            write!(f, "{}", c)?
        }

        write!(f, ")")
    }
}

///
/// A Code defined as a set of codewords.
///
/// In most cases all codewords should have the same length, but this
/// is not a nessecary invariant of this struct.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Code<A> {
    pub words: Vec<Word<A>>,
}

impl<A: PartialEq> Code<A> {
    pub fn len(&self) -> usize {
        self.words.len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn new(alphabeth: &[A], k: usize) -> Self
    where
        A: Clone,
    {
        let mut words = Vec::new();
        fn r<A>(alphabeth: &[A], rem: Word<A>, k: usize, words: &mut Vec<Word<A>>)
        where
            A: Clone,
        {
            if k == 0 {
                words.push(rem)
            } else {
                for a in alphabeth {
                    let mut cp = rem.clone();
                    cp.0.push(a.clone());
                    r(alphabeth, cp, k - 1, words)
                }
            }
        }

        r(alphabeth, Word(Vec::new()), k, &mut words);

        Self { words }
    }

    pub fn repeat(self, r: usize) -> Self
    where
        A: Copy,
    {
        Self {
            words: self
                .words
                .into_iter()
                .map(|w| Word(w.0.repeat(r)))
                .collect(),
        }
    }

    pub fn distance(&self) -> usize {
        let mut min = usize::MAX;
        for i in 0..self.words.len() {
            for j in (i + 1)..self.words.len() {
                let d = self.words[i].distance(&self.words[j]);
                if d < min {
                    min = d;
                }
            }
        }

        min
    }

    pub fn min_distance(&self) -> (&Word<A>, &Word<A>, usize) {
        let mut min = usize::MAX;
        let mut a = 0;
        let mut b = 0;
        for i in 0..self.words.len() {
            for j in (i + 1)..self.words.len() {
                let d = self.words[i].distance(&self.words[j]);
                if d < min {
                    min = d;
                    a = i;
                    b = j;
                }
            }
        }

        (&self.words[a], &self.words[b], min)
    }

    pub fn match_for(&self, word: impl Borrow<Word<A>>) -> (&Word<A>, usize) {
        let word = word.borrow();

        let mut d = usize::MAX;
        let mut r = None;
        for v in self.words.iter() {
            let n_d = v.distance(word);
            if n_d < d {
                d = n_d;
                r = Some(v);
            }
        }

        (r.unwrap(), d)
    }

    #[allow(clippy::same_item_push)] // Clippy is stupid again
    pub fn each_subset(&self, size: usize, callback: impl Fn(Self))
    where
        A: Clone + Debug,
    {
        let mut mask = Vec::with_capacity(self.words.len());
        for _ in 0..size {
            mask.push(true);
        }
        for _ in size..self.words.len() {
            mask.push(false);
        }

        let mut i = 0;
        super::permut::permut(mask, |set| {
            let mut subset = Vec::with_capacity(size);
            for (i, item) in set.iter().enumerate() {
                if *item {
                    subset.push(self.words[i].clone());
                }
            }

            i += 1;
            if i % 1000 == 0 {
                println!("{i}: {:?}", subset)
            }

            callback(Code { words: subset })
        })
    }
}

impl<A: Display> Display for Code<A> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Code {{")?;
        for word in self.words.iter() {
            write!(f, "{},", word)?;
        }
        write!(f, "}}")
    }
}

impl<A> IntoIterator for Code<A> {
    type Item = Word<A>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.words.into_iter()
    }
}

///
/// The hamming bound for a code over an
/// alphabeth with r chacracters with codewords length n
/// and an error correction of t bits.
///
pub fn hamming_bound(r: usize, n: usize, t: usize) -> f64 {
    r.pow(n as u32) as f64
        / (0..=t)
            .into_iter()
            .fold(0, |sum, k| sum + binom(n, k) * (r - 1).pow(k as u32)) as f64
}

///
/// The singelton bound for a code if r charcters with codewords
/// length n that should have a codeword distance greater ir equal
/// than d.
pub fn singelton_bound(r: usize, n: usize, d: usize) -> usize {
    r.pow((n - d + 1) as u32)
}
