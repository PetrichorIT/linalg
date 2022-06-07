//! Definitions for sources in the context of source encoding.

///
/// A generic distribution of
/// from a generic alphabeth.
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GenericSource<T> {
    alphabeth: Vec<T>,
    distr: Vec<usize>,
    len: usize,
}

#[allow(unused)]
impl<T> GenericSource<T> {
    ///
    /// The number of chracters sampled.
    ///
    pub fn len(&self) -> usize {
        self.len
    }

    ///
    /// Whether the code is empty.
    ///
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    ///
    /// The entropy of the source base r.
    ///
    pub fn entropy(&self, r: f64) -> f64 {
        let mut sum = 0.0;
        for p in self.distr.iter() {
            if *p != 0 {
                let p = *p as f64 / self.len as f64;
                sum -= p * p.log(r);
            }
        }
        sum
    }

    ///
    /// All relevant characters with an occurence probability > 0.
    ///
    pub fn alph(&self) -> impl Iterator<Item = &T> {
        self.alphabeth.iter()
    }

    ///
    /// All relevant characters with an occurence probability > 0 and their absoulute occurence counts.
    ///
    pub fn alph_with_p(&self) -> impl Iterator<Item = (&T, usize)> {
        self.alphabeth.iter().zip(self.distr.iter().copied())
    }

    ///
    /// Generates a distribution from a alphabeth and manual occurence counts.
    ///
    pub fn from_distr<F>(distr: Vec<usize>, items: F) -> Self
    where
        F: Fn(usize) -> T,
    {
        let alphabeth = (0..distr.len()).map(items).collect();
        let len = distr.iter().sum();
        Self {
            alphabeth,
            distr,
            len,
        }
    }

    ///
    /// Manually creates a new distribution.
    ///
    pub fn new(alphabeth: Vec<T>, distr: Vec<usize>) -> Self {
        assert_eq!(alphabeth.len(), distr.len());
        let len = distr.iter().sum();

        Self {
            alphabeth,
            distr,
            len,
        }
    }

    ///
    /// Generates the product source to self*self.
    ///
    pub fn product_source(&self) -> GenericSource<(T, T)>
    where
        T: Clone,
    {
        // Assumes that only the first n values are given
        let n = self
            .distr
            .iter()
            .enumerate()
            .find(|(_, e)| **e == 0)
            .unwrap()
            .0;
        assert!(n * n <= 256);

        let mut alphabeth = Vec::new();
        let mut new_dist = Vec::new();
        let mut len = 0;
        for i in 0..n {
            for j in 0..n {
                let occur = self.distr[i] * self.distr[j];
                alphabeth.push((self.alphabeth[i].clone(), self.alphabeth[j].clone()));
                new_dist.push(occur);
                len += occur;
            }
        }

        GenericSource {
            alphabeth,
            distr: new_dist,
            len,
        }
    }
}
