use std::ops::Mul;

///
/// A byte distribution.
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct FullSourceDistr {
    distr: [usize; 256],
    len: usize,
}

impl FullSourceDistr {
    ///
    /// The number of bytes that were sampled.
    ///
    pub fn len(&self) -> usize {
        self.len
    }

    ///
    /// The entropy base r of the source.
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
    /// The alphabeth of used bytes.
    ///
    pub fn alph(&self) -> impl Iterator<Item = usize> + '_ {
        self.distr
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| if *e != 0 { Some(idx) } else { None })
    }

    ///
    /// The alphabeth of used bytes with their respective occurences.
    ///
    pub fn alph_with_p(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.distr
            .iter()
            .enumerate()
            .filter_map(|(idx, e)| if *e != 0 { Some((idx, *e)) } else { None })
    }

    ///
    /// Creates a byte distribution by using the given raw counts of each byte.
    ///
    /// Note that the input should not be greater than 256. If it is smaller all higher
    /// bytes will be asuumed to have an occurence count of 0.
    ///
    pub fn from_distr(probs: &[usize]) -> Self {
        assert!(probs.len() <= 256);
        let mut this = Self {
            distr: [0; 256],
            len: 0,
        };
        for i in 0..probs.len() {
            this.distr[i] += probs[i];
            this.len += probs[i];
        }
        this
    }

    ///
    /// Creates a new empty distribution. If full_sample is false each character is assumed to
    /// occure at least once so that the code build from this distribution knows it may has
    /// to expect a byte that will not be sampled.
    ///
    pub fn new(full_sample: bool) -> Self {
        if full_sample {
            Self {
                distr: [0; 256],
                len: 0,
            }
        } else {
            Self {
                distr: [1; 256],
                len: 256,
            }
        }
    }

    ///
    /// Samples a bytesstream.
    ///
    pub fn add_bytes(&mut self, bytes: &[u8]) {
        for byte in bytes {
            self.add_byte(*byte)
        }
    }

    ///
    /// Adds a byte to the sample.
    ///
    #[inline]
    pub fn add_byte(&mut self, byte: u8) {
        self.distr[byte as usize] += 1;
        self.len += 1;
    }

    ///
    /// Generates a new source by combining `lhs` and `rhs` into a product source.
    ///  
    pub fn product_source(lhs: &FullSourceDistr, rhs: &FullSourceDistr) -> Self {
        // Assumes that only the first n values are given
        let n1 = lhs
            .distr
            .iter()
            .enumerate()
            .find(|(_, e)| **e == 0)
            .unwrap_or((256, &usize::MAX))
            .0;
        let n2 = rhs
            .distr
            .iter()
            .enumerate()
            .find(|(_, e)| **e == 0)
            .unwrap_or((256, &usize::MAX))
            .0;

        assert!(n1 * n2 <= 256);

        let mut new_dist = [0; 256];
        let mut k = 0;
        let mut len = 0;
        for i in 0..n1 {
            for j in 0..n2 {
                let occur = lhs.distr[i] * rhs.distr[j];
                new_dist[k] = occur;
                len += occur;
                k += 1;
            }
        }

        assert!(k == n1 * n2);

        Self {
            distr: new_dist,
            len,
        }
    }
}

impl Mul for FullSourceDistr {
    type Output = FullSourceDistr;
    fn mul(self, rhs: Self) -> Self::Output {
        FullSourceDistr::product_source(&self, &rhs)
    }
}

impl Mul<&FullSourceDistr> for FullSourceDistr {
    type Output = FullSourceDistr;
    fn mul(self, rhs: &FullSourceDistr) -> Self::Output {
        FullSourceDistr::product_source(&self, rhs)
    }
}

impl Mul for &FullSourceDistr {
    type Output = FullSourceDistr;
    fn mul(self, rhs: Self) -> Self::Output {
        FullSourceDistr::product_source(self, rhs)
    }
}

impl Mul<FullSourceDistr> for &FullSourceDistr {
    type Output = FullSourceDistr;
    fn mul(self, rhs: FullSourceDistr) -> Self::Output {
        FullSourceDistr::product_source(self, &rhs)
    }
}
