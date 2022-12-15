use crate::{
    algebra::_gauss_jordan_impl,
    prelude::{Matrix, GF},
    sets::Finite,
    vector::{Vector, VectorSpace},
};

use num_traits::Zero;

pub struct LinearCode<const Q: usize, const N: usize, const K: usize> {
    def: LinearCodeDef<Q>,
}

pub enum LinearCodeDef<const Q: usize> {
    GeneratorMatrix(Matrix<GF<Q>>),
    ControllMatrix(Matrix<GF<Q>>),
}

impl<const Q: usize, const N: usize, const K: usize> LinearCode<Q, N, K> {
    pub fn controll_matrix(&self) -> Option<&Matrix<GF<Q>>> {
        match self.def {
            LinearCodeDef::ControllMatrix(ref m) => Some(m),
            _ => None,
        }
    }

    pub fn generator_matrix(&self) -> Option<&Matrix<GF<Q>>> {
        match self.def {
            LinearCodeDef::GeneratorMatrix(ref m) => Some(m),
            _ => None,
        }
    }

    pub fn distance(&self) -> usize {
        self.all_codewords()
            .into_iter()
            .filter(|w| !w.data.iter().all(|v| *v == GF(0)))
            .map(|w| {
                w.data
                    .into_iter()
                    .map(|w| w.is_zero())
                    .filter(|b| !*b)
                    .count()
            })
            .min()
            .unwrap()
    }

    pub fn with_generator(generator: Matrix<GF<Q>>) -> Self {
        assert_eq!(generator.layout().cols(), N);
        assert_eq!(generator.layout().rows(), K);
        Self {
            def: LinearCodeDef::GeneratorMatrix(generator),
        }
    }

    pub fn with_control(control: Matrix<GF<Q>>) -> Self {
        assert_eq!(control.layout().rows(), N);
        Self {
            def: LinearCodeDef::ControllMatrix(control),
        }
    }

    pub fn hamming() -> Option<Self> {
        // Q, N, K allready computed
        let r = N - K;

        // Check Hamming 1950
        if N != (Q.pow(r as u32) - 1) / (Q - 1) {
            println!(
                "Hamming 1950 error: N = {} K = {} R = {} Q = {}",
                N, K, r, Q
            );
            return None;
        }

        let vr = VectorSpace::<GF<Q>>::new(r).equivalence_classes();

        let mut h = Matrix::fill((vr.len(), r), GF(0));
        for (i, v) in vr.into_iter().enumerate() {
            h.insert_vector(
                i..,
                0..,
                v.orientation(crate::vector::VectorOrientation::Row),
            )
        }

        Some(Self {
            def: LinearCodeDef::ControllMatrix(h),
        })
    }

    pub fn reed_solomon(d: usize) -> Option<Self> {
        if d >= Q || Q < 2 {
            return None;
        }

        let n = Q + 1;
        let mut h = Matrix::fill((n, d - 1), GF(0));

        // (1) Set first row (zeros allready set)
        h[(0, d - 2)] = GF(1);

        // (2) Set first col
        for i in 1..n {
            h[(i, 0)] = GF(1);
        }

        // (3) Create partial vandermonde
        for (r, g) in GF::<Q>::all().into_iter().enumerate() {
            let mut acc = g;
            h[(r, 1)] = acc;
            for c in 2..(d - 1) {
                acc *= g;
                h[(r, c)] = acc;
            }
        }

        Some(Self {
            def: LinearCodeDef::ControllMatrix(h),
        })
    }

    pub fn all_codewords(&self) -> Vec<Vector<GF<Q>>> {
        match self.def {
            LinearCodeDef::GeneratorMatrix(_) => VectorSpace::<GF<Q>>::new(K)
                .all_vectors()
                .into_iter()
                .map(|v| self.encode(v).unwrap())
                .collect(),
            LinearCodeDef::ControllMatrix(ref h) => VectorSpace::<GF<Q>>::new(N)
                .all_vectors()
                .into_iter()
                .filter(|w| (Matrix::from(w.clone()) * h).is_zero())
                .collect(),
        }
    }

    pub fn encode(&self, word: Vector<GF<Q>>) -> Option<Vector<GF<Q>>> {
        assert_eq!(word.len(), K);

        let word = Matrix::from(word);
        let result = match self.def {
            LinearCodeDef::GeneratorMatrix(ref g) => Matrix::mmul(&word, g),
            _ => return None,
        };

        Some(Vector::from(result))
    }
}

impl<const Q: usize, const N: usize, const K: usize> PartialEq for LinearCode<Q, N, K> {
    fn eq(&self, other: &Self) -> bool {
        let lhs = self.all_codewords();
        let mut rhs = other.all_codewords();

        'outer: for word in lhs {
            for i in 0..rhs.len() {
                if rhs[i] == word {
                    rhs.swap_remove(i);
                    continue 'outer;
                }
            }

            return false;
        }

        return rhs.is_empty();
    }
}

pub fn generator_to_controll<const Q: usize>(mut generator: Matrix<GF<Q>>) -> Matrix<GF<Q>> {
    let k = generator.layout().rows();
    let n = generator.layout().cols();

    _gauss_jordan_impl(&mut generator).expect("Failed gauss jordan. Code was not linear");
    let p = generator.extract(0.., k..);
    assert_eq!(p.layout(), (k, n - k).into());

    // -P^T
    let mut h = Matrix::fill((n, n - k), GF(0));
    h.insert(0.., 0.., &p);
    h.insert(n - k.., 0.., &Matrix::eye(k));

    h
}
