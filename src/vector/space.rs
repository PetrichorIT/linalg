use super::Vector;
use crate::sets::{CField, Field, Finite};
use std::marker::PhantomData;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VectorOrientation {
    Row,
    Colum,
    Undefined,
}

pub struct VectorSpace<K: Field, const N: usize> {
    _phantom: PhantomData<K>,
}

impl<K: CField, const N: usize> VectorSpace<K, N> {
    pub fn unit_vectors() -> [Vector<K, N>; N] {
        let mut zero = [Vector::new([K::zero(); N]); N];
        for (i, item) in zero.iter_mut().enumerate() {
            item[i] = K::one()
        }
        zero
    }
}

impl<K: CField + Finite, const N: usize> VectorSpace<K, N> {
    pub fn all_vectors() -> Vec<Vector<K, N>> {
        let mut result = Vec::new();
        let mut buf = Vector::new([K::zero(); N]);

        fn permut<K: CField, const N: usize>(
            i: usize,
            rem: &[K],
            result: &mut Vec<Vector<K, N>>,
            buf: &mut Vector<K, N>,
        ) {
            if i >= N {
                result.push(*buf)
            } else {
                for k in 0..rem.len() {
                    let e = rem[k];
                    buf[i] = e;
                    permut(i + 1, rem, result, buf);
                }
            }
        }

        permut(0, &K::all(), &mut result, &mut buf);

        result
    }

    pub fn equivalence_classes() -> Vec<Vector<K, N>> {
        let scalars = K::all();
        let mut all = Self::all_vectors();
        let mut classes = Vec::new();

        // start with unity vectos -- they cant be in the same class.
        let mut i = 0;
        while i < all.len() {
            if all[i].is_unit() {
                let unit = all.swap_remove(i);
                // Generate class members
                all.retain(|e| !scalars.iter().any(|s| unit * s == *e));
                classes.push(unit);
            } else {
                i += 1;
            }
        }

        while !all.is_empty() {
            let next = all.swap_remove(0);
            all.retain(|e| !scalars.iter().any(|s| next * s == *e));
            classes.push(next);
        }

        classes
    }
}
