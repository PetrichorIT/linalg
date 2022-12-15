use super::Vector;
use crate::sets::{CField, Field, Finite};
use std::marker::PhantomData;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum VectorOrientation {
    Row,
    Colum,
    Undefined,
}

pub struct VectorSpace<K: Field> {
    n: usize,
    _phantom: PhantomData<K>,
}

impl<K: CField> VectorSpace<K> {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            _phantom: PhantomData,
        }
    }

    pub fn unit_vectors(&self) -> Vec<Vector<K>> {
        let mut zero =
            vec![Vector::from_vec(vec![K::zero(); self.n], VectorOrientation::Undefined); self.n];
        for (i, item) in zero.iter_mut().enumerate() {
            item[i] = K::one()
        }
        zero
    }
}

impl<K: CField + Finite> VectorSpace<K> {
    pub fn all_vectors(&self) -> Vec<Vector<K>> {
        let mut result = Vec::new();
        let mut buf = Vector::from_vec(vec![K::zero(); self.n], VectorOrientation::Undefined);

        fn permut<K: CField>(
            i: usize,
            rem: &[K],
            result: &mut Vec<Vector<K>>,
            buf: &mut Vector<K>,
        ) {
            if i >= buf.len() {
                result.push(buf.clone())
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

    pub fn equivalence_classes(&self) -> Vec<Vector<K>> {
        let scalars = K::all();
        let mut all = self.all_vectors();
        let mut classes = Vec::new();

        all.reverse();

        // start with unity vectos -- they cant be in the same class.
        let mut i = 0;
        while i < all.len() {
            if all[i].is_unit() {
                let unit = all.remove(i);
                // Generate class members
                all.retain(|e| !scalars.iter().any(|s| &unit * s == *e));
                classes.push(unit);
            } else {
                i += 1;
            }
        }

        while !all.is_empty() {
            let next = all.remove(0);
            all.retain(|e| !scalars.iter().any(|s| &next * s == *e));
            classes.push(next);
        }

        classes
    }
}
