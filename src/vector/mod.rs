use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Index, IndexMut, Mul, Sub, SubAssign},
};

mod space;
use num_traits::{One, Zero};
pub use space::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Vector<K, const N: usize> {
    data: [K; N],
    orientation: VectorOrientation,
}

impl<K, const N: usize> Vector<K, N> {
    pub const fn new(data: [K; N]) -> Self {
        Self {
            data,
            orientation: VectorOrientation::Undefined,
        }
    }

    pub const fn rowvec(data: [K; N]) -> Self {
        Self {
            data,
            orientation: VectorOrientation::Row,
        }
    }

    pub const fn colvec(data: [K; N]) -> Self {
        Self {
            data,
            orientation: VectorOrientation::Colum,
        }
    }
}

impl<K: Copy + Zero + One + PartialEq, const N: usize> Vector<K, N> {
    pub fn unit<const I: usize>() -> Self {
        assert!(I < N,);

        let mut data = [K::zero(); N];
        data[I] = K::one();
        Self::new(data)
    }

    pub fn is_unit(&self) -> bool {
        let mut f = false;
        for i in 0..N {
            if !self[i].is_zero() {
                if self[i].is_one() {
                    if f {
                        return false;
                    } else {
                        f = true;
                    }
                } else {
                    return false;
                }
            }
        }
        f
    }
}

impl<K, const N: usize> Index<usize> for Vector<K, N> {
    type Output = K;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<K, const N: usize> IndexMut<usize> for Vector<K, N> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

// ADD

impl<K, const N: usize> Add<Vector<K, N>> for Vector<K, N>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn add(mut self, rhs: Vector<K, N>) -> Self::Output {
        for i in 0..N {
            self[i] = self[i].clone() + rhs[i].clone()
        }
        self
    }
}

impl<K, const N: usize> Add<&'_ Vector<K, N>> for Vector<K, N>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn add(mut self, rhs: &'_ Vector<K, N>) -> Self::Output {
        for i in 0..N {
            self[i] = self[i].clone() + rhs[i].clone()
        }
        self
    }
}

impl<K, const N: usize> Add<Vector<K, N>> for &'_ Vector<K, N>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn add(self, rhs: Vector<K, N>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<K, const N: usize> Add<&'_ Vector<K, N>> for &'_ Vector<K, N>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn add(self, rhs: &'_ Vector<K, N>) -> Self::Output {
        self.clone() + rhs
    }
}

// ADD ASSIGN

impl<K, const N: usize> AddAssign<Vector<K, N>> for Vector<K, N>
where
    K: AddAssign + Clone,
{
    fn add_assign(&mut self, rhs: Vector<K, N>) {
        for i in 0..N {
            self[i] += rhs[i].clone()
        }
    }
}

impl<K, const N: usize> AddAssign<&'_ Vector<K, N>> for Vector<K, N>
where
    K: AddAssign + Clone,
{
    fn add_assign(&mut self, rhs: &'_ Vector<K, N>) {
        for i in 0..N {
            self[i] += rhs[i].clone()
        }
    }
}

// SUB

impl<K, const N: usize> Sub<Vector<K, N>> for Vector<K, N>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn sub(mut self, rhs: Vector<K, N>) -> Self::Output {
        for i in 0..N {
            self[i] = self[i].clone() - rhs[i].clone()
        }
        self
    }
}

impl<K, const N: usize> Sub<&'_ Vector<K, N>> for Vector<K, N>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn sub(mut self, rhs: &'_ Vector<K, N>) -> Self::Output {
        for i in 0..N {
            self[i] = self[i].clone() - rhs[i].clone()
        }
        self
    }
}

impl<K, const N: usize> Sub<Vector<K, N>> for &'_ Vector<K, N>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn sub(self, rhs: Vector<K, N>) -> Self::Output {
        self.clone() - rhs
    }
}

impl<K, const N: usize> Sub<&'_ Vector<K, N>> for &'_ Vector<K, N>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn sub(self, rhs: &'_ Vector<K, N>) -> Self::Output {
        self.clone() - rhs
    }
}

// SUB ASSIGN

impl<K, const N: usize> SubAssign<Vector<K, N>> for Vector<K, N>
where
    K: SubAssign + Clone,
{
    fn sub_assign(&mut self, rhs: Vector<K, N>) {
        for i in 0..N {
            self[i] -= rhs[i].clone()
        }
    }
}

impl<K, const N: usize> SubAssign<&'_ Vector<K, N>> for Vector<K, N>
where
    K: SubAssign + Clone,
{
    fn sub_assign(&mut self, rhs: &'_ Vector<K, N>) {
        for i in 0..N {
            self[i] -= rhs[i].clone()
        }
    }
}

// MUL (SCALAR)

impl<K, const N: usize> Mul<K> for Vector<K, N>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn mul(mut self, rhs: K) -> Self::Output {
        for i in 0..N {
            self[i] = self[i].clone() * rhs.clone();
        }
        self
    }
}

impl<K, const N: usize> Mul<&'_ K> for Vector<K, N>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn mul(mut self, rhs: &'_ K) -> Self::Output {
        for i in 0..N {
            self[i] = self[i].clone() * rhs.clone();
        }
        self
    }
}

impl<K, const N: usize> Mul<K> for &'_ Vector<K, N>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn mul(self, rhs: K) -> Self::Output {
        self.clone() * rhs
    }
}

impl<K, const N: usize> Mul<&'_ K> for &'_ Vector<K, N>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K, N>;
    fn mul(self, rhs: &'_ K) -> Self::Output {
        self.clone() * rhs
    }
}
