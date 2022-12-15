use std::{
    fmt::Debug,
    ops::{Add, AddAssign, Index, IndexMut, Mul, MulAssign, Sub, SubAssign},
};

mod space;
use num_traits::{One, Zero};
pub use space::*;

use crate::prelude::{Matrix, MatrixLayout};

#[derive(Clone, Hash)]
pub struct Vector<K> {
    pub(crate) data: Vec<K>,
    pub(crate) orientation: VectorOrientation,
}

impl<K> Vector<K> {
    pub fn orientation(mut self, orientation: VectorOrientation) -> Self {
        self.orientation = orientation;
        self
    }

    pub fn new<const N: usize>(data: [K; N]) -> Self {
        Self {
            data: data.into(),
            orientation: VectorOrientation::Undefined,
        }
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }

    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    pub fn rowvec<const N: usize>(data: [K; N]) -> Self {
        Self {
            data: data.into(),
            orientation: VectorOrientation::Row,
        }
    }

    pub fn colvec<const N: usize>(data: [K; N]) -> Self {
        Self {
            data: data.into(),
            orientation: VectorOrientation::Colum,
        }
    }

    pub fn from_vec(data: Vec<K>, orientation: VectorOrientation) -> Self {
        Self { data, orientation }
    }
}

impl<K: Copy + Zero + One + PartialEq> Vector<K> {
    pub fn unit<const N: usize, const I: usize>() -> Self {
        assert!(I < N,);

        let mut data = [K::zero(); N];
        data[I] = K::one();
        Self::new(data)
    }

    pub fn is_unit(&self) -> bool {
        let mut f = false;
        for i in 0..self.data.len() {
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

impl<K> Index<usize> for Vector<K> {
    type Output = K;
    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<K> IndexMut<usize> for Vector<K> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}

impl<K: Debug> Debug for Vector<K> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.data.iter()).finish()
    }
}

// ADD

impl<K> Add<Vector<K>> for Vector<K>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn add(mut self, rhs: Vector<K>) -> Self::Output {
        for i in 0..self.data.len() {
            self[i] = self[i].clone() + rhs[i].clone()
        }
        self
    }
}

impl<K> Add<&'_ Vector<K>> for Vector<K>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn add(mut self, rhs: &'_ Vector<K>) -> Self::Output {
        for i in 0..self.data.len() {
            self[i] = self[i].clone() + rhs[i].clone()
        }
        self
    }
}

impl<K> Add<Vector<K>> for &'_ Vector<K>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn add(self, rhs: Vector<K>) -> Self::Output {
        self.clone() + rhs
    }
}

impl<K> Add<&'_ Vector<K>> for &'_ Vector<K>
where
    K: Add<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn add(self, rhs: &'_ Vector<K>) -> Self::Output {
        self.clone() + rhs
    }
}

// ADD ASSIGN

impl<K> AddAssign<Vector<K>> for Vector<K>
where
    K: AddAssign + Clone,
{
    fn add_assign(&mut self, rhs: Vector<K>) {
        for i in 0..self.data.len() {
            self[i] += rhs[i].clone()
        }
    }
}

impl<K> AddAssign<&'_ Vector<K>> for Vector<K>
where
    K: AddAssign + Clone,
{
    fn add_assign(&mut self, rhs: &'_ Vector<K>) {
        for i in 0..self.data.len() {
            self[i] += rhs[i].clone()
        }
    }
}

// SUB

impl<K> Sub<Vector<K>> for Vector<K>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn sub(mut self, rhs: Vector<K>) -> Self::Output {
        for i in 0..self.data.len() {
            self[i] = self[i].clone() - rhs[i].clone()
        }
        self
    }
}

impl<K> Sub<&'_ Vector<K>> for Vector<K>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn sub(mut self, rhs: &'_ Vector<K>) -> Self::Output {
        for i in 0..self.data.len() {
            self[i] = self[i].clone() - rhs[i].clone()
        }
        self
    }
}

impl<K> Sub<Vector<K>> for &'_ Vector<K>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn sub(self, rhs: Vector<K>) -> Self::Output {
        self.clone() - rhs
    }
}

impl<K> Sub<&'_ Vector<K>> for &'_ Vector<K>
where
    K: Sub<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn sub(self, rhs: &'_ Vector<K>) -> Self::Output {
        self.clone() - rhs
    }
}

// SUB ASSIGN

impl<K> SubAssign<Vector<K>> for Vector<K>
where
    K: SubAssign + Clone,
{
    fn sub_assign(&mut self, rhs: Vector<K>) {
        for i in 0..self.data.len() {
            self[i] -= rhs[i].clone()
        }
    }
}

impl<K> SubAssign<&'_ Vector<K>> for Vector<K>
where
    K: SubAssign + Clone,
{
    fn sub_assign(&mut self, rhs: &'_ Vector<K>) {
        for i in 0..self.data.len() {
            self[i] -= rhs[i].clone()
        }
    }
}

// MUL (SCALAR)

impl<K> Mul<K> for Vector<K>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn mul(mut self, rhs: K) -> Self::Output {
        for i in 0..self.data.len() {
            self[i] = self[i].clone() * rhs.clone();
        }
        self
    }
}

impl<K> Mul<&'_ K> for Vector<K>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn mul(mut self, rhs: &'_ K) -> Self::Output {
        for i in 0..self.data.len() {
            self[i] = self[i].clone() * rhs.clone();
        }
        self
    }
}

impl<K> Mul<K> for &'_ Vector<K>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn mul(self, rhs: K) -> Self::Output {
        self.clone() * rhs
    }
}

impl<K> Mul<&'_ K> for &'_ Vector<K>
where
    K: Mul<K, Output = K> + Clone,
{
    type Output = Vector<K>;
    fn mul(self, rhs: &'_ K) -> Self::Output {
        self.clone() * rhs
    }
}

impl<K> MulAssign<K> for Vector<K>
where
    K: MulAssign<K> + Clone,
{
    fn mul_assign(&mut self, rhs: K) {
        for i in 0..self.data.len() {
            self[i] *= rhs.clone()
        }
    }
}
impl<K> MulAssign<&'_ K> for Vector<K>
where
    K: MulAssign<K> + Clone,
{
    fn mul_assign(&mut self, rhs: &'_ K) {
        for i in 0..self.data.len() {
            self[i] *= rhs.clone()
        }
    }
}

// FROM (MATRIX)

impl<K> From<Matrix<K>> for Vector<K> {
    fn from(vector: Matrix<K>) -> Self {
        assert!(vector.layout().is_vec());
        if vector.layout().is_colvec() {
            Vector::from_vec(vector.into_raw(), VectorOrientation::Colum)
        } else {
            Vector::from_vec(vector.into_raw(), VectorOrientation::Row)
        }
    }
}

impl<K> From<Vector<K>> for Matrix<K> {
    fn from(vector: Vector<K>) -> Self {
        Matrix::new(
            match vector.orientation {
                VectorOrientation::Colum => MatrixLayout::new(vector.data.len(), 1),
                VectorOrientation::Row => MatrixLayout::new(1, vector.data.len()),
                VectorOrientation::Undefined => MatrixLayout::new(1, vector.data.len()),
            },
            Vec::from(vector.data),
        )
    }
}

impl<K: PartialEq> PartialEq for Vector<K> {
    fn eq(&self, other: &Self) -> bool {
        self.data == other.data
    }
}

impl<K: Eq> Eq for Vector<K> {}
