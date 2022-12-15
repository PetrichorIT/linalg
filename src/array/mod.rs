use std::ops::{Index, IndexMut};

pub struct Array<T> {
    data: Vec<T>,
    shape: Vec<usize>,
}

impl<T: Copy> Array<T> {
    pub fn fill<const N: usize>(shape: [usize; N], filler: T) -> Self {
        Self {
            data: vec![filler; Self::len_of_shape(&shape)],
            shape: shape.to_vec(),
        }
    }
}

impl<T> Array<T> {
    fn len_of_shape(shape: &[usize]) -> usize {
        let mut r = 1;
        for dim in shape {
            r *= dim;
        }
        r
    }

    fn offset_for_index<const N: usize>(&self, index: [usize; N]) -> usize {
        assert_eq!(N, self.shape.len());
        index
            .iter()
            .zip(&self.shape)
            .map(|(idx, dim_size)| idx * dim_size)
            .sum()
    }
}

impl<T, const N: usize> Index<[usize; N]> for Array<T> {
    type Output = T;
    fn index(&self, index: [usize; N]) -> &Self::Output {
        let offset = self.offset_for_index(index);
        &self.data[offset]
    }
}

impl<T, const N: usize> IndexMut<[usize; N]> for Array<T> {
    fn index_mut(&mut self, index: [usize; N]) -> &mut Self::Output {
        let offset = self.offset_for_index(index);
        &mut self.data[offset]
    }
}
