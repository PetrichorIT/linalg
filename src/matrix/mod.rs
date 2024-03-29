//!
//! A module defining the core structures of [linalg](crate).
//!
//! # Matrices
//!
//! A [Matrix] is defined through a raw buffer and a [Matrix Layout](MatrixLayout)
//! defining the cell interpretation of the buffer.
//! Generic functions and mathematical operations are defined on the generic
//! type a piori, more complex operations will be added in other modules.
//!

use std::{
    convert::TryFrom,
    fmt::{Debug, Display},
    hash::Hash,
    mem::swap,
    ops::{
        Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Deref,
        DerefMut, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, RangeInclusive, Sub,
        SubAssign,
    },
};

use num_traits::{Num, Zero};

use crate::{num::ClippableRange, vector::Vector};

mod tests;

///
/// A description of a 2d-matrix dimensions.
///
/// This type represents the layout of a two-dimensional matrix
/// defined by rows and collums. This type should be used in conjunction
/// with the generic `Matrix` struct.
///
/// # Example
///
/// ```
/// use linalg::prelude::*;
///
/// let layout = MatrixLayout::new(2, 3);
/// let matrix = Matrix::fill(layout, 0usize);
/// assert!(matrix.layout().rows() == 2 && matrix.layout().cols() == 3);
/// assert!(matrix[(0, 0)] == 0usize);
/// ```
///
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub struct MatrixLayout {
    rows: usize,
    cols: usize,
}

impl MatrixLayout {
    /// The number of rows in the matrix.
    #[inline(always)]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// The number of collums in the matrix.
    #[inline(always)]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// The total number of cells in the matrix.
    #[inline(always)]
    pub fn size(&self) -> usize {
        self.rows * self.cols
    }

    /// Indicator if the layout matches a square matrix layout.
    #[inline(always)]
    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }

    /// Indicator if the layout macthes a row vector.
    #[inline(always)]
    pub fn is_rowvec(&self) -> bool {
        self.rows == 1
    }

    /// Indicator if the layout macthes a collum vector.
    #[inline(always)]
    pub fn is_colvec(&self) -> bool {
        self.cols == 1
    }

    /// Indicator if the layout is either a rowvec or colvec.
    #[inline(always)]
    pub fn is_vec(&self) -> bool {
        self.cols == 1 || self.rows == 1
    }

    #[inline]
    pub fn bounds(&self) -> (RangeInclusive<usize>, RangeInclusive<usize>) {
        (0..=(self.rows - 1), 0..=(self.cols - 1))
    }
    ///
    /// A conversion function from a two-dimensional index to a
    /// one-dimensional index for the raw buffer.
    /// This conversion follows the following formular:
    /// BufferIndex = Index_Row * Collum_Size + Index_Collum
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let layout = MatrixLayout::new(2, 3);
    /// assert!(layout.index((1, 1)) == 4);
    /// ```
    ///
    #[inline(always)]
    pub fn index(&self, index: (usize, usize)) -> usize {
        index.0 * self.cols + index.1
    }

    /// Performs an inplace transposition of the layout.
    pub fn transpose(&mut self) {
        swap(&mut self.rows, &mut self.cols)
    }

    /// Performs an out-of-place transposition of the layout and returns the result.
    pub fn transposed(&self) -> Self {
        MatrixLayout {
            rows: self.cols,
            cols: self.rows,
        }
    }
}

impl MatrixLayout {
    /// Creates a new matrix layout.
    #[inline]
    pub fn new(rows: usize, cols: usize) -> Self {
        Self { rows, cols }
    }

    ///
    /// Tries to generate a layout based on a given two-dimensonal vector of elements.
    /// This layout will try to match the given dataset.
    /// This conversion will fail if the stuctured data has an irregular format,
    /// meaning there a rows with a different number of elements in them.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let data = vec![
    ///     vec![1, 2, 3],
    ///     vec![4, 5, 6],
    /// ];
    ///
    /// let layout = MatrixLayout::assumed_from(&data);
    /// assert!(layout.is_ok());
    /// assert!(layout.unwrap() == MatrixLayout::new(2, 3));
    /// ```
    ///
    pub fn assumed_from<T>(structured: &[Vec<T>]) -> Result<Self, &'static str> {
        let rows = structured.len();
        if rows == 0 {
            return Ok(MatrixLayout { rows: 0, cols: 0 });
        }

        let cols = structured.first().unwrap().len();
        for row in structured.iter().skip(1) {
            if row.len() != cols {
                return Err("Structured data has irregulare format");
            }
        }

        Ok(MatrixLayout { rows, cols })
    }

    /// Generates a square matrix layout, given the edge size.
    #[inline]
    pub fn square(size: usize) -> Self {
        Self {
            rows: size,
            cols: size,
        }
    }
}

impl From<(usize, usize)> for MatrixLayout {
    fn from(tupel: (usize, usize)) -> Self {
        Self {
            rows: tupel.0,
            cols: tupel.1,
        }
    }
}

///
/// A two-dimensional matrix of generic elements.
///
/// This struct contains a layout element of type [MatrixLayout] and
/// a raw buffer to store matrix cells. The size of the row-major raw buffer thereby is
/// defined by the layouts [MatrixLayout::size()] function.
/// Note that the [Matrix] has technicly no generic constraint but most functions are
/// only implemented over the [Num] trait.
///
/// # Example
///
/// ```
/// use linalg::prelude::*;
///
/// let matrix: Matrix<usize> = Matrix::fill(MatrixLayout::new(2, 3), 0);
/// assert!(matrix.size() == 2*3);
/// assert!(matrix[(0, 0)] == 0);
/// ```
///
#[derive(Debug, Clone)]
pub struct Matrix<T> {
    layout: MatrixLayout,
    raw: Vec<T>,
}

impl<T> Matrix<T> {
    pub(crate) fn into_raw(self) -> Vec<T> {
        self.raw
    }

    /// The layout constraints definig the matrix.
    #[inline(always)]
    pub fn layout(&self) -> MatrixLayout {
        self.layout
    }

    /// The raw buffer where the cells are stored.
    #[inline(always)]
    pub fn raw(&self) -> &Vec<T> {
        &self.raw
    }

    /// The total number of cells in the raw buffer (used).
    #[inline(always)]
    pub fn size(&self) -> usize {
        self.layout.size()
    }

    ///
    /// Creates a new matrix with a given layout over the given buffer.
    ///  
    /// # Panics
    ///
    /// This function panics should the layout not fit the buffer.
    /// This can be archived should the layouts size match the buffers length.
    ///
    /// ```should_panic
    /// use linalg::prelude::*;
    ///
    /// // Creates a 6 element buffer
    /// // -> layouts (1, 6) (2, 3) (3, 2) (6, 1) are possible
    /// let buffer = vec![1, 2, 3, 4, 5, 6];
    ///
    /// let matrix = Matrix::new((3, 3), buffer);
    /// ```
    ///
    #[inline]
    pub fn new<U>(layout: U, raw: Vec<T>) -> Self
    where
        U: Into<MatrixLayout>,
    {
        let layout = layout.into();
        assert_eq!(layout.size(), raw.len());
        Self { layout, raw }
    }

    /// Creates a row vector with the given buffer.
    #[inline]
    pub fn rowvec(raw: Vec<T>) -> Self {
        let layout = MatrixLayout::new(1, raw.len());
        Self { layout, raw }
    }

    /// Create a collum vector with the given buffer.
    #[inline]
    pub fn colvec(raw: Vec<T>) -> Self {
        let layout = MatrixLayout::new(raw.len(), 1);
        Self { layout, raw }
    }

    /// Creates a new matrix by combining n columvectors of size m into
    /// a matrix with m rows and n colums
    ///
    pub fn cloned_from_parts_vertical(colvecs: &[Matrix<T>]) -> Self
    where
        T: Clone,
    {
        assert!(!colvecs.is_empty());
        // assert!(colvecs[0].layout().is_colvec());

        let col_depth = colvecs[0].len();
        let layout = MatrixLayout::new(col_depth, colvecs.len());
        let mut raw = Vec::with_capacity(layout.size());

        for i in 0..col_depth {
            for colvec in colvecs.iter() {
                raw.push(colvec.raw[i].clone())
            }
        }

        Self::new(layout, raw)
    }

    pub fn referenced_from_parts_vertical(colvecs: &[Matrix<T>]) -> Matrix<&T> {
        assert!(!colvecs.is_empty());
        // assert!(colvecs[0].layout().is_colvec());

        let col_depth = colvecs[0].len();
        let layout = MatrixLayout::new(col_depth, colvecs.len());
        let mut raw = Vec::with_capacity(layout.size());

        for i in 0..col_depth {
            for colvec in colvecs.iter() {
                raw.push(&colvec.raw[i])
            }
        }

        Matrix::new(layout, raw)
    }

    ///
    /// Creates a new matrix with an allocated but uninialized
    /// memory buffer.
    ///
    /// # Safety
    ///
    /// Ensure that all cells of the allocated buffer (self.raw().len())
    /// are either filled with valid instances of the given type,
    /// or never used.
    ///
    #[allow(clippy::uninit_vec)] // This is literally the meaning of this fn.
    pub unsafe fn uninitalized<U>(layout: U) -> Self
    where
        U: Into<MatrixLayout>,
    {
        let layout = layout.into();
        let mut raw = Vec::with_capacity(layout.size());

        // SAFTY:
        // This function will only be called internally in its 'unsafe' but not
        // outwardly unsafe shape, to create buffers and such.
        raw.set_len(layout.size());

        Self { layout, raw }
    }

    ///
    /// Resets all allocated memory to zero bytes.
    ///
    /// # Safety
    ///
    /// This should only be done if zero-memory produces valid instances of type T
    /// or matrix will be reinitalized fully afterwards.
    ///
    pub unsafe fn memreset(&mut self) {
        std::ptr::write_bytes(self.raw.as_mut_ptr(), 0, self.raw.len())
    }

    ///
    /// Extacts the submatrix constainted by the ranges, and returns
    /// an out-of-place copy of its elements
    ///
    pub fn extract<U, V>(&self, rows: U, cols: V) -> Matrix<T>
    where
        T: Copy,
        U: ClippableRange<usize>,
        V: ClippableRange<usize>,
    {
        if self.is_empty() {
            return Matrix::new((0, 0), Vec::new());
        }
        let (row_bounds, col_bounds) = self.layout().bounds();
        let rows = rows.into_clipped(row_bounds);
        let cols = cols.into_clipped(col_bounds);

        let rc = rows.end() - rows.start() + 1;
        let cc = cols.end() - cols.start() + 1;

        let mut result = unsafe { Matrix::uninitalized((rc, cc)) };
        for i in 0..rc {
            for j in 0..cc {
                result[(i, j)] = self[(*rows.start() + i, *cols.start() + j)]
            }
        }

        result
    }

    pub fn insert_vector(
        &mut self,
        rows: impl ClippableRange<usize>,
        cols: impl ClippableRange<usize>,
        vector: Vector<T>,
    ) where
        T: Copy,
    {
        self.insert(rows, cols, &vector.into())
    }

    ///
    /// Inserts elements of a given matrix into the self using the provided
    /// ranges as anchor point and constrains.
    ///
    pub fn insert<U, V>(&mut self, rows: U, cols: V, matrix: &Matrix<T>)
    where
        T: Copy,
        U: ClippableRange<usize>,
        V: ClippableRange<usize>,
    {
        if matrix.len() == 0 {
            return;
        }
        let (row_bounds, col_bounds) = self.layout().bounds();

        // Bounds clipped to parent indices
        let rows = rows.into_clipped(row_bounds);
        let cols = cols.into_clipped(col_bounds);

        let anchor = (*rows.start(), *cols.start());

        // Reduce bounds to fit matrix
        let rows = rows.into_clipped(anchor.0..=(anchor.0 + matrix.layout().rows() - 1));
        let cols = cols.into_clipped(anchor.1..=(anchor.1 + matrix.layout().cols() - 1));

        let rc = rows.end() - rows.start() + 1;
        let cc = cols.end() - cols.start() + 1;

        for i in 0..rc {
            for j in 0..cc {
                self[(*rows.start() + i, *cols.start() + j)] = matrix[(i, j)]
            }
        }
    }
}

impl<T: Num> Matrix<T>
where
    T: Copy,
{
    /// Indicates whether a matrix is symetric square matrix or not.
    pub fn is_symetric(&self) -> bool {
        if !self.layout.is_square() {
            return false;
        }

        for i in 0..self.layout.rows {
            for j in 0..=i {
                if self.raw[i * self.layout.rows + j] != self.raw[j * self.layout.rows + i] {
                    return false;
                }
            }
        }

        true
    }

    ///
    /// Indicates whether all all elments over the main diagonal are zero.
    /// Note that this functions can be applied to matrices of all forms not
    /// just square matrices.
    ///
    /// # Panics
    ///
    /// This function will panic if applied to matrices of size 0.
    ///
    /// ```should_panic
    /// use linalg::{matrix, matrix::Matrix};
    ///
    /// let m: Matrix<usize> = matrix![];
    ///
    /// println!("{}", m.is_lower_triag());
    /// ````
    ///
    pub fn is_lower_triag(&self) -> bool {
        let lcr = (self.layout.cols - 1).min(self.layout.rows);

        for i in 0..lcr {
            for j in (i + 1)..self.layout.cols {
                if !self[(i, j)].is_zero() {
                    return false;
                }
            }
        }

        true
    }

    ///
    /// Indicates whether all all elments below the main diagonal are zero.
    /// Note that this functions can be applied to matrices of all forms not
    /// just square matrices.
    ///
    pub fn is_upper_triag(&self) -> bool {
        for i in 1..self.layout.rows {
            for j in 0..i.min(self.layout.cols) {
                if !self[(i, j)].is_zero() {
                    return false;
                }
            }
        }

        true
    }

    ///
    /// Indicates whether the matrix only posses elements on the main diagonal
    /// or not.
    /// Note that this functions is not limited to square matrices.
    /// This function will always return true if the matrix is created through the
    /// [Matrix::diag()] constructor.
    ///
    /// # Examples
    ///
    /// ```
    /// use linalg::matrix;
    ///
    /// let m = matrix![
    ///     1, 0, 0;
    ///     0, 2, 0;
    ///     0, 0, 3;
    ///     0, 0, 0;
    /// ];
    ///
    /// assert!(m.is_diag());
    /// ```
    ///
    pub fn is_diag(&self) -> bool {
        for i in 0..self.layout.rows {
            for j in 0..self.layout.cols {
                if !self[(i, j)].is_zero() && i != j {
                    return false;
                }
            }
        }

        true
    }

    ///
    /// Creates a new matrix copying values from the given
    /// slice into a new buffer, using the given layout.
    ///
    pub fn from_slice<U>(layout: U, slice: &[T]) -> Self
    where
        U: Into<MatrixLayout>,
    {
        let layout = layout.into();
        assert_eq!(layout.size(), slice.len());
        Self {
            layout,
            raw: slice.to_vec(),
        }
    }

    ///
    /// Creates a new matrix copying values from the given slices in order
    /// into a new buffer, using the given layout.
    ///
    pub fn from_slices<U>(layout: U, slices: Vec<&[T]>) -> Self
    where
        U: Into<MatrixLayout>,
    {
        let layout = layout.into();
        let mut buffer = Vec::with_capacity(layout.size());
        for slice in slices {
            buffer.extend_from_slice(slice);
        }
        assert_eq!(layout.size(), buffer.len());
        Self {
            layout,
            raw: buffer,
        }
    }

    ///
    /// Creates a new matrix with the given layout filling all cells
    /// with the given filler element.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::fill(MatrixLayout::new(2, 3), 1usize);
    /// assert!(matrix[(0, 0)] == 1usize);
    /// assert!(matrix[(1, 2)] == 1usize);
    /// ```
    ///
    pub fn fill<U>(layout: U, filler: T) -> Self
    where
        U: Into<MatrixLayout>,
    {
        // SAFTY:
        // Can be used since the underling vector will be filled according to its
        // own len, thus all cells used are guarnteed to be there
        // since the vector has the given size AND capacity
        let mut matrix = unsafe { Self::uninitalized(layout.into()) };
        for cell in &mut matrix.raw {
            *cell = filler;
        }
        matrix
    }

    ///
    /// Creates a new matrix with the given layout filling all cells
    /// with the given filler element.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::fill(MatrixLayout::new(2, 3), 0usize);
    /// assert!(matrix[(0, 0)] == 0usize);
    /// assert!(matrix[(1, 2)] == 0usize);
    /// ```
    ///
    #[inline]
    pub fn zeroed<U>(layout: U) -> Self
    where
        U: Into<MatrixLayout>,
    {
        Self::fill(layout, T::zero())
    }

    ///
    /// Create a new diagonal matrix using the given vector, filling all
    /// other cells with the filler element.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::diag(vec![1, 2, 3]);
    /// assert!(matrix.layout().is_square());
    /// assert!(matrix.layout().rows() == 3);
    /// assert!(matrix[(0, 0)] == 1);
    /// assert!(matrix[(1, 1)] == 2);
    /// ```
    ///
    pub fn diag(vec: Vec<T>) -> Self {
        let layout = MatrixLayout::square(vec.len());

        let mut matrix = Matrix::zeroed(layout);
        for i in 0..vec.len() {
            matrix[(i, i)] = vec[i];
        }

        matrix
    }

    ///
    /// Creates a new eye-matrix with the given eye element
    /// filling the remaing cells with the filler element.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let eye = Matrix::<usize>::eye(3);
    /// assert!(eye.layout().is_square());
    /// assert!(eye.layout().rows() == 3);
    /// assert!(eye[(1, 1)] == 1);
    /// ```
    ///
    pub fn eye(size: usize) -> Self {
        let layout = MatrixLayout::square(size);

        let mut matrix = Matrix::zeroed(layout);
        for i in 0..size {
            matrix[(i, i)] = T::one();
        }

        matrix
    }

    ///
    /// Resizes the matrix according to the new layout given.
    /// Note that this operation happens in-place should the number of collums
    /// not change, or out-of-place should they change.
    ///
    /// The new cells in the new matrix will be filled with zero elemtents.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let mut matrix = matrix![
    ///     1, 2;
    ///     3, 4;
    /// ];
    /// matrix.resize((3, 3));
    ///
    /// assert_eq!(matrix.layout(), MatrixLayout::new(3, 3));
    /// assert_eq!(matrix[(0, 0)], 1);
    /// assert_eq!(matrix[(1, 1)], 4);
    /// assert_eq!(matrix[(2, 2)], 0usize);
    /// ```
    ///
    pub fn resize<U>(&mut self, new_layout: U)
    where
        U: Into<MatrixLayout>,
    {
        let new_layout = new_layout.into();
        let rows = new_layout.rows();
        let cols = new_layout.cols();

        if self.layout.cols == cols {
            if self.layout.rows > rows {
                self.raw.truncate(rows * cols);
            } else {
                self.raw
                    .extend(vec![T::zero(); cols * (rows - self.layout.rows)]);
            }
            self.layout.rows = rows;
        } else {
            let mut matrix = Self::zeroed(new_layout);
            let rows = self.layout.rows.min(rows);
            let columns = self.layout.cols.min(cols);
            for j in 0..columns {
                for i in 0..rows {
                    matrix[(i, j)] = self[(i, j)];
                }
            }
            *self = matrix;
        }
    }
}

///
/// Conversion from structured two-dimensional vector.
///
/// This conversion will try to assume the matrix layout based on the given
/// two-dimensional vector. If the vector has an irregular shape, meaning rows
/// with different number of elements the conversion will fail.
/// If not the vector will be flattend and saved as a concrete slice.
///
/// The Errors are represented as static strings (&'static str).
///

impl<T: Num> TryFrom<Vec<Vec<T>>> for Matrix<T> {
    type Error = &'static str;

    fn try_from(structured: Vec<Vec<T>>) -> Result<Self, &'static str> {
        let layout = MatrixLayout::assumed_from(&structured)?;
        Ok(Matrix {
            layout,
            raw: structured.into_iter().flatten().collect(),
        })
    }
}

///
/// Conversion from raw vector.
///
/// This conversion will create a defaulted collum vector matrix, using the
/// given vector as elements.
///

impl<T> From<Vec<T>> for Matrix<T> {
    fn from(raw: Vec<T>) -> Self {
        Matrix {
            layout: MatrixLayout {
                rows: raw.len(),
                cols: 1,
            },
            raw,
        }
    }
}

impl<T> From<Matrix<T>> for Vec<T> {
    fn from(matrix: Matrix<T>) -> Self {
        matrix.raw
    }
}

///
/// Matrix: Custom Trait Implementations
///

impl<T: Num> Default for Matrix<T>
where
    T: Default,
{
    fn default() -> Self {
        Self {
            layout: MatrixLayout::default(),
            raw: Vec::new(),
        }
    }
}

impl<T: Num> IntoIterator for Matrix<T> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.raw.into_iter()
    }
}

impl<T> Deref for Matrix<T> {
    type Target = [T];

    fn deref(&self) -> &Self::Target {
        self.raw.deref()
    }
}

impl<T> Index<usize> for Matrix<T> {
    type Output = T;

    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.raw[index]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.raw[index]
    }
}

impl<T> DerefMut for Matrix<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.raw.deref_mut()
    }
}

impl<T> Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    #[inline(always)]
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.raw[self.layout.index(index)]
    }
}

impl<T> IndexMut<(usize, usize)> for Matrix<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.raw[self.layout.index(index)]
    }
}

impl<T> Index<[usize; 2]> for Matrix<T> {
    type Output = T;

    #[inline(always)]
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        &self.raw[self.layout.index((index[0], index[1]))]
    }
}

impl<T> IndexMut<[usize; 2]> for Matrix<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        &mut self.raw[self.layout.index((index[0], index[1]))]
    }
}

impl<T> PartialEq for Matrix<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.raw == other.raw && self.layout == other.layout
    }
}

impl<T: Num> Eq for Matrix<T> where T: Eq {}

impl<T: Num> Hash for Matrix<T>
where
    T: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.layout.hash(state);
        self.raw.hash(state);
    }
}

impl<T: Num> PartialOrd for Matrix<T>
where
    T: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        assert!(self.layout == other.layout);
        assert!(self.layout.is_vec());

        let mut err = false;
        let mut eq = true;
        let mut lhs_g = false;

        // SAFTY: Loop is safe since (layout.size == raw.len)
        for i in 0..self.layout.size() {
            match self.raw[i].partial_cmp(&other.raw[i]) {
                Some(std::cmp::Ordering::Greater) => {
                    lhs_g = true;
                    eq = false;
                }
                Some(std::cmp::Ordering::Less) => {
                    eq = false;
                }
                Some(_) => {}
                None => {
                    err = true;
                }
            }
        }

        match (err, eq) {
            (true, _) => None,
            (false, true) => Some(std::cmp::Ordering::Equal),
            (false, false) => {
                if lhs_g {
                    Some(std::cmp::Ordering::Greater)
                } else {
                    Some(std::cmp::Ordering::Less)
                }
            }
        }
    }

    fn lt(&self, other: &Self) -> bool {
        assert_eq!(self.layout(), other.layout());
        assert!(self.layout().is_vec());

        for k in 0..self.layout().size() {
            if self.raw[k] >= other.raw[k] {
                return false;
            }
        }

        true
    }

    fn le(&self, other: &Self) -> bool {
        assert_eq!(self.layout(), other.layout());
        assert!(self.layout().is_vec());

        for k in 0..self.layout().size() {
            if self.raw[k] > other.raw[k] {
                return false;
            }
        }

        true
    }

    fn gt(&self, other: &Self) -> bool {
        assert_eq!(self.layout(), other.layout());
        assert!(self.layout().is_vec());

        for k in 0..self.layout().size() {
            if self.raw[k] <= other.raw[k] {
                return false;
            }
        }

        true
    }

    fn ge(&self, other: &Self) -> bool {
        assert_eq!(self.layout(), other.layout());
        assert!(self.layout().is_vec());

        for k in 0..self.layout().size() {
            if self.raw[k] < other.raw[k] {
                return false;
            }
        }

        true
    }
}

///
/// Matrix: Addition
///

impl<T> Add<Matrix<T>> for Matrix<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn add(self, rhs: Matrix<T>) -> Self::Output {
        Add::add(self, &rhs)
    }
}

impl<T> Add<&'_ Matrix<T>> for Matrix<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn add(mut self, rhs: &'_ Matrix<T>) -> Self::Output {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i] = self.raw[i].add(rhs.raw[i]);
        }
        self
    }
}

impl<T> Add<Matrix<T>> for &'_ Matrix<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn add(self, rhs: Matrix<T>) -> Self::Output {
        Add::add(self.clone(), &rhs)
    }
}

impl<T> Add<&'_ Matrix<T>> for &'_ Matrix<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn add(self, rhs: &'_ Matrix<T>) -> Self::Output {
        Add::add(self.clone(), rhs)
    }
}

impl<T> AddAssign<Matrix<T>> for Matrix<T>
where
    T: AddAssign + Copy,
{
    fn add_assign(&mut self, rhs: Matrix<T>) {
        AddAssign::add_assign(self, &rhs)
    }
}

impl<T> AddAssign<&'_ Matrix<T>> for Matrix<T>
where
    T: AddAssign + Copy,
{
    fn add_assign(&mut self, rhs: &'_ Matrix<T>) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].add_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Subtraction
///

impl<T> Sub<Matrix<T>> for Matrix<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn sub(self, rhs: Matrix<T>) -> Self::Output {
        Sub::sub(self, &rhs)
    }
}

impl<T> Sub<&'_ Matrix<T>> for Matrix<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn sub(mut self, rhs: &'_ Matrix<T>) -> Self::Output {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i] = self.raw[i].sub(rhs.raw[i]);
        }
        self
    }
}

impl<T> Sub<Matrix<T>> for &'_ Matrix<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn sub(self, rhs: Matrix<T>) -> Self::Output {
        Sub::sub(self.clone(), &rhs)
    }
}

impl<T> Sub<&'_ Matrix<T>> for &'_ Matrix<T>
where
    T: Sub<Output = T> + Copy,
{
    type Output = Matrix<T>;

    fn sub(self, rhs: &'_ Matrix<T>) -> Self::Output {
        Sub::sub(self.clone(), rhs)
    }
}

impl<T> SubAssign<Matrix<T>> for Matrix<T>
where
    T: SubAssign + Copy,
{
    fn sub_assign(&mut self, rhs: Matrix<T>) {
        SubAssign::sub_assign(self, &rhs)
    }
}

impl<T> SubAssign<&'_ Matrix<T>> for Matrix<T>
where
    T: SubAssign + Copy,
{
    fn sub_assign(&mut self, rhs: &'_ Matrix<T>) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].sub_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Bitwise AND
///

impl<T: Num> BitAnd for Matrix<T>
where
    T: BitAnd + Copy,
{
    type Output = Matrix<<T as BitAnd>::Output>;

    fn bitand(self, rhs: Self) -> Self::Output {
        assert!(self.layout == rhs.layout);

        // SAFTY: Since the iteration will follow the raw vector thus will fill
        // all elements. Since layouts are equal, all operands are indeed provided
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for i in 0..result.raw.len() {
            result.raw[i] = self.raw[i].bitand(rhs.raw[i]);
        }

        result
    }
}

impl<T: Num> BitAndAssign for Matrix<T>
where
    T: BitAndAssign + Copy,
{
    fn bitand_assign(&mut self, rhs: Self) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].bitand_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Bitwise OR
///

impl<T: Num> BitOr for Matrix<T>
where
    T: BitOr + Copy,
{
    type Output = Matrix<<T as BitOr>::Output>;

    fn bitor(self, rhs: Self) -> Self::Output {
        assert!(self.layout == rhs.layout);

        // SAFTY: Since the iteration will follow the raw vector thus will fill
        // all elements. Since layouts are equal, all operands are indeed provided
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for i in 0..result.raw.len() {
            result.raw[i] = self.raw[i].bitor(rhs.raw[i]);
        }

        result
    }
}

impl<T: Num> BitOrAssign for Matrix<T>
where
    T: BitOrAssign + Copy,
{
    fn bitor_assign(&mut self, rhs: Self) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].bitor_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Bitwise XOR
///

impl<T: Num> BitXor for Matrix<T>
where
    T: BitXor + Copy,
{
    type Output = Matrix<<T as BitXor>::Output>;

    fn bitxor(self, rhs: Self) -> Self::Output {
        assert!(self.layout == rhs.layout);

        // SAFTY: Since the iteration will follow the raw vector thus will fill
        // all elements. Since layouts are equal, all operands are indeed provided
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for i in 0..result.raw.len() {
            result.raw[i] = self.raw[i].bitxor(rhs.raw[i]);
        }

        result
    }
}

impl<T: Num> BitXorAssign for Matrix<T>
where
    T: BitXorAssign + Copy,
{
    fn bitxor_assign(&mut self, rhs: Self) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].bitxor_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Misc Traits
///

impl<T: Num> Display for Matrix<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Matrix<{}, {}> {{", self.layout.rows, self.layout.cols)?;

        for i in 0..self.layout.rows {
            for j in 0..self.layout.cols {
                write!(f, " {}", self[(i, j)])?;
                // str.push_str(&format!(" {}", self[(i, j)]))
            }
            writeln!(f)?;
        }

        write!(f, "}}")
    }
}

///
/// # Generic algebraic operations
///
/// A collection of algebraic operation that are generic and only bound by
/// internal Rust-Traits.
/// This grouping does not include more specific operations (e.g. det() or eig())
/// which will be implemented in a specific way to improve computational efficenicy.
///

impl<T: Num> Matrix<T>
where
    T: Add<Output = T> + Copy + Default,
{
    ///
    /// Calculates the sum of diagonal elements (the trace) of a square matrix
    /// using the default element as zero.
    ///
    /// # Panics
    ///
    /// This function will panic if applied to a non-square matrix.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let diag = Matrix::diag(vec![1, 2, 3]);
    /// assert!(diag.layout() == MatrixLayout::new(3, 3));
    /// assert!(diag.trace() == 6);
    /// ```
    ///
    pub fn trace(&self) -> T {
        assert!(self.layout.is_square());
        let mut sum = T::default();
        for i in 0..self.layout.rows {
            sum = sum.add(self[(i, i)])
        }
        sum
    }
}

impl<T: Num> Matrix<T>
where
    T: Copy,
{
    ///
    /// Transposes the matrix with a out-of-place new buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// use linalg::prelude::*;
    /// use std::convert::TryFrom;
    ///
    /// let mut matrix = Matrix::<usize>::try_from(vec![
    ///     vec![1, 2, 3usize],
    ///     vec![4, 5, 6],
    /// ]).unwrap();
    /// matrix.transpose();
    ///
    /// assert!(matrix.layout() == MatrixLayout::new(3, 2));
    /// assert!(matrix[(2, 1)] == 6);
    /// ```
    ///
    pub fn transpose(&mut self) {
        if self.layout().is_colvec() || self.layout().is_rowvec() {
            // For vectorized matrices transposition will not change the raw buffer
            self.layout.transpose()
        } else {
            *self = self.transposed()
        }
    }

    ///
    /// Returns a new out-of-place transposed matrix.
    ///
    /// # Examples
    ///
    /// ```
    /// use linalg::prelude::*;
    /// use std::convert::TryFrom;
    ///
    /// let matrix = Matrix::<usize>::try_from(vec![
    ///     vec![1, 2, 3usize],
    ///     vec![4, 5, 6],
    /// ]).unwrap();
    /// let matrix = matrix.transposed();
    ///
    /// assert!(matrix.layout() == MatrixLayout::new(3, 2));
    /// assert!(matrix[(2, 1)] == 6);
    /// ```
    ///
    pub fn transposed(&self) -> Self {
        // SAFTY:
        // Since the tranposed matrix has the transposed layout
        // and thus the same size requirements for the raw vector
        // there must be a eqivalent value for each cell in the original matrix
        let mut transposed = unsafe { Matrix::uninitalized(self.layout.transposed()) };
        self.transpose_into(&mut transposed);

        transposed
    }

    ///
    /// Writes the tranposed matrix into the given target,
    /// assumming the buffer is already pre-allocted, either as allready used
    /// memory or as uninitialied memory.
    ///
    /// # Panics
    ///
    /// This function assumms that the target has a layout that is
    /// at least sized equally (layout.size()) to provide te correct amount of
    /// memory.
    ///
    pub fn transpose_into(&self, target: &mut Self) {
        assert!(target.layout().size() == self.layout().transposed().size());
        target.layout = self.layout().transposed();

        for i in 0..target.layout.rows {
            for j in 0..target.layout.cols {
                target[(i, j)] = self[(j, i)]
            }
        }
    }
}

impl<T: Num> Matrix<T>
where
    T: Mul + Copy,
{
    ///
    /// Performs a scalar multiplication with the given scalar
    /// returning a new matrix as result.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::diag(vec![1, 2, 3]);
    /// let double = matrix.scalar(2);
    /// assert!(matrix.layout() == double.layout());
    /// assert!(2 * matrix[(0, 0)] == double[(0, 0)]);
    /// assert!(2 * matrix[(2, 0)] == double[(2, 0)]);
    /// ```
    ///
    pub fn scalar(&self, scalar: T) -> Matrix<<T as Mul>::Output> {
        // SAFTY:
        // Matix is identical in layout, so all cells will be filled,
        // cause iteration over raw values
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for k in 0..self.layout.size() {
            result.raw[k] = scalar * self.raw[k];
        }

        result
    }
}

impl<T: Num> Mul<T> for Matrix<T>
where
    T: Mul + Copy,
{
    type Output = Matrix<<T as Mul>::Output>;
    fn mul(mut self, rhs: T) -> Self::Output {
        self.scale(rhs);
        self
    }
}

impl<T: Num> Matrix<T>
where
    T: Mul<Output = T> + Copy,
{
    ///
    /// Performs a scalar multiplication in-place with all elements of
    /// the given matrix.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::diag(vec![1, 2, 3]);
    /// let mut double = matrix.clone();
    /// double.scale(2);
    /// assert!(matrix.layout() == double.layout());
    /// assert!(2 * matrix[(0, 0)] == double[(0, 0)]);
    /// assert!(2 * matrix[(2, 0)] == double[(2, 0)]);
    /// ```
    ///
    pub fn scale(&mut self, scalar: T) {
        for k in 0..self.layout.size() {
            self.raw[k] = scalar * self.raw[k];
        }
    }
}

impl<T: Num> MulAssign<T> for Matrix<T>
where
    T: Mul<Output = T> + Copy,
{
    fn mul_assign(&mut self, rhs: T) {
        self.scale(rhs)
    }
}

impl<T: Num> Matrix<T>
where
    T: Copy + Div,
{
    ///
    /// Performs a scalar multiplication with the given scalar
    /// returning a new matrix as result.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::diag(vec![2, 4, 6]);
    /// let half = matrix.scalar_div(2);
    /// assert!(matrix.layout() == half.layout());
    /// assert!(matrix[(0, 0)] / 2 == half[(0, 0)]);
    /// assert!(matrix[(2, 0)] / 2 == half[(2, 0)]);
    /// ```
    ///
    pub fn scalar_div(&self, divider: T) -> Matrix<<T as Div>::Output> {
        // SAFTY:
        // Matix is identical in layout, so all cells will be filled,
        // cause iteration over raw values
        let mut result = unsafe { Matrix::uninitalized(self.layout) };
        for k in 0..self.layout.size() {
            result.raw[k] = self.raw[k] / divider;
        }

        result
    }
}

impl<T: Num> Div<T> for Matrix<T>
where
    T: Div + Copy,
{
    type Output = Matrix<<T as Div>::Output>;
    fn div(mut self, rhs: T) -> Self::Output {
        self.scale_div(rhs);
        self
    }
}

impl<T: Num> Matrix<T>
where
    T: Mul<Output = T> + Copy,
{
    ///
    /// Performs a scalar multiplication in-place with all elements of
    /// the given matrix.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let matrix = Matrix::diag(vec![2, 4, 6]);
    /// let mut half = matrix.clone();
    /// half.scale_div(2);
    /// assert_eq!(matrix.layout(), half.layout());
    /// assert_eq!(matrix[(0, 0)] / 2, half[(0, 0)]);
    /// assert_eq!(matrix[(2, 0)] / 2, half[(2, 0)]);
    /// ```
    ///
    pub fn scale_div(&mut self, scalar: T) {
        for k in 0..self.layout.size() {
            self.raw[k] = self.raw[k] / scalar;
        }
    }
}

impl<T: Num> DivAssign<T> for Matrix<T>
where
    T: Div<Output = T> + Copy,
{
    fn div_assign(&mut self, rhs: T) {
        self.scale_div(rhs)
    }
}

///
/// Matrix: Multiplication
///

impl<T> Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    ///
    /// Performs a matrix mutiplication with the two operands,
    /// and returns the result.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::prelude::*;
    ///
    /// let a = matrix![
    ///     8, -2,  9;
    ///     2,  1, -8;
    ///     4,  -5, 1;
    /// ];
    ///
    /// let b = Matrix::eye(3);
    /// let c = Matrix::mmul(&a, &b);
    ///
    /// assert_eq!(a, c);
    /// ```
    ///
    pub fn mmul(lhs: &Matrix<T>, rhs: &Matrix<T>) -> Matrix<T> {
        let layout = MatrixLayout {
            rows: lhs.layout.rows,
            cols: rhs.layout.cols,
        };

        // SAFTY:
        // Matrix Multiplication gurantees the setting of all
        let mut result = unsafe { Matrix::uninitalized(layout) };
        Matrix::mmul_into(lhs, rhs, &mut result);
        result
    }

    ///
    /// Performs a matrix multiplication using the first, and second arguments
    /// as operands and the thrid operand as target container.
    ///
    /// # Panics
    ///
    /// This function panics if either the operands are incompatible (lhs.cols != rhs.rows),
    /// or the result container does not conform to the expected layout (lhs.rows, rhs.cols).
    ///
    pub fn mmul_into(lhs: &Matrix<T>, rhs: &Matrix<T>, result: &mut Matrix<T>) {
        assert!(lhs.layout.cols == rhs.layout.rows);
        assert!(result.layout.rows == lhs.layout.rows);
        assert!(result.layout.cols == rhs.layout.cols);

        for i in 0..result.layout.rows {
            for j in 0..result.layout.cols {
                let mut sum = T::zero();
                for k in 0..lhs.layout.cols {
                    sum = sum + lhs[(i, k)] * rhs[(k, j)];
                }
                result[(i, j)] = sum;
            }
        }
    }

    //
    // TODO: Optimize mmul under the consideration that lhs is consumed
    // in Mul (sometimes rhs as well)
}

// Matrix x Matrix mutiplication

impl<T> Mul<Matrix<T>> for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    type Output = Matrix<T>;

    fn mul(self, rhs: Matrix<T>) -> Self::Output {
        Matrix::mmul(&self, &rhs)
    }
}

impl<T> Mul<&'_ Matrix<T>> for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    type Output = Matrix<T>;

    fn mul(self, rhs: &'_ Matrix<T>) -> Self::Output {
        Matrix::mmul(&self, rhs)
    }
}

impl<T> Mul<Matrix<T>> for &'_ Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    type Output = Matrix<T>;

    fn mul(self, rhs: Matrix<T>) -> Self::Output {
        Matrix::mmul(self, &rhs)
    }
}

impl<T> Mul<&'_ Matrix<T>> for &'_ Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    type Output = Matrix<T>;

    fn mul(self, rhs: &'_ Matrix<T>) -> Self::Output {
        Matrix::mmul(self, rhs)
    }
}

// Matrix x Vector mutiplication

// TODO:
// To correct impl Matrix * vector the output should be a vector aswell
// -> const parameter cannot be determined since matrix does not support generic size

// impl<T, const N: usize> Mul<Vector<T, N>> for Matrix<T> where
//     T: Mul<Output = T> + Add<Output = T> + Zero + Copy
// {
//     type Output = ;
// }

// Matrix x Matrix Assignment

impl<T: Num> MulAssign<Matrix<T>> for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    fn mul_assign(&mut self, rhs: Matrix<T>) {
        *self = Matrix::mmul(self, &rhs)
    }
}

impl<T: Num> MulAssign<&Matrix<T>> for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Zero + Copy,
{
    fn mul_assign(&mut self, rhs: &Matrix<T>) {
        *self = Matrix::mmul(self, rhs)
    }
}

//
// Matrix Negation
//

impl<T: Num> Neg for Matrix<T>
where
    T: Neg<Output = T> + Copy,
{
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for k in 0..self.layout.size() {
            self.raw[k] = self.raw[k].neg();
        }

        self
    }
}

impl<T: Num> Neg for &Matrix<T>
where
    T: Neg + Copy,
{
    type Output = Matrix<<T as Neg>::Output>;

    fn neg(self) -> Self::Output {
        // SAFTY: All cells will be filled through a raw itr
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for k in 0..self.layout.size() {
            result.raw[k] = self.raw[k].neg();
        }

        result
    }
}

//
// Matrix: Power  TODO
//

impl<T: Copy + Zero> num_traits::Zero for Matrix<T> {
    fn zero() -> Self {
        Matrix::new((1, 1), vec![T::zero()])
    }

    fn is_zero(&self) -> bool {
        self.raw.iter().all(|v| v.is_zero())
    }
}
