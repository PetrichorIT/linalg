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
    fmt::Display,
    hash::Hash,
    mem::swap,
    ops::{
        Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Index,
        IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
    },
};

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
/// use linalg::matrix::*;
///
/// let layout = MatrixLayout::new(2, 3);
/// let matrix = Matrix::fill(layout, 0usize);
/// assert!(matrix.layout().rows() == 2 && matrix.layout().cols() == 3);
/// assert!(matrix[(0, 0)] == 0usize);
/// ```
///
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MatrixLayout {
    rows: usize,
    cols: usize,
}

impl MatrixLayout {
    /// The number of rows in the matrix.
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// The number of collums in the matrix.
    #[inline]
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

    ///
    /// A conversion function from a two-dimensional index to a
    /// one-dimensional index for the raw buffer.
    /// This conversion follows the following formular:
    /// BufferIndex = Index_Row * Collum_Size + Index_Collum
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::matrix::*;
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
    /// use linalg::matrix::*;
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
    pub fn assumed_from<T>(structured: &Vec<Vec<T>>) -> Result<Self, &'static str> {
        let rows = structured.len();
        if rows == 0 {
            return Ok(MatrixLayout { rows: 0, cols: 0 });
        }

        let cols = structured.first().unwrap().len();
        for row in structured.iter().skip(1) {
            if row.len() != cols {
                return Err(&"Structured data has irregulare format");
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

impl Default for MatrixLayout {
    fn default() -> Self {
        Self { rows: 0, cols: 0 }
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
/// a raw buffer to store matrix cells. The size of the raw buffer thereby is
/// defined by the layouts [MatrixLayout::size()] function.
///
/// # Example
///
/// ```
/// use linalg::matrix::*;
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
    /// The layout constraints definig the matrix.
    #[inline(always)]
    pub fn layout(&self) -> &MatrixLayout {
        &self.layout
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
}

impl<T> Matrix<T> {
    ///
    /// Creates a new matrix with a given layout over the given buffer.
    ///  
    /// # Panics
    ///
    /// This function panics should the layout not fit the buffer.
    /// This can be archived should the layouts size match the buffers length.
    ///
    /// ```should_panic
    /// use linalg::matrix::*;
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

    ///
    /// Creates a new matrix with an allocated but uninialized
    /// memory buffer.
    ///
    /// # Safty
    ///
    /// Ensure that all cells of the allocated buffer (self.raw().len())
    /// are either filled with valid instances of the given type,
    /// or never used.
    ///
    pub unsafe fn uninitalized(layout: MatrixLayout) -> Self {
        let mut raw = Vec::with_capacity(layout.size());

        // SAFTY:
        // This function will only be called internally in its 'unsafe' but not
        // outwardly unsafe shape, to create buffers and such.
        raw.set_len(layout.size());

        Self { layout, raw }
    }
}

impl<T> Matrix<T>
where
    T: Copy,
{
    ///
    /// Creates a new matrix with the given layout filling all cells
    /// with the given filler element.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::matrix::*;
    ///
    /// let matrix = Matrix::fill(MatrixLayout::new(2, 3), 0usize);
    /// assert!(matrix[(0, 0)] == 0usize);
    /// assert!(matrix[(1, 2)] == 0usize);
    /// ```
    ///
    pub fn fill(layout: MatrixLayout, filler: T) -> Self {
        // SAFTY:
        // Can be used since the underling vector will be filled according to its
        // own len, thus all cells used are guarnteed to be there
        // since the vector has the given size AND capacity
        let mut matrix = unsafe { Self::uninitalized(layout) };
        for cell in &mut matrix.raw {
            *cell = filler;
        }
        matrix
    }

    ///
    /// Create a new diagonal matrix using the given vector, filling all
    /// other cells with the filler element.
    ///
    /// # Example
    ///
    /// ```
    /// use linalg::matrix::*;
    ///
    /// let matrix = Matrix::diag(vec![1, 2, 3], 0usize);
    /// assert!(matrix.layout().is_square());
    /// assert!(matrix.layout().rows() == 3);
    /// assert!(matrix[(0, 0)] == 1);
    /// assert!(matrix[(1, 1)] == 2);
    /// ```
    ///
    pub fn diag(vec: Vec<T>, filler: T) -> Self {
        let layout = MatrixLayout::square(vec.len());

        let mut matrix = Matrix::fill(layout, filler);
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
    /// use linalg::matrix::*;
    ///
    /// let eye = Matrix::eye(3, 1, 0usize);
    /// assert!(eye.layout().is_square());
    /// assert!(eye.layout().rows() == 3);
    /// assert!(eye[(1, 1)] == 1);
    /// ```
    ///
    pub fn eye(size: usize, eye: T, filler: T) -> Self {
        let layout = MatrixLayout::square(size);

        let mut matrix = Matrix::fill(layout, filler);
        for i in 0..size {
            matrix[(i, i)] = eye;
        }

        matrix
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

impl<T> TryFrom<Vec<Vec<T>>> for Matrix<T> {
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

///
/// Matrix: Custom Trait Implementations
///

impl<T> Default for Matrix<T>
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

impl<T> IntoIterator for Matrix<T> {
    type Item = T;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.raw.into_iter()
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

impl<T> PartialEq for Matrix<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.raw == other.raw && self.layout == other.layout
    }
}

impl<T> Eq for Matrix<T> where T: Eq {}

impl<T> Hash for Matrix<T>
where
    T: Hash,
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.layout.hash(state);
        self.raw.hash(state);
    }
}

impl<T> PartialOrd for Matrix<T>
where
    T: PartialOrd,
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        assert!(self.layout == other.layout);
        assert!(self.layout.is_colvec() || self.layout.is_rowvec());

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
}

///
/// Matrix: Addition
///

impl<T> Add for Matrix<T>
where
    T: Add + Copy,
{
    type Output = Matrix<T::Output>;

    fn add(self, rhs: Self) -> Self::Output {
        assert!(self.layout == rhs.layout);

        // SAFTY: Since the iteration will follow the raw vector thus will fill
        // all elements. Since layouts are equal, all operands are indeed provided
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for i in 0..result.raw.len() {
            result.raw[i] = self.raw[i].add(rhs.raw[i]);
        }

        result
    }
}

impl<T> AddAssign for Matrix<T>
where
    T: AddAssign + Copy,
{
    fn add_assign(&mut self, rhs: Self) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].add_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Subtraction
///

impl<T> Sub for Matrix<T>
where
    T: Sub + Copy,
{
    type Output = Matrix<T::Output>;

    fn sub(self, rhs: Self) -> Self::Output {
        assert!(self.layout == rhs.layout);

        // SAFTY: Since the iteration will follow the raw vector thus will fill
        // all elements. Since layouts are equal, all operands are indeed provided
        let mut result = unsafe { Matrix::uninitalized(self.layout) };

        for i in 0..result.raw.len() {
            result.raw[i] = self.raw[i].sub(rhs.raw[i]);
        }

        result
    }
}

impl<T> SubAssign for Matrix<T>
where
    T: SubAssign + Copy,
{
    fn sub_assign(&mut self, rhs: Self) {
        assert!(self.layout == rhs.layout);
        for i in 0..self.raw.len() {
            self.raw[i].sub_assign(rhs.raw[i])
        }
    }
}

///
/// Matrix: Bitwise AND
///

impl<T> BitAnd for Matrix<T>
where
    T: BitAnd + Copy,
{
    type Output = Matrix<T::Output>;

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

impl<T> BitAndAssign for Matrix<T>
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

impl<T> BitOr for Matrix<T>
where
    T: BitOr + Copy,
{
    type Output = Matrix<T::Output>;

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

impl<T> BitOrAssign for Matrix<T>
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

impl<T> BitXor for Matrix<T>
where
    T: BitXor + Copy,
{
    type Output = Matrix<T::Output>;

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

impl<T> BitXorAssign for Matrix<T>
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
/// Matrix: Multiplication
///

impl<T> Mul for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Copy + Default,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        assert!(self.layout.cols == rhs.layout.rows);

        let layout = MatrixLayout {
            rows: self.layout.rows,
            cols: rhs.layout.cols,
        };

        // SAFTY:
        // Matrix Multiplication gurantees the setting of all
        let mut result = unsafe { Matrix::uninitalized(layout) };

        for i in 0..result.layout.rows {
            for j in 0..result.layout.cols {
                let mut sum = T::default();
                for k in 0..self.layout.cols {
                    sum = sum + self[(i, k)] * rhs[(k, j)];
                }
                result[(i, j)] = sum;
            }
        }

        result
    }
}

impl<T> MulAssign for Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Copy + Default,
{
    fn mul_assign(&mut self, rhs: Self) {
        assert!(self.layout.cols == rhs.layout.rows);

        let layout = MatrixLayout {
            rows: self.layout.rows,
            cols: rhs.layout.cols,
        };

        // SAFTY:
        // Matrix Multiplication gurantees the setting of all
        let mut result = unsafe { Matrix::uninitalized(layout) };

        for i in 0..result.layout.rows {
            for j in 0..result.layout.cols {
                let mut sum = T::default();
                for k in 0..self.layout.cols {
                    sum = sum + self[(i, k)] * rhs[(k, j)];
                }
                result[(i, j)] = sum;
            }
        }

        *self = result
    }
}

///
/// Matrix: Misc Traits
///

impl<T> Display for Matrix<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut str = String::with_capacity(self.layout.size() * 2);
        str.push_str(&format!(
            "Matrix<{}, {}> {{\n",
            self.layout.rows, self.layout.cols
        ));

        for i in 0..self.layout.rows {
            for j in 0..self.layout.cols {
                str.push_str(&format!(" {}", self[(i, j)]))
            }
            str.push('\n');
        }

        str.push_str("}}");

        write!(f, "{}", str)
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

impl<T> Matrix<T>
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
    /// use linalg::matrix::*;
    ///
    /// let diag = Matrix::diag(vec![1, 2, 3], 42usize);
    /// assert!(*diag.layout() == MatrixLayout::new(3, 3));
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

impl<T> Matrix<T>
where
    T: Copy,
{
    ///
    /// Transposes the matrix with a out-of-place new buffer.
    ///
    /// # Examples
    ///
    /// ```
    /// use linalg::matrix::*;
    /// use std::convert::TryFrom;
    ///
    /// let mut matrix = Matrix::<usize>::try_from(vec![
    ///     vec![1, 2, 3usize],
    ///     vec![4, 5, 6],
    /// ]).unwrap();
    /// matrix.transpose();
    ///
    /// assert!(*matrix.layout() == MatrixLayout::new(3, 2));
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
    /// use linalg::matrix::*;
    /// use std::convert::TryFrom;
    ///
    /// let matrix = Matrix::<usize>::try_from(vec![
    ///     vec![1, 2, 3usize],
    ///     vec![4, 5, 6],
    /// ]).unwrap();
    /// let matrix = matrix.transposed();
    ///
    /// assert!(*matrix.layout() == MatrixLayout::new(3, 2));
    /// assert!(matrix[(2, 1)] == 6);
    /// ```
    ///
    pub fn transposed(&self) -> Self {
        // SAFTY:
        // Since the tranposed matrix has the transposed layout
        // and thus the same size requirements for the raw vector
        // there must be a eqivalent value for each cell in the original matrix
        let mut transposed = unsafe { Matrix::uninitalized(self.layout.transposed()) };

        for i in 0..transposed.layout.rows {
            for j in 0..transposed.layout.cols {
                transposed[(i, j)] = self[(j, i)]
            }
        }

        transposed
    }
}

impl<T> Matrix<T>
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
    /// use linalg::matrix::*;
    ///
    /// let matrix = Matrix::diag(vec![1, 2, 3], 0usize);
    /// let double = matrix.scalar(2);
    /// assert!(matrix.layout() == double.layout());
    /// assert!(2 * matrix[(0, 0)] == double[(0, 0)]);
    /// assert!(2 * matrix[(2, 0)] == double[(2, 0)]);
    /// ```
    ///
    pub fn scalar(&self, scalar: T) -> Matrix<T::Output> {
        // SAFTY:
        // Matix is identical in layout, so all cells will be filled,
        // cause iteration over raw values
        let mut result = unsafe { Matrix::uninitalized(self.layout.clone()) };

        for k in 0..self.layout.size() {
            result.raw[k] = scalar * self.raw[k];
        }

        result
    }
}

impl<T> Matrix<T>
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
    /// use linalg::matrix::*;
    ///
    /// let matrix = Matrix::diag(vec![1, 2, 3], 0usize);
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

impl<T> Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + Copy + Default,
{
    pub fn mmul(&self, rhs: Self) -> Self {
        assert!(self.layout.cols == rhs.layout.rows);

        let layout = MatrixLayout {
            rows: self.layout.rows,
            cols: rhs.layout.cols,
        };

        // SAFTY:
        // Matrix Multiplication gurantees the setting of all
        let mut result = unsafe { Matrix::uninitalized(layout) };

        for i in 0..result.layout.rows {
            for j in 0..result.layout.cols {
                let mut sum = T::default();
                for k in 0..self.layout.cols {
                    sum = sum + self[(i, k)] * rhs[(k, j)];
                }
                result[(i, j)] = sum;
            }
        }

        result
    }
}

impl<T> Neg for Matrix<T>
where
    T: Neg + Copy,
{
    type Output = Matrix<T::Output>;

    fn neg(self) -> Self::Output {
        // SAFTY:
        // Matrix shares same layout thus all cells will be initalized with a raw iteration.
        let mut result = unsafe { Matrix::uninitalized(self.layout.clone()) };

        for k in 0..self.layout.size() {
            result.raw[k] = self.raw[k].neg();
        }

        result
    }
}

//
// Matrix: Power  TODO
//
