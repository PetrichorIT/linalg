use std::{
    convert::TryFrom,
    fmt::Display,
    hash::Hash,
    mem::swap,
    ops::{
        Add, AddAssign, BitAnd, BitAndAssign, BitOr, BitOrAssign, BitXor, BitXorAssign, Index,
        IndexMut, Mul,
    },
};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct MatrixLayout {
    rows: usize,
    cols: usize,
}

impl MatrixLayout {
    #[inline(always)]
    pub fn size(&self) -> usize {
        self.rows * self.cols
    }

    #[inline(always)]
    pub fn is_square(&self) -> bool {
        self.rows == self.cols
    }

    #[inline(always)]
    pub fn is_rowvec(&self) -> bool {
        self.cols == 1
    }

    #[inline(always)]
    pub fn is_colvec(&self) -> bool {
        self.rows == 1
    }

    #[inline(always)]
    pub fn index(&self, index: (usize, usize)) -> usize {
        index.0 * self.cols + index.1
    }

    pub fn transpose(&mut self) {
        swap(&mut self.rows, &mut self.cols)
    }

    pub fn transposed(&self) -> Self {
        MatrixLayout {
            rows: self.cols,
            cols: self.rows,
        }
    }
}

impl MatrixLayout {
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

///
/// Matrix
///

#[derive(Debug, Clone)]
pub struct Matrix<T> {
    pub layout: MatrixLayout,
    pub raw: Vec<T>,
}

///
/// Matrix: Construction
///

impl<T> Matrix<T> {
    fn uninitalized(layout: MatrixLayout) -> Self {
        let mut raw = Vec::with_capacity(layout.size());

        // SAFTY:
        // This function will only be called internally in its 'unsafe' but not
        // outwardly unsafe shape, to create buffers and such.
        unsafe { raw.set_len(layout.size()) }

        Self { layout, raw }
    }
}

impl<T> Matrix<T>
where
    T: Copy,
{
    pub fn fill(layout: MatrixLayout, filler: T) -> Self {
        // SAFTY:
        // Can be used since the underling vector will be filled according to its
        // own len, thus all cells used are guarnteed to be there
        // since the vector has the given size AND capacity
        let mut matrix = Self::uninitalized(layout);
        for cell in &mut matrix.raw {
            *cell = filler;
        }
        matrix
    }

    pub fn diag(vec: Vec<T>, filler: T) -> Self {
        let layout = MatrixLayout::square(vec.len());

        let mut matrix = Matrix::fill(layout, filler);
        for i in 0..vec.len() {
            matrix[(i, i)] = vec[i];
        }

        matrix
    }

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
/// Matrix: Operations
///

impl<T> Matrix<T>
where
    T: Copy,
{
    pub fn transpose(&mut self) {
        *self = self.transposed()
    }

    pub fn transposed(&self) -> Self {
        // SAFTY:
        // Since the tranposed matrix has the transposed layout
        // and thus the same size requirements for the raw vector
        // there must be a eqivalent value for each cell in the original matrix
        let mut transposed = Matrix::uninitalized(self.layout.transposed());

        for i in 0..transposed.layout.rows {
            for j in 0..transposed.layout.cols {
                transposed[(i, j)] = self[(j, i)]
            }
        }

        transposed
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
        let mut result = Matrix::uninitalized(self.layout);

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
        let mut result = Matrix::uninitalized(self.layout);

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
        let mut result = Matrix::uninitalized(self.layout);

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
        let mut result = Matrix::uninitalized(self.layout);

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
        let mut result = Matrix::uninitalized(layout);

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
