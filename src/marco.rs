///
/// A macro to create a row-major matrix using the structs in [matrix](mod@crate::matrix).
///
/// # Syntax
///
/// In this marco elements are sperated by commans, rows are terminated by semicolons.
///
/// # Safty & Panics
///
/// Note that the length of rows should be equal but this is in fact not checked.
/// If you create an invalid matrix the macro will try to find a possible layout
/// matching the number of elements and the given number of rows.
/// If that is possible this layout will be used and a valid, but unintended
/// matrix will be created. If this is not possible the programm will panic.
///
/// # Example
///
/// ```
/// use linalg::prelude::*;
///
/// let macro_gen = matrix![
///     1, 2, 3;
///     4, 5, 6;
/// ];
///
/// let expl = Matrix::new((2, 3), vec![1,2,3,4,5,6]);
///
/// assert_eq!(macro_gen, expl);
/// ```
///
#[macro_export]
macro_rules! matrix {

    () => (
        {
            use linalg::matrix::{MatrixLayout, Matrix};
            Matrix::new(MatrixLayout::new(0, 0), Vec::new())
        }
    );

    ($($item:expr,)* ; -> ($c:expr; [ $($b:expr,)* ])) => (
        {
            let buffer = vec![$($b,)* $($item,)*];

            use linalg::matrix::MatrixLayout;
            let rows = $c + 1;
            let cols = buffer.len() / rows;
            let layout = MatrixLayout::new(rows, cols);

            use linalg::matrix::Matrix;
            Matrix::new(layout, buffer)
        }
    );

    ($($item:expr,)* ; $($($tail:expr,)*;)+ -> ($c:expr; [ $($b:expr,)* ])) => (
        matrix!($($($tail,)*;)+ -> ($c + 1; [ $($b,)* $($item,)* ]))
    );

    ($($($item:expr),*;)*) => (
        matrix!($($($item,)*;)* -> (0; []))
    );
}

// This marco is only used internally
// and thus uses the crate prefix instead of the linalg prefix

#[doc(hidden)]
#[macro_export]
macro_rules! mat {

    () => (
        {
            use $crate::matrix::{MatrixLayout, Matrix};
            Matrix::new(MatrixLayout::new(0, 0), Vec::new())
        }
    );

    ($($item:expr,)* ; -> ($c:expr; [ $($b:expr,)* ])) => (
        {
            let buffer = vec![$($b,)* $($item,)*];

            use $crate::matrix::MatrixLayout;
            let rows = $c + 1;
            let cols = buffer.len() / rows;
            let layout = MatrixLayout::new(rows, cols);

            use $crate::matrix::Matrix;
            Matrix::new(layout, buffer)
        }
    );

    ($($item:expr,)* ; $($($tail:expr,)*;)+ -> ($c:expr; [ $($b:expr,)* ])) => (
        mat!($($($tail,)*;)+ -> ($c + 1; [ $($b,)* $($item,)* ]))
    );

    ($($($item:expr),*;)*) => (
        mat!($($($item,)*;)* -> (0; []))
    );
}
