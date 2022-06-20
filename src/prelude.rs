//! A reexport of all functions

/// # Marcos
pub use crate::matrix;

/// # Matrix Core
pub use crate::matrix::Matrix;
pub use crate::matrix::MatrixLayout;

/// # Discret mathematics functions
pub use crate::discret::bell;
pub use crate::discret::binom;
pub use crate::discret::fac;
pub use crate::discret::permut_of_type;
pub use crate::discret::prt_over;
pub use crate::discret::prt_under;
pub use crate::discret::stirling1;
pub use crate::discret::stirling2;
pub use crate::discret::Stirling1Tbl;
pub use crate::discret::Stirling2Tbl;

/// # Linear system of equationss
pub use crate::algebra::det;
pub use crate::algebra::eig;
pub use crate::algebra::eigv_powitr;
pub use crate::algebra::eigv_powitr_inv;
pub use crate::algebra::inv;
pub use crate::algebra::lr;
pub use crate::algebra::lr_dry;
pub use crate::algebra::qr;
pub use crate::algebra::trace;
pub use crate::algebra::LrDecomposition;
pub use crate::algebra::QrDecomposition;

/// # Linear optimization
pub use crate::lop::LOPIncrement;
pub use crate::lop::LOPOptions;
pub use crate::lop::LOPSolution;
pub use crate::lop::LinearOpimizationProblem;
pub use crate::lop::LOP;

/// # Polynom utils
pub use crate::poly::pdiv;
pub use crate::poly::pinterpol;
pub use crate::poly::Polynom;

/// # Sets
pub use crate::sets::GF;
