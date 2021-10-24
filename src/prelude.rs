//! A reexport of all functions

/// # Marcos
pub use crate::matrix;

/// # Matrix Core
pub use crate::core::Matrix;
pub use crate::core::MatrixLayout;

/// # Discret mathematics functions
pub use crate::discret::bell;
pub use crate::discret::binom;
pub use crate::discret::binom_tbl;
pub use crate::discret::fac;
pub use crate::discret::permut_of_type;
pub use crate::discret::prt_over;
pub use crate::discret::prt_under;
pub use crate::discret::stirling1;
pub use crate::discret::stirling1_tbl;
pub use crate::discret::stirling2;
pub use crate::discret::stirling2_tbl;

/// # Linear system of equationss
pub use crate::lse::det;
pub use crate::lse::eig;
pub use crate::lse::eigv_powitr;
pub use crate::lse::eigv_powitr_inv;
pub use crate::lse::inv;
pub use crate::lse::lr;
pub use crate::lse::lr_dry;
pub use crate::lse::qr;
pub use crate::lse::trace;
pub use crate::lse::LrDecomposition;
pub use crate::lse::QrDecomposition;

/// # Linear optimization
pub use crate::lop::LOPIncrement;
pub use crate::lop::LOPOptions;
pub use crate::lop::LOPSolution;
pub use crate::lop::LinearOpimizationProblem;
pub use crate::lop::LOP;

/// # Polynom utils
pub use crate::poly::pinterpol;
pub use crate::poly::Polynom;
