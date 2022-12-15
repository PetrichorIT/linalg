//! Tools for the creation, evaluation and usage of codes.

pub mod huffman;
pub mod permut;

mod source;
pub use source::*;

mod code;
pub use code::*;

mod linear;
pub use linear::*;

#[cfg(test)]
mod tests;
