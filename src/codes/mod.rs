//! Tools for the creation, evaluation and usage of codes.

pub mod huffman;
pub mod permut;

mod source;
pub use source::*;

mod code;
pub use code::*;

#[cfg(test)]
mod tests;
