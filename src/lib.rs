//! Rust implementation of the public `lrcalc` library surface.
//!
//! The near-term target is ABI compatibility with Anders Buch's C library.
//! Native Rust APIs should remain separate from the C-facing allocation and
//! iterator types in [`abi`].

#![deny(unsafe_op_in_unsafe_fn)]

pub mod abi;
pub mod kostka;
pub mod kostka_fast;
pub mod lr_ehrhart;
pub mod lr_gt;
pub mod lr_signed;
pub mod lrcoef;
pub mod partition;
pub mod schubert;
pub mod schur;

pub use abi::{iv_free, iv_hash, iv_new, iv_new_copy, iv_new_zero, ivlc_free_all, IVector};
