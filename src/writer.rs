mod arc;

use std::error::Error;

use nalgebra::RealField;
use uom::Conversion;
use uom::si::angle::radian;
use uom::si::length::angstrom;
use uom::si::{SI, Units};

use super::structure::Structure;

/// A trait for serialising a single [`Structure<T>`] to a string using a
/// specific file format.
///
/// Mirrors the [`StructParser`](super::parser::StructParser) trait for the
/// writing direction: one implementor per file format, one method that does
/// the work.
pub trait WriteOne {
    type WriteError: Error;

    /// Render a single structure as a string in this format.
    fn write_one<T>(&self, structure: &Structure<T>) -> Result<String, Self::WriteError>
    where
        T: RealField + Copy + Conversion<T, T = T> + std::fmt::Display,
        SI<T>: Units<T>,
        angstrom: Conversion<T, T = T>,
        radian: Conversion<T, T = T>;
}

/// A trait for serialising multiple structures to a single string using a
/// specific file format.
///
/// Multi-frame formats (such as ARC) typically emit a file-level header or
/// other framing that [`WriteOne`] alone cannot produce, so this trait exists
/// alongside `WriteOne` rather than being built on top of it.
pub trait WriteAll {
    type WriteError: Error;

    /// Render a slice of structures as a string in this format.
    fn write_all<T>(&self, structures: &[Structure<T>]) -> Result<String, Self::WriteError>
    where
        T: RealField + Copy + Conversion<T, T = T> + std::fmt::Display,
        SI<T>: Units<T>,
        angstrom: Conversion<T, T = T>,
        radian: Conversion<T, T = T>;
}

pub use arc::ArcWriter;
