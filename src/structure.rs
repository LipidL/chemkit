use super::atom::Atom;
use super::cell::CellParameters;
use nalgebra::RealField;
use std::collections::HashMap;
use uom::Conversion;
use uom::si::{SI, Units};

/// A periodic structure in 3D space.
#[derive(Debug, Clone)]
pub struct Structure<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub atoms: Vec<Atom<T>>,
    pub cell: Option<CellParameters<T>>,
    /// Format-specific metadata (e.g. energy, structure index) stored as
    /// key-value pairs of strings.
    pub properties: HashMap<String, String>,
}
