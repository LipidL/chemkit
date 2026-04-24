mod arc;

use std::collections::HashMap;
use std::marker::PhantomData;

use nalgebra::RealField;
use thiserror::Error;
use uom::si::angle::radian;
use uom::{
    Conversion,
    si::{SI, Units},
};

use super::atom::Atom;
use super::cell::CellParameters;
use super::periodic_table::Element;
use super::structure::Structure;
use crate::coordinate::{Cartesian, CoordinateSystem};

/// How atom positions are expressed in the raw file
pub(crate) enum PositionKind {
    /// Atom positions are given in Cartesian coordinates
    Cartesian,
    /// Atom positions are given in fractional coordinates
    Fractional,
}

/// Raw lattice – each column of `matrix` is a lattice vector
pub(crate) struct RawLattice {
    pub matrix: [[f64; 3]; 3],
}

impl RawLattice {
    /// Construct a `RawLattice` from conventional cell parameters.
    ///
    /// `a`, `b`, `c` are in the caller's length unit.
    /// `alpha`, `beta`, `gamma` are in **degrees**.
    pub fn from_cell_parameters(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        let alpha = alpha.to_radians();
        let beta = beta.to_radians();
        let gamma = gamma.to_radians();

        let cos_alpha = alpha.cos();
        let cos_beta = beta.cos();
        let cos_gamma = gamma.cos();
        let sin_gamma = gamma.sin();

        let cx = cos_beta;
        let cy = (cos_alpha - cos_beta * cos_gamma) / sin_gamma;
        let cz = (1.0 - cx * cx - cy * cy).sqrt();

        let matrix = [
            [a, b * cos_gamma, c * cx],
            [0.0, b * sin_gamma, c * cy],
            [0.0, 0.0, c * cz],
        ];

        Self { matrix }
    }

    /// Return a nalgebra matrix whose **columns** are the lattice vectors.
    ///
    /// Multiplying this matrix by a fractional-coordinate column vector gives
    /// the corresponding Cartesian position:
    ///
    /// ```text
    /// cart = M * frac
    /// ```
    pub fn to_matrix(&self) -> nalgebra::Matrix3<f64> {
        nalgebra::Matrix3::from_columns(&[
            nalgebra::Vector3::new(self.matrix[0][0], self.matrix[0][1], self.matrix[0][2]),
            nalgebra::Vector3::new(self.matrix[1][0], self.matrix[1][1], self.matrix[1][2]),
            nalgebra::Vector3::new(self.matrix[2][0], self.matrix[2][1], self.matrix[2][2]),
        ])
    }

    /// Extract conventional cell parameters (a, b, c, α, β, γ) from the
    /// raw lattice matrix.
    pub fn to_cell_parameters<T, LenUnit>(&self) -> CellParameters<T>
    where
        T: RealField + Copy + Conversion<T, T = T>,
        SI<T>: Units<T>,
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        let mat = self.to_matrix().cast::<T>();

        CellParameters::from_matrix::<LenUnit>(mat)
    }
}

/// A single atom record as it appears in the raw file
pub(crate) struct AtomRecord {
    pub element: Element,
    /// Raw positional values in the file's own length unit
    pub position: [f64; 3],
}

/// The raw output of a format-specific parser
/// The `LenUnit` type parameter selects the unit of length of the file format.
pub(crate) struct RawStructure<LenUnit>
where
    LenUnit: uom::si::length::Unit,
{
    /// `None` for non-periodic formats
    lattice: Option<RawLattice>,
    /// Atom records as they appear in the raw file
    atoms: Vec<AtomRecord>,
    /// The kind of position data (cartesian or fractional)
    position_kind: PositionKind,
    /// Format-specific metadata (e.g. structure index, energy, Q-value)
    /// stored as key-value string pairs.
    properties: HashMap<String, String>,
    /// The unit of length used in the file format
    _unit: PhantomData<LenUnit>,
}

#[derive(Debug, Error)]
pub enum BuildError {
    #[error("fractional without lattice")]
    FractionalWithoutLattice,
}

impl<LenUnit> RawStructure<LenUnit>
where
    LenUnit: uom::si::length::Unit,
{
    /// Assemble a canonical [`Structure<T>`] from the raw parsed data.
    ///
    /// The type parameter `T` selects the floating-point precision of the
    /// output
    ///
    /// * Returns [`BuildError::FractionalWithoutLattice`] when
    ///   `position_kind` is `Fractional` but no lattice was provided.
    pub fn build<T>(self) -> Result<Structure<T>, BuildError>
    where
        T: RealField + Copy + Conversion<T, T = T>,
        SI<T>: Units<T>,
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        // Extract fields before `self.atoms` is consumed by `into_iter`,
        // avoiding a partial-move borrow-check error inside the closure.
        let position_kind = self.position_kind;
        let lattice = self.lattice;
        let properties = self.properties;

        let cartesian_atoms = self
            .atoms
            .into_iter()
            .map(|record| {
                let cart = match position_kind {
                    PositionKind::Cartesian => {
                        let res = nalgebra::Vector3::new(
                            record.position[0],
                            record.position[1],
                            record.position[2],
                        )
                        .cast::<T>();
                        res
                    }
                    PositionKind::Fractional => {
                        let lat = lattice
                            .as_ref()
                            .ok_or(BuildError::FractionalWithoutLattice)?
                            .to_matrix();
                        let frac = nalgebra::Vector3::new(
                            record.position[0],
                            record.position[1],
                            record.position[2],
                        );
                        let res = (lat * frac).cast::<T>();
                        res
                    }
                };

                Ok(Atom {
                    element: record.element,
                    position: CoordinateSystem::Cartesian(Cartesian::from_vector::<LenUnit>(cart)),
                })
            })
            .collect::<Result<Vec<_>, BuildError>>()?;

        // The borrow of `lattice` inside the closure has ended; we can now
        // move it to produce the optional cell parameters.
        let cell = lattice.map(|l| l.to_cell_parameters::<T, LenUnit>());

        Ok(Structure {
            atoms: cartesian_atoms,
            cell,
            properties,
        })
    }
}

/// A trait for parsing structures from text input.
/// The `LenUnit` type parameter is the unit of length used for the parsed structures.
pub(crate) trait RawParser<LenUnit>
where
    LenUnit: uom::si::length::Unit,
{
    /// Parse one or more structures from the input text.
    /// Returns a `Vec` to naturally support multi-frame formats
    fn parse_raw(&self, input: &str) -> Result<Vec<RawStructure<LenUnit>>, String>;
}

#[derive(Debug, Error)]
pub enum ParseStructureError {
    /// The format-specific parser rejected the input
    #[error("format error: {0}")]
    Format(String),
    /// The raw data was syntactically valid but could not be parsed
    #[error("build error: {0}")]
    Build(BuildError),
}

pub trait StructureParser<LenUnit>
where
    LenUnit: uom::si::length::Unit,
{
    fn parse<T>(&self, input: &str) -> Result<Vec<Structure<T>>, ParseStructureError>
    where
        T: RealField + Copy + Conversion<T, T = T>,
        SI<T>: Units<T>,
        LenUnit: Conversion<T, T = T>,
        radian: Conversion<T, T = T>;
}

impl<P, LenUnit> StructureParser<LenUnit> for P
where
    P: RawParser<LenUnit>,
    LenUnit: uom::si::length::Unit,
{
    fn parse<T>(&self, input: &str) -> Result<Vec<Structure<T>>, ParseStructureError>
    where
        T: RealField + Copy + Conversion<T, T = T>,
        SI<T>: Units<T>,
        LenUnit: Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        self.parse_raw(input)
            .map_err(ParseStructureError::Format)?
            .into_iter()
            .map(|raw| raw.build::<T>().map_err(ParseStructureError::Build))
            .collect()
    }
}

pub use arc::{ArcParseError, ArcParser};
