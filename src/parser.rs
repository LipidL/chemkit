mod arc;

use std::collections::HashMap;
use std::error::Error;

use nalgebra::RealField;
use uom::si::angle::{Angle, radian};
use uom::si::length::{Length, angstrom};
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
pub enum PositionKind {
    /// Atom positions are given in Cartesian coordinates
    Cartesian,
    /// Atom positions are given in fractional coordinates
    Fractional,
}

/// Raw lattice – each row of `matrix` is a lattice vector
pub struct RawLattice {
    pub matrix: [[f64; 3]; 3],
}

impl RawLattice {
    pub fn new(matrix: [[f64; 3]; 3]) -> Self {
        Self { matrix }
    }

    /// Construct a `RawLattice` from conventional cell parameters.
    ///
    /// `a`, `b`, `c` are in the caller's length unit; `alpha`, `beta`,
    /// `gamma` are in **degrees**.
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
            [a, 0.0, 0.0],
            [b * cos_gamma, b * sin_gamma, 0.0],
            [c * cx, c * cy, c * cz],
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
        let a = self.matrix[0];
        let b = self.matrix[1];
        let c = self.matrix[2];
        // Matrix3::new arguments are in row-major (m_row_col) order.
        // We want lattice vectors a, b, c to be the three *columns*.
        nalgebra::Matrix3::new(a[0], b[0], c[0], a[1], b[1], c[1], a[2], b[2], c[2])
    }

    /// Generic version of [`to_matrix`]: returns a `Matrix3<T>` by casting each
    /// stored `f64` element to `T` via [`ComplexField::from_f64`].
    ///
    /// This is used by the generic [`RawStructure::build`] to avoid mixing
    /// `f64` matrix arithmetic with the caller's chosen precision `T`.
    pub fn to_matrix_as<T>(&self) -> nalgebra::Matrix3<T>
    where
        T: RealField + Copy,
    {
        let conv =
            |v: f64| T::from_f64(v).expect("RawLattice::to_matrix_as: f64 → T conversion failed");
        let a = self.matrix[0];
        let b = self.matrix[1];
        let c = self.matrix[2];
        nalgebra::Matrix3::new(
            conv(a[0]),
            conv(b[0]),
            conv(c[0]),
            conv(a[1]),
            conv(b[1]),
            conv(c[1]),
            conv(a[2]),
            conv(b[2]),
            conv(c[2]),
        )
    }

    /// Extract conventional cell parameters (a, b, c, α, β, γ) from the
    /// raw lattice matrix.
    ///
    /// `scale` converts raw matrix values to ångströms (pass `1.0` when the
    /// matrix is already in ångströms).
    pub fn to_cell_parameters(&self, scale: f64) -> CellParameters<f64> {
        let va = nalgebra::Vector3::new(
            self.matrix[0][0] * scale,
            self.matrix[0][1] * scale,
            self.matrix[0][2] * scale,
        );
        let vb = nalgebra::Vector3::new(
            self.matrix[1][0] * scale,
            self.matrix[1][1] * scale,
            self.matrix[1][2] * scale,
        );
        let vc = nalgebra::Vector3::new(
            self.matrix[2][0] * scale,
            self.matrix[2][1] * scale,
            self.matrix[2][2] * scale,
        );

        let a = va.norm();
        let b = vb.norm();
        let c = vc.norm();

        let alpha = (vb.dot(&vc) / (b * c)).acos();
        let beta = (va.dot(&vc) / (a * c)).acos();
        let gamma = (va.dot(&vb) / (a * b)).acos();

        CellParameters::new(
            Length::new::<angstrom>(a),
            Length::new::<angstrom>(b),
            Length::new::<angstrom>(c),
            Angle::new::<radian>(alpha),
            Angle::new::<radian>(beta),
            Angle::new::<radian>(gamma),
        )
    }

    /// Generic version of [`to_cell_parameters`]: all arithmetic is performed
    /// in `T`, so the returned [`CellParameters<T>`] matches the precision of
    /// the caller's chosen output type.
    ///
    /// `scale` must already be of type `T`; obtain it by converting the f64
    /// value from [`uom`] with [`ComplexField::from_f64`].
    pub fn to_cell_parameters_as<T>(&self, scale: T) -> CellParameters<T>
    where
        T: RealField + Copy + Conversion<T, T = T>,
        SI<T>: Units<T>,
        angstrom: Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        let conv = |v: f64| {
            T::from_f64(v).expect("RawLattice::to_cell_parameters_as: f64 → T conversion failed")
        };

        let va = nalgebra::Vector3::new(
            conv(self.matrix[0][0]) * scale,
            conv(self.matrix[0][1]) * scale,
            conv(self.matrix[0][2]) * scale,
        );
        let vb = nalgebra::Vector3::new(
            conv(self.matrix[1][0]) * scale,
            conv(self.matrix[1][1]) * scale,
            conv(self.matrix[1][2]) * scale,
        );
        let vc = nalgebra::Vector3::new(
            conv(self.matrix[2][0]) * scale,
            conv(self.matrix[2][1]) * scale,
            conv(self.matrix[2][2]) * scale,
        );

        let a = va.norm();
        let b = vb.norm();
        let c = vc.norm();

        let alpha = (vb.dot(&vc) / (b * c)).acos();
        let beta = (va.dot(&vc) / (a * c)).acos();
        let gamma = (va.dot(&vb) / (a * b)).acos();

        CellParameters::new(
            Length::new::<angstrom>(a),
            Length::new::<angstrom>(b),
            Length::new::<angstrom>(c),
            Angle::new::<radian>(alpha),
            Angle::new::<radian>(beta),
            Angle::new::<radian>(gamma),
        )
    }
}

/// A single atom record as it appears in the raw file
pub struct AtomRecord {
    pub element: Element,
    /// Raw positional values in the file's own length unit
    pub position: [f64; 3],
}

/// The raw output of a format-specific parser
pub struct RawStructure {
    /// `None` for non-periodic formats
    pub lattice: Option<RawLattice>,
    pub atoms: Vec<AtomRecord>,
    pub position_kind: PositionKind,
    /// The unit of length used in the raw file.
    /// A raw positional value `x` represents `x * length_unit` in SI.
    ///
    /// Stored as `f64` because raw file data is always parsed at full
    /// precision regardless of the output type `T` chosen in [`build`].
    pub length_unit: Length<SI<f64>, f64>,
    /// Format-specific metadata (e.g. structure index, energy, Q-value)
    /// stored as key-value string pairs.
    pub properties: HashMap<String, String>,
}

#[derive(Debug)]
pub enum BuildError {
    FractionalWithoutLattice,
}

impl RawStructure {
    /// Assemble a canonical [`Structure<T>`] from the raw parsed data.
    ///
    /// The type parameter `T` selects the floating-point precision of the
    /// output:
    ///
    /// ```text
    /// raw.build::<f32>()  →  Structure<f32>
    /// raw.build::<f64>()  →  Structure<f64>
    /// ```
    ///
    /// Raw file values (stored internally as `f64`) are cast to `T` via
    /// [`ComplexField::from_f64`], which performs a lossless round-trip for
    /// `f64` and a precision-reducing but infallible cast for `f32`.
    ///
    /// * Cartesian positions are scaled from file units to ångströms.
    /// * Fractional positions are converted to Cartesian via the lattice
    ///   matrix, then scaled to angstroms.
    /// * Returns [`BuildError::FractionalWithoutLattice`] when
    ///   `position_kind` is `Fractional` but no lattice was provided.
    pub fn build<T>(self) -> Result<Structure<T>, BuildError>
    where
        T: RealField + Copy + Conversion<T, T = T>,
        SI<T>: Units<T>,
        angstrom: Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        // Obtain the scale factor in f64 from uom, then cast to T.
        let scale_f64 = self.length_unit.get::<angstrom>();
        let scale = T::from_f64(scale_f64)
            .expect("RawStructure::build: length_unit scale → T conversion failed");

        // Extract fields before `self.atoms` is consumed by `into_iter`,
        // avoiding a partial-move borrow-check error inside the closure.
        let position_kind = self.position_kind;
        let lattice = self.lattice;
        let properties = self.properties;

        let cartesian_atoms = self
            .atoms
            .into_iter()
            .map(|record| {
                // Closure that casts a raw f64 position value to T.
                let conv = |v: f64| {
                    T::from_f64(v)
                        .expect("RawStructure::build: position value → T conversion failed")
                };

                let cart = match position_kind {
                    PositionKind::Cartesian => nalgebra::Vector3::new(
                        conv(record.position[0]) * scale,
                        conv(record.position[1]) * scale,
                        conv(record.position[2]) * scale,
                    ),
                    PositionKind::Fractional => {
                        let lat = lattice
                            .as_ref()
                            .ok_or(BuildError::FractionalWithoutLattice)?;
                        let frac = nalgebra::Vector3::new(
                            conv(record.position[0]),
                            conv(record.position[1]),
                            conv(record.position[2]),
                        );
                        // M * frac gives Cartesian in file units; scale to Å.
                        lat.to_matrix_as::<T>() * frac * scale
                    }
                };

                Ok(Atom {
                    element: record.element,
                    position: CoordinateSystem::Cartesian(Cartesian::from_vector(cart)),
                })
            })
            .collect::<Result<Vec<_>, BuildError>>()?;

        // The borrow of `lattice` inside the closure has ended; we can now
        // move it to produce the optional cell parameters.
        let cell = lattice.map(|l| l.to_cell_parameters_as::<T>(scale));

        Ok(Structure {
            atoms: cartesian_atoms,
            cell,
            properties,
        })
    }
}

pub trait StructParser {
    type ParseError: Error;

    /// Parse one or more structures from the input text.
    /// Returns a `Vec` to naturally support multi-frame formats
    fn parse(&self, input: &str) -> Result<Vec<RawStructure>, Self::ParseError>;
}

pub use arc::{ArcParseError, ArcParser};
