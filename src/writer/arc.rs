use std::convert::Infallible;
use std::fmt::Write;

use nalgebra::RealField;
use uom::Conversion;
use uom::si::angle::radian;
use uom::si::length::angstrom;
use uom::si::{SI, Units};

use super::{WriteAll, WriteOne};
use crate::coordinate::CoordinateSystem;
use crate::periodic_table::PeriodicTable;
use crate::structure::Structure;

/// A writer for the Biosym/Materials Studio ARC multi-frame archive format.
pub struct ArcWriter;

impl ArcWriter {
    pub fn new() -> Self {
        Self
    }
}

impl WriteOne for ArcWriter {
    type WriteError = Infallible;

    /// Render a single structure as one ARC frame block (Energy + PBC + atoms + end×2).
    ///
    /// No file-level header is emitted; use [`WriteAll`] to obtain a
    /// complete, self-contained ARC file with the `!BIOSYM archive 2` header.
    fn write_one<T>(&self, structure: &Structure<T>) -> Result<String, Self::WriteError>
    where
        T: RealField + Copy + Conversion<T, T = T> + std::fmt::Display,
        SI<T>: Units<T>,
        angstrom: Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        write_structure_block(structure)
    }
}

impl WriteAll for ArcWriter {
    type WriteError = Infallible;

    /// Render a slice of structures as a complete ARC file.
    ///
    /// Emits the `!BIOSYM archive 2` / `PBC=ON` file header, then appends
    /// one frame block per structure in order.
    fn write_all<T>(&self, structures: &[Structure<T>]) -> Result<String, Self::WriteError>
    where
        T: RealField + Copy + Conversion<T, T = T> + std::fmt::Display,
        SI<T>: Units<T>,
        angstrom: Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        let mut out = String::from("!BIOSYM archive 2\nPBC=ON\n");
        for structure in structures {
            out.push_str(&write_structure_block(structure)?);
        }
        Ok(out)
    }
}

// ── Serialisation helper ──────────────────────────────────────────────────────

/// Render a single structure block:
/// `Energy` line → optional `PBC` line → atom records → `end`×2.
///
/// This function is shared by both [`WriteOne`] and [`WriteAll`].
fn write_structure_block<T>(structure: &Structure<T>) -> Result<String, Infallible>
where
    T: RealField + Copy + Conversion<T, T = T> + std::fmt::Display,
    SI<T>: Units<T>,
    angstrom: Conversion<T, T = T>,
    radian: Conversion<T, T = T>,
{
    let table = PeriodicTable::new();

    // Precompute radians → degrees factor in type T.
    let pi = T::from_f64(std::f64::consts::PI).expect("π → T conversion must not fail for f32/f64");
    let deg_per_rad = T::from_f64(180.0).expect("180 → T conversion must not fail") / pi;

    // Pull ARC-specific metadata from the properties map; fall back to safe
    // defaults when the structure was not originally read from an ARC file.
    let index: u64 = structure
        .properties
        .get("arc_index")
        .and_then(|s| s.parse().ok())
        .unwrap_or(0);
    let q_value: &str = structure
        .properties
        .get("arc_q_value")
        .map(String::as_str)
        .unwrap_or("0.0");
    let energy: &str = structure
        .properties
        .get("arc_energy")
        .map(String::as_str)
        .unwrap_or("0.0");

    let mut out = String::new();

    // ── Energy header line ────────────────────────────────────────────────
    // Biosym convention: 22 leading spaces, then "Energy <index> <q> <E>".
    // The parser only requires "Energy" as the first whitespace-separated
    // token followed by three numeric values, so exact column widths are
    // flexible.
    writeln!(
        out,
        "                      Energy {:>9} {:>14} {:>16}",
        index, q_value, energy
    )
    .unwrap();

    // ── PBC line (cell parameters in Å and degrees) ───────────────────────
    // Omitted when the structure has no periodic cell (molecule-in-a-box).
    if let Some(cell) = &structure.cell {
        let a = cell.a.get::<angstrom>();
        let b = cell.b.get::<angstrom>();
        let c = cell.c.get::<angstrom>();
        let alpha = cell.alpha.get::<radian>() * deg_per_rad;
        let beta = cell.beta.get::<radian>() * deg_per_rad;
        let gamma = cell.gamma.get::<radian>() * deg_per_rad;

        // Format: "PBC" followed by six right-justified 14.8f fields.
        // Matches the column layout produced by Materials Studio.
        writeln!(
            out,
            "PBC{:>14.8}{:>14.8}{:>14.8}{:>14.8}{:>14.8}{:>14.8}",
            a, b, c, alpha, beta, gamma
        )
        .unwrap();
    }

    // Atom records
    // Layout (Biosym atom-record format):
    //   %-9s %14.9f %14.9f %14.9f CORE %4d %-2s %-2s %10.4f %4d
    for (idx, atom) in structure.atoms.iter().enumerate() {
        let serial = idx + 1; // 1-based

        // Resolve the symbol via a reverse atomic-number lookup.  The
        // PeriodicTable stores elements by symbol string, so this is an O(118)
        // linear scan – perfectly acceptable for typical structure sizes.
        let symbol = table.symbol_of(atom.element.atomic_number).unwrap_or("?");

        // Accept both Cartesian and Spherical storage; convert to a raw
        // Cartesian vector in ångströms for output.
        let vec = match &atom.position {
            CoordinateSystem::Cartesian(c) => c.to_vector::<angstrom>(),
            CoordinateSystem::Spherical(s) => s.to_vector::<angstrom>(),
        };
        let x = vec.x;
        let y = vec.y;
        let z = vec.z;

        writeln!(
            out,
            "{:<9}{:>14.9}{:>14.9}{:>14.9} CORE {:>4} {:<2} {:<2} {:>10.4} {:>4}",
            symbol, x, y, z, serial, symbol, symbol, 0.0_f64, serial
        )
        .unwrap();
    }

    out.push_str("end\nend\n");

    Ok(out)
}
