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
            CoordinateSystem::Cartesian(c) => c.to_vector(),
            CoordinateSystem::Spherical(s) => s.to_vector(),
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

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{ArcParser, StructParser};
    use crate::writer::{WriteAll, WriteOne};

    const SAMPLE: &str = "\
!BIOSYM archive 2
PBC=ON
                      Energy         0          3.8755     -35146.406250
!DATE rotcell   0.000000000000000  -0.000000000000000   0.000000000000000
PBC   54.05200000   85.09600000   75.83690000   90.00000000   90.00000000   90.00000000
O        6.066246000   63.544500000    3.058141000 CORE    1 O  O    0.0000    1
O        5.993490000   80.551210000    2.889045000 CORE    2 O  O    0.0000    2
end
end
                      Energy       100          7.1159     -34976.703125
!DATE rotcell   0.000000000000000  -0.000000000000000   0.000000000000000
PBC   54.05200000   85.09600000   75.83690000   90.00000000   90.00000000   90.00000000
O        5.985308177   63.919811755    2.958429785 CORE    1 O  O    0.0000    1
O        6.103548946   80.754955910    2.709105556 CORE    2 O  O    0.0000    2
end
end";

    /// Parse the SAMPLE text and build concrete f64 [`Structure`]s.
    fn parse_sample() -> Vec<crate::structure::Structure<f64>> {
        ArcParser::new()
            .parse(SAMPLE)
            .expect("parse failed")
            .into_iter()
            .map(|raw| raw.build::<f64>().expect("build failed"))
            .collect()
    }

    // WriteOne tests

    #[test]
    fn write_one_round_trips_structure_count() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_one(&structures[0]).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");
        assert_eq!(reparsed.len(), 1, "expected exactly one frame block");
    }

    #[test]
    fn write_one_round_trips_atom_count() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_one(&structures[0]).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");
        assert_eq!(reparsed[0].atoms.len(), 2);
    }

    #[test]
    fn write_one_round_trips_position() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_one(&structures[0]).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");

        // Extract the original x-coordinate via the Cartesian variant.
        let CoordinateSystem::Cartesian(ref c) = structures[0].atoms[0].position else {
            panic!("expected Cartesian storage");
        };
        let orig_x: f64 = c.to_vector().x;
        let rep_x: f64 = reparsed[0].atoms[0].position[0];

        assert!(
            (orig_x - rep_x).abs() < 1e-6,
            "x after round-trip differs: orig={orig_x} reparsed={rep_x}"
        );
    }

    #[test]
    fn write_one_without_cell_omits_pbc_line() {
        use crate::atom::Atom;
        use crate::coordinate::{Cartesian, CoordinateSystem};
        use crate::periodic_table::PeriodicTable;
        use uom::si::f64::Length;
        use uom::si::length::angstrom as ang;

        let table = PeriodicTable::new();
        let structure = crate::structure::Structure {
            atoms: vec![Atom {
                element: table.get_by_symbol("H").unwrap(),
                position: CoordinateSystem::Cartesian(Cartesian::new(
                    Length::new::<ang>(1.0),
                    Length::new::<ang>(2.0),
                    Length::new::<ang>(3.0),
                )),
            }],
            cell: None,
            properties: Default::default(),
        };

        let text = ArcWriter::new().write_one(&structure).unwrap();
        assert!(
            !text.contains("PBC"),
            "no PBC line expected when cell is None"
        );
        assert!(text.contains("CORE"), "atom record must be present");
        assert!(text.contains("end"), "end markers must be present");
    }

    // WriteAll tests

    #[test]
    fn write_all_emits_file_header() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_all(&structures).unwrap();
        assert!(
            text.starts_with("!BIOSYM archive 2\nPBC=ON\n"),
            "expected Biosym file header at start"
        );
    }

    #[test]
    fn write_all_round_trips_structure_count() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_all(&structures).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");
        assert_eq!(reparsed.len(), 2, "expected two frame blocks");
    }

    #[test]
    fn write_all_round_trips_atom_counts() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_all(&structures).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");
        assert_eq!(reparsed[0].atoms.len(), 2);
        assert_eq!(reparsed[1].atoms.len(), 2);
    }

    #[test]
    fn write_all_preserves_cell_parameters() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_all(&structures).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");

        let lattice = reparsed[0].lattice.as_ref().expect("expected a lattice");
        let a_mag = (lattice.matrix[0][0].powi(2)
            + lattice.matrix[0][1].powi(2)
            + lattice.matrix[0][2].powi(2))
        .sqrt();
        assert!(
            (a_mag - 54.052_f64).abs() < 1e-3,
            "a-axis magnitude after round-trip: {a_mag}"
        );
    }

    #[test]
    fn write_all_preserves_arc_index_property() {
        let structures = parse_sample();
        let text = ArcWriter::new().write_all(&structures).unwrap();
        let reparsed = ArcParser::new().parse(&text).expect("re-parse failed");

        assert_eq!(
            reparsed[0].properties.get("arc_index").map(String::as_str),
            Some("0"),
            "first frame arc_index"
        );
        assert_eq!(
            reparsed[1].properties.get("arc_index").map(String::as_str),
            Some("100"),
            "second frame arc_index"
        );
    }
}
