use std::collections::HashMap;

use thiserror::Error;
use uom::si::f64::Length;
use uom::si::length::angstrom;

use super::{AtomRecord, PositionKind, RawLattice, RawStructure, StructParser};
use crate::periodic_table::PeriodicTable;

// Error type

/// The only hard error the lenient state-machine parser can produce:
/// a line was unambiguously identified as an atom line (field 4 == `"CORE"`)
/// but the element symbol in field 0 is not in the periodic table.
#[derive(Error, Debug)]
pub enum ArcParseError {
    #[error("line {line}: unknown element '{symbol}'")]
    UnknownElement { line: usize, symbol: String },
}

// Parser

pub struct ArcParser;

impl ArcParser {
    pub fn new() -> Self {
        Self
    }
}

impl StructParser for ArcParser {
    type ParseError = ArcParseError;

    /// Parse all complete structures from an ARC file using a four-state machine.
    ///
    /// ```text
    /// SeekEnergy ──(Energy line)──> SeekPbc
    ///            <─(any other)────
    ///
    /// SeekPbc ──(PBC line)──> SeekAtoms
    ///         <─(any other)──
    ///
    /// SeekAtoms ──(atom line)──> SeekAtoms   (accumulate)
    ///           ──("end")──────> SeekSecondEnd
    ///           <─(any other)───
    ///
    /// SeekSecondEnd ──("end")──> SeekEnergy  (emit structure)
    ///               <─(any other)───
    /// ```
    ///
    /// Unrecognised lines are silently skipped at every transition.
    /// The only hard error is encountering a definite atom line whose element
    /// symbol is not in the periodic table.
    fn parse(&self, input: &str) -> Result<Vec<RawStructure>, Self::ParseError> {
        let table = PeriodicTable::new();
        let mut structures = Vec::new();

        enum State {
            SeekEnergy,
            SeekPbc {
                index: u64,
                q_value: f64,
                energy: f64,
            },
            SeekAtoms {
                index: u64,
                q_value: f64,
                energy: f64,
                lattice: RawLattice,
                atoms: Vec<AtomRecord>,
            },
            SeekSecondEnd {
                index: u64,
                q_value: f64,
                energy: f64,
                lattice: RawLattice,
                atoms: Vec<AtomRecord>,
            },
        }

        let mut state = State::SeekEnergy;

        for (line_no, line) in input.lines().enumerate().map(|(i, l)| (i + 1, l)) {
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }

            state = match state {
                // Looking for the start of a structure block
                State::SeekEnergy => {
                    if let Some((index, q_value, energy)) = try_parse_energy(trimmed) {
                        State::SeekPbc {
                            index,
                            q_value,
                            energy,
                        }
                    } else {
                        State::SeekEnergy
                    }
                }

                // Looking for the cell parameters
                State::SeekPbc {
                    index,
                    q_value,
                    energy,
                } => {
                    if let Some(lattice) = try_parse_pbc(trimmed) {
                        State::SeekAtoms {
                            index,
                            q_value,
                            energy,
                            lattice,
                            atoms: Vec::new(),
                        }
                    } else {
                        State::SeekPbc {
                            index,
                            q_value,
                            energy,
                        }
                    }
                }

                // Accumulating atom lines
                State::SeekAtoms {
                    index,
                    q_value,
                    energy,
                    lattice,
                    mut atoms,
                } => {
                    if trimmed == "end" {
                        State::SeekSecondEnd {
                            index,
                            q_value,
                            energy,
                            lattice,
                            atoms,
                        }
                    } else if let Some(record) = try_parse_atom(trimmed, &table, line_no)? {
                        atoms.push(record);
                        State::SeekAtoms {
                            index,
                            q_value,
                            energy,
                            lattice,
                            atoms,
                        }
                    } else {
                        // Unrecognised line — skip
                        State::SeekAtoms {
                            index,
                            q_value,
                            energy,
                            lattice,
                            atoms,
                        }
                    }
                }

                // Waiting for the closing "end" to complete the block
                State::SeekSecondEnd {
                    index,
                    q_value,
                    energy,
                    lattice,
                    atoms,
                } => {
                    if trimmed == "end" {
                        let mut properties = HashMap::new();
                        properties.insert("arc_index".to_string(), index.to_string());
                        properties.insert("arc_q_value".to_string(), q_value.to_string());
                        properties.insert("arc_energy".to_string(), energy.to_string());

                        structures.push(RawStructure {
                            lattice: Some(lattice),
                            atoms,
                            position_kind: PositionKind::Cartesian,
                            length_unit: Length::new::<angstrom>(1.0),
                            properties,
                        });

                        State::SeekEnergy
                    } else {
                        // Unrecognised line between the two "end"s — skip
                        State::SeekSecondEnd {
                            index,
                            q_value,
                            energy,
                            lattice,
                            atoms,
                        }
                    }
                }
            };
        }

        Ok(structures)
    }
}

// Line recognisers

/// Recognise an `Energy <index> <Q> <energy>` line.
/// Returns `None` for any line that does not match.
fn try_parse_energy(line: &str) -> Option<(u64, f64, f64)> {
    let mut fields = line.split_whitespace();
    if fields.next()? != "Energy" {
        return None;
    }
    let index = fields.next()?.parse().ok()?;
    let q_value = fields.next()?.parse().ok()?;
    let energy = fields.next()?.parse().ok()?;
    Some((index, q_value, energy))
}

/// Recognise a `PBC <a> <b> <c> <α> <β> <γ>` line.
/// Returns `None` for any line that does not match.
fn try_parse_pbc(line: &str) -> Option<RawLattice> {
    let mut fields = line.split_whitespace();
    if fields.next()? != "PBC" {
        return None;
    }
    let a: f64 = fields.next()?.parse().ok()?;
    let b: f64 = fields.next()?.parse().ok()?;
    let c: f64 = fields.next()?.parse().ok()?;
    let alpha: f64 = fields.next()?.parse().ok()?;
    let beta: f64 = fields.next()?.parse().ok()?;
    let gamma: f64 = fields.next()?.parse().ok()?;
    Some(RawLattice::from_cell_parameters(
        a, b, c, alpha, beta, gamma,
    ))
}

/// Recognise an atom line, identified by `"CORE"` at whitespace-field index 4.
///
/// Returns:
/// * `Ok(Some(record))` — valid atom line, element found in the periodic table
/// * `Ok(None)` — not an atom line (or coordinates are malformed); skip it
/// * `Err(UnknownElement)` — line is definitely an atom line but the element
///   symbol at field 0 is unknown — this is the only hard parse error
fn try_parse_atom(
    line: &str,
    table: &PeriodicTable,
    line_no: usize,
) -> Result<Option<AtomRecord>, ArcParseError> {
    let fields: Vec<&str> = line.split_whitespace().collect();

    // Use the "CORE" sentinel at index 4 to identify atom lines.
    if fields.len() < 10 || fields[4] != "CORE" {
        return Ok(None);
    }

    let symbol = fields[0];
    let element = table
        .get_by_symbol(symbol)
        .ok_or_else(|| ArcParseError::UnknownElement {
            line: line_no,
            symbol: symbol.to_string(),
        })?;

    // If coordinates are somehow malformed, treat the line as unrecognised.
    let (x, y, z) = match (
        fields[1].parse::<f64>(),
        fields[2].parse::<f64>(),
        fields[3].parse::<f64>(),
    ) {
        (Ok(x), Ok(y), Ok(z)) => (x, y, z),
        _ => return Ok(None),
    };

    Ok(Some(AtomRecord {
        element,
        position: [x, y, z],
    }))
}

// Tests

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::StructParser;

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

    #[test]
    fn test_parse_structure_count() {
        let parser = ArcParser;
        let structures = parser.parse(SAMPLE).expect("parse failed");
        assert_eq!(structures.len(), 2, "expected 2 structures");
    }

    #[test]
    fn test_parse_atom_count() {
        let parser = ArcParser;
        let structures = parser.parse(SAMPLE).expect("parse failed");
        assert_eq!(
            structures[0].atoms.len(),
            2,
            "first structure should have 2 atoms"
        );
        assert_eq!(
            structures[1].atoms.len(),
            2,
            "second structure should have 2 atoms"
        );
    }

    #[test]
    fn test_parse_first_atom_position() {
        let parser = ArcParser;
        let structures = parser.parse(SAMPLE).expect("parse failed");
        let pos = &structures[0].atoms[0].position;
        assert!((pos[0] - 6.066246000_f64).abs() < 1e-9);
        assert!((pos[1] - 63.544500000_f64).abs() < 1e-9);
        assert!((pos[2] - 3.058141000_f64).abs() < 1e-9);
    }

    #[test]
    fn test_parse_lattice() {
        let parser = ArcParser;
        let structures = parser.parse(SAMPLE).expect("parse failed");
        let lattice = structures[0].lattice.as_ref().expect("expected a lattice");
        // The a-vector is aligned with x; its norm should equal a = 54.052 Å.
        let a_mag = (lattice.matrix[0][0].powi(2)
            + lattice.matrix[0][1].powi(2)
            + lattice.matrix[0][2].powi(2))
        .sqrt();
        assert!((a_mag - 54.052_f64).abs() < 1e-4, "a = {a_mag}");
    }

    #[test]
    fn test_parse_properties() {
        let parser = ArcParser;
        let structures = parser.parse(SAMPLE).expect("parse failed");

        let p0 = &structures[0].properties;
        assert_eq!(p0.get("arc_index").map(String::as_str), Some("0"));
        assert_eq!(p0.get("arc_q_value").map(String::as_str), Some("3.8755"));
        assert_eq!(
            p0.get("arc_energy").map(String::as_str),
            Some("-35146.40625")
        );

        let p1 = &structures[1].properties;
        assert_eq!(p1.get("arc_index").map(String::as_str), Some("100"));
    }

    #[test]
    fn test_position_kind_is_cartesian() {
        let parser = ArcParser;
        let structures = parser.parse(SAMPLE).expect("parse failed");
        assert!(matches!(
            structures[0].position_kind,
            PositionKind::Cartesian
        ));
    }

    #[test]
    fn test_no_header() {
        // A file with no !BIOSYM / PBC=ON header should parse normally.
        let input = "\
                      Energy         0          0.0     -100.0
PBC   10.0   10.0   10.0   90.0   90.0   90.0
H        1.0   2.0   3.0 CORE    1 H  H    0.0000    1
end
end";
        let parser = ArcParser;
        let structures = parser.parse(input).expect("parse failed");
        assert_eq!(structures.len(), 1);
        assert_eq!(structures[0].atoms.len(), 1);
    }

    #[test]
    fn test_error_unknown_element() {
        let input = "\
                      Energy         0          0.0     -100.0
PBC   10.0   10.0   10.0   90.0   90.0   90.0
Xx        1.0   2.0   3.0 CORE    1 Xx  Xx    0.0000    1
end
end";
        let parser = ArcParser;
        let result = parser.parse(input);
        assert!(matches!(result, Err(ArcParseError::UnknownElement { .. })));
    }

    #[test]
    fn test_truncated_file_returns_empty() {
        // A truncated file is not an error — the state machine simply never
        // completes a structure block, so nothing is emitted.
        let input = "\
                      Energy         0          0.0     -100.0
PBC   10.0   10.0   10.0   90.0   90.0   90.0";
        let parser = ArcParser;
        let result = parser
            .parse(input)
            .expect("truncated input should not error");
        assert!(result.is_empty());
    }

    #[test]
    fn test_interleaved_noise_is_skipped() {
        // Unrecognised lines appearing between any two sections are silently
        // discarded by the state machine.
        let input = "\
Some random header line
Another line that is not a keyword
                      Energy         0          0.0     -100.0
This line is between Energy and PBC and should be skipped
PBC   10.0   10.0   10.0   90.0   90.0   90.0
H        1.0   2.0   3.0 CORE    1 H  H    0.0000    1
end
A line squeezed between the two end markers
end";
        let parser = ArcParser;
        let structures = parser.parse(input).expect("parse failed");
        assert_eq!(structures.len(), 1);
        assert_eq!(structures[0].atoms.len(), 1);
    }
}
