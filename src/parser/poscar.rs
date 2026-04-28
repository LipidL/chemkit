use std::{collections::HashMap, marker::PhantomData};

use nalgebra::{Matrix3, Vector3};
use uom::si::length::angstrom;

use super::{AtomRecord, PositionKind, RawLattice, RawParser, RawStructure};
use crate::periodic_table::{Element, PeriodicTable};

// Parser
pub struct PoscarParser;

impl PoscarParser {
    pub fn new() -> Self {
        Self
    }
}

/// The mode of scaling factors used to specify the lattice.
enum ScalingMode {
    /// The scaling factor specifies the desired cell volume
    CellVolume,
    /// The scaling factor specifies the factors directly
    Direct,
}

impl ScalingMode {
    /// Parses a trimmed scaling factor string into a [`ScalingMode`].
    fn from_str(s: &str) -> Self {
        if s.starts_with('-') {
            Self::CellVolume
        } else {
            Self::Direct
        }
    }
}

fn parse_potcar(input: &str, table: PeriodicTable) -> Vec<Element> {
    input
        .lines()
        .into_iter()
        .filter_map(|s| {
            s.trim()
                .strip_prefix("TITEL")?
                .trim()
                .strip_prefix("=")?
                .split(&[' ', '_'])
                .find_map(|s| table.get_by_symbol(s))
        })
        .collect()
}

impl RawParser<angstrom, (&str, Option<&str>)> for PoscarParser {
    fn parse_raw(
        &self,
        input: (&str, Option<&str>),
    ) -> Result<Vec<RawStructure<angstrom>>, String> {
        let poscar_input = input.0;
        let potcar_input = input.1;
        let mut properties = HashMap::new();
        let table = PeriodicTable::new();
        let mut poscar_iter = poscar_input.lines().into_iter();

        // Comment line
        let comment = poscar_iter.next().ok_or("No comment line")?;
        properties.insert("comment".to_string(), comment.to_string());

        // Scaling factor
        let scaling_factor_line = poscar_iter
            .next()
            .ok_or("No scaling factor line provided")?
            .trim();
        let scaling_mode = ScalingMode::from_str(scaling_factor_line);

        // Lattice (3 lines)
        let mut lattice = Vec::new();
        for _ in 0..3 {
            let line = poscar_iter.next().ok_or("No lattice line")?;
            let row: [f64; 3] = line
                .split_whitespace()
                .map(|s| {
                    s.parse::<f64>()
                        .map_err(|e| format!("Error parsing lattice: {e}"))
                })
                .collect::<Result<Vec<_>, _>>()?
                .try_into()
                .map_err(|v: Vec<_>| {
                    format!("Expected 3 values per lattice line, got {}", v.len())
                })?;
            lattice.push(row);
        }
        let lattice_array: [f64; 9] = lattice
            .into_iter()
            .flatten()
            .collect::<Vec<_>>()
            .try_into()
            .map_err(|v: Vec<_>| format!("Expected 9 values for lattice, got {}", v.len()))?;
        let lattice = Matrix3::from_row_slice(&lattice_array);

        // Calculate scaling factor
        let (lattice, scaling_factor) = match scaling_mode {
            ScalingMode::Direct => {
                let scaling_factor = scaling_factor_line
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<f64>()
                            .map_err(|e| format!("Error parsing scaling factor: {e}"))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let scaling_factor = match scaling_factor.len() {
                    1 => [scaling_factor[0]; 3],
                    3 => scaling_factor.try_into().unwrap(),
                    _ => return Err("Expected 1 or 3 values for scaling factor".to_string()),
                };
                (lattice, scaling_factor)
            }
            ScalingMode::CellVolume => {
                let expected_volume = scaling_factor_line
                    .parse::<f64>()
                    .map_err(|e| format!("Error parsing scaling factor: {e}"))?
                    .abs(); // The expected volume is given in the scaling factor line with leading "-"
                let actual_volume = lattice.determinant();
                let scaling_factor = expected_volume / actual_volume;
                // Scale the lattice to get the expected volume
                let lattice = lattice * scaling_factor;
                (lattice, [scaling_factor; 3])
            }
        };
        let lattice = RawLattice::from_matrix(lattice);
        let scaling_factor = Vector3::from(scaling_factor);

        // Species names and/or ion number
        let species_names_or_ion_number_line = poscar_iter
            .next()
            .ok_or("No species names or ion number line")?;
        let species: Vec<(Element, usize)> = match species_names_or_ion_number_line
            .split_whitespace()
            .map(|s| {
                table
                    .get_by_symbol(s)
                    .ok_or(format!("Unknown element: {s}"))
            })
            .collect::<Result<Vec<_>, _>>()
        {
            // This poscar file has species names followed by ion numbers
            Ok(species) => {
                let ion_numbers = poscar_iter
                    .next()
                    .ok_or("No ion number line")?
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<usize>()
                            .map_err(|e| format!("Error parsing ion number: {e}"))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                species.into_iter().zip(ion_numbers).collect()
            }
            // This poscar file has only ion numbers line, and the species names should be read from potcar
            Err(_) => {
                let ion_numbers = species_names_or_ion_number_line
                    .split_whitespace()
                    .map(|s| {
                        s.parse::<usize>()
                            .map_err(|e| format!("Error parsing ion number: {e}"))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let potcar_input = potcar_input
                    .ok_or("No species names line in POSCAR while no POTCAR provided")?;
                let elements = parse_potcar(potcar_input, table);
                elements.into_iter().zip(ion_numbers).collect()
            }
        };

        // Coordinate system
        let selective_dynamics_or_coordinate_system_line = poscar_iter
            .next()
            .ok_or("No Selective Dynamics or Coordinate System line")?;
        let coordinate_system_line = if selective_dynamics_or_coordinate_system_line
            .trim()
            .to_lowercase()
            .starts_with("s")
        {
            poscar_iter
                .next()
                .ok_or("No Coordinate System line after Selective Dynamics line")?
        } else {
            selective_dynamics_or_coordinate_system_line
        };
        let position_kind = if coordinate_system_line
            .trim()
            .to_lowercase()
            .starts_with(&['c', 'k'])
        {
            PositionKind::Cartesian
        } else {
            PositionKind::Fractional
        };

        // Read atomic positions
        let num_atoms = species.iter().map(|s| s.1).sum::<usize>();
        let mut atoms = Vec::new();
        fn current_species(species: Vec<(Element, usize)>, i: usize) -> Element {
            let mut total_cnt = 0;
            for (_, spe) in species.iter().enumerate() {
                total_cnt += spe.1;
                if i < total_cnt {
                    return spe.0.clone();
                }
            }
            unreachable!()
        }
        for i in 0..num_atoms {
            let line = poscar_iter
                .next()
                .ok_or("Not enough atomic position lines")?;

            let coords: [f64; 3] = line
                .split_whitespace()
                .filter_map(|s| s.parse().ok())
                .collect::<Vec<_>>()
                .try_into()
                .map_err(|v: Vec<_>| {
                    format! {"Invalid atomic position line: expected 3 values, got {}",
                    v.len()}
                })?;
            let coords = Vector3::from(coords);
            // If the coordination is Cartesian, scale the coordinates
            let coords = match position_kind {
                PositionKind::Cartesian => coords.component_mul(&scaling_factor),
                PositionKind::Fractional => coords,
            };
            // convert back to [f64; 3]
            let coords: [f64; 3] = coords.into();
            let element = current_species(species.clone(), i);
            atoms.push((element, coords));
        }
        let atoms = atoms
            .into_iter()
            .map(|(e, c)| AtomRecord {
                element: e,
                position: c,
            })
            .collect::<Vec<_>>();
        let raw_structure = RawStructure {
            lattice: Some(lattice),
            atoms,
            position_kind,
            properties,
            _unit: PhantomData,
        };
        Ok(vec![raw_structure])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::RawParser;

    fn check_position(expected: &[f64; 3], actual: &[f64; 3], tolerance: f64) {
        assert!(
            (Vector3::from_column_slice(expected) - Vector3::from_column_slice(actual)).norm()
                < tolerance
        );
    }

    const SAMPLE_FRACTIONAL: &str = "\
Comment line
1.0
10.0 0.0 0.0
0.0 10.0 0.0
0.0 0.0 10.0
O H
1 2
Direct
0.0 0.0 0.0
0.0 0.0 0.5
0.0 0.5 0.5";

    #[test]
    fn test_parse_fractional() {
        let parser = PoscarParser;
        let structures = parser
            .parse_raw((SAMPLE_FRACTIONAL, None))
            .expect("parse failed");
        assert_eq!(structures.len(), 1);
        let structure = &structures[0];
        assert_eq!(structure.atoms.len(), 3);
        assert_eq!(structure.atoms[0].element.atomic_number, 8);
        assert_eq!(structure.atoms[1].element.atomic_number, 1);
        assert_eq!(structure.atoms[2].element.atomic_number, 1);
        let expected_frac_positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.5, 0.5]];
        for (expected, actual) in expected_frac_positions.iter().zip(structure.atoms.iter()) {
            check_position(expected, &actual.position, 1e-6);
        }
        assert!(matches!(structure.position_kind, PositionKind::Fractional));
    }
}
