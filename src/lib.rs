pub mod atom;
pub mod cell;
pub mod coordinate;
pub mod parser;
pub mod periodic_table;
pub mod structure;
pub mod writer;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::writer::WriteAll;

    use uom::si::length::angstrom;

    // Tests that a structure can be read from arc and write to arc file without significant loss
    #[test]
    fn test_arc_roundtrip() {
        use parser::RawParser;
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
        let parser = parser::ArcParser::new();
        let structures = parser
            .parse_raw(SAMPLE)
            .unwrap()
            .into_iter()
            .map(|s| s.build::<f64>().unwrap())
            .collect::<Vec<_>>();
        let writer = writer::ArcWriter::new();
        let written = writer.write_all(&structures).unwrap();
        let parsed = parser
            .parse_raw(&written)
            .unwrap()
            .into_iter()
            .map(|s| s.build::<f64>().unwrap())
            .collect::<Vec<_>>();
        assert_eq!(
            parsed.len(),
            structures.len(),
            "parsed arc should have same number of structures"
        );
        // Check for all structures
        let tolerance = 1e-6;
        for (parsed, expected) in parsed.iter().zip(structures.iter()) {
            // Check for cell
            let parsed_cell = parsed.cell.as_ref().unwrap().to_matrix::<angstrom>();
            let expected_cell = expected.cell.as_ref().unwrap().to_matrix::<angstrom>();
            assert!(
                (parsed_cell - expected_cell).abs().max() < tolerance,
                "parsed cell should match expected cell"
            );

            // Check for atoms
            for (parsed_atom, expected_atom) in parsed.atoms.iter().zip(expected.atoms.iter()) {
                let parsed_pos = match &parsed_atom.position {
                    coordinate::CoordinateSystem::Cartesian(c) => c.to_vector::<angstrom>(),
                    coordinate::CoordinateSystem::Spherical(s) => s.to_vector::<angstrom>(),
                };
                let expected_pos = match &expected_atom.position {
                    coordinate::CoordinateSystem::Cartesian(c) => c.to_vector::<angstrom>(),
                    coordinate::CoordinateSystem::Spherical(s) => s.to_vector::<angstrom>(),
                };
                assert!(
                    (parsed_pos - expected_pos).abs().max() < tolerance,
                    "parsed atom position should match expected atom position"
                );
            }
        }
    }
}
