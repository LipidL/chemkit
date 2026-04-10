use std::collections::HashMap;
use uom::si::f64::{Length, Mass};
use uom::si::length::angstrom;
use uom::si::mass::dalton;

#[derive(Clone, Debug)]
pub struct Element {
    pub name: String,
    pub atomic_number: u16,
    pub mass: Mass,
    pub valence_radius: Length,
    pub valence_electrons: u32,
    pub atom_radius: Length,
}

pub struct PeriodicTable {
    elements: HashMap<String, Element>,
}

impl PeriodicTable {
    pub fn new() -> Self {
        let mut elements = HashMap::new();
        // Note: The Element struct needs an ion_radius field to pass tests
        elements.insert(
            "H".to_string(),
            Element {
                name: "Hydrogen".to_string(),
                atomic_number: 1,
                mass: Mass::new::<dalton>(1.008),
                valence_radius: Length::new::<angstrom>(0.32),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(0.25),
            },
        );

        elements.insert(
            "He".to_string(),
            Element {
                name: "Helium".to_string(),
                atomic_number: 2,
                mass: Mass::new::<dalton>(4.003),
                valence_radius: Length::new::<angstrom>(0.93),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(0.31),
            },
        );

        elements.insert(
            "Li".to_string(),
            Element {
                name: "Lithium".to_string(),
                atomic_number: 3,
                mass: Mass::new::<dalton>(6.941),
                valence_radius: Length::new::<angstrom>(1.23),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(1.67),
            },
        );

        elements.insert(
            "Be".to_string(),
            Element {
                name: "Beryllium".to_string(),
                atomic_number: 4,
                mass: Mass::new::<dalton>(9.0122),
                valence_radius: Length::new::<angstrom>(0.90),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(1.12),
            },
        );

        elements.insert(
            "B".to_string(),
            Element {
                name: "Boron".to_string(),
                atomic_number: 5,
                mass: Mass::new::<dalton>(10.811),
                valence_radius: Length::new::<angstrom>(0.82),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(0.87),
            },
        );

        elements.insert(
            "C".to_string(),
            Element {
                name: "Carbon".to_string(),
                atomic_number: 6,
                mass: Mass::new::<dalton>(12.011),
                valence_radius: Length::new::<angstrom>(0.67),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(0.77),
            },
        );

        elements.insert(
            "N".to_string(),
            Element {
                name: "Nitrogen".to_string(),
                atomic_number: 7,
                mass: Mass::new::<dalton>(14.007),
                valence_radius: Length::new::<angstrom>(0.56),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(0.75),
            },
        );

        elements.insert(
            "O".to_string(),
            Element {
                name: "Oxygen".to_string(),
                atomic_number: 8,
                mass: Mass::new::<dalton>(16.00),
                valence_radius: Length::new::<angstrom>(0.60),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(0.73),
            },
        );

        elements.insert(
            "F".to_string(),
            Element {
                name: "Fluorine".to_string(),
                atomic_number: 9,
                mass: Mass::new::<dalton>(18.998),
                valence_radius: Length::new::<angstrom>(0.50),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(0.71),
            },
        );

        elements.insert(
            "Ne".to_string(),
            Element {
                name: "Neon".to_string(),
                atomic_number: 10,
                mass: Mass::new::<dalton>(20.180),
                valence_radius: Length::new::<angstrom>(1.12),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(0.69),
            },
        );

        elements.insert(
            "Na".to_string(),
            Element {
                name: "Sodium".to_string(),
                atomic_number: 11,
                mass: Mass::new::<dalton>(22.990),
                valence_radius: Length::new::<angstrom>(1.54),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(1.90),
            },
        );

        elements.insert(
            "Mg".to_string(),
            Element {
                name: "Magnesium".to_string(),
                atomic_number: 12,
                mass: Mass::new::<dalton>(24.305),
                valence_radius: Length::new::<angstrom>(1.36),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(1.45),
            },
        );

        elements.insert(
            "Al".to_string(),
            Element {
                name: "Aluminum".to_string(),
                atomic_number: 13,
                mass: Mass::new::<dalton>(26.982),
                valence_radius: Length::new::<angstrom>(1.18),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.18),
            },
        );

        elements.insert(
            "Si".to_string(),
            Element {
                name: "Silicon".to_string(),
                atomic_number: 14,
                mass: Mass::new::<dalton>(28.086),
                valence_radius: Length::new::<angstrom>(1.11),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.11),
            },
        );

        elements.insert(
            "P".to_string(),
            Element {
                name: "Phosphorus".to_string(),
                atomic_number: 15,
                mass: Mass::new::<dalton>(30.974),
                valence_radius: Length::new::<angstrom>(1.06),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.06),
            },
        );

        elements.insert(
            "S".to_string(),
            Element {
                name: "Sulfur".to_string(),
                atomic_number: 16,
                mass: Mass::new::<dalton>(32.066),
                valence_radius: Length::new::<angstrom>(1.02),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.02),
            },
        );

        elements.insert(
            "Cl".to_string(),
            Element {
                name: "Chlorine".to_string(),
                atomic_number: 17,
                mass: Mass::new::<dalton>(35.453),
                valence_radius: Length::new::<angstrom>(0.99),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(0.99),
            },
        );

        elements.insert(
            "Ar".to_string(),
            Element {
                name: "Argon".to_string(),
                atomic_number: 18,
                mass: Mass::new::<dalton>(39.948),
                valence_radius: Length::new::<angstrom>(1.54),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(0.97),
            },
        );

        elements.insert(
            "K".to_string(),
            Element {
                name: "Potassium".to_string(),
                atomic_number: 19,
                mass: Mass::new::<dalton>(39.098),
                valence_radius: Length::new::<angstrom>(2.03),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(2.43),
            },
        );

        elements.insert(
            "Ca".to_string(),
            Element {
                name: "Calcium".to_string(),
                atomic_number: 20,
                mass: Mass::new::<dalton>(40.078),
                valence_radius: Length::new::<angstrom>(1.74),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(1.94),
            },
        );

        elements.insert(
            "Sc".to_string(),
            Element {
                name: "Scandium".to_string(),
                atomic_number: 21,
                mass: Mass::new::<dalton>(44.956),
                valence_radius: Length::new::<angstrom>(1.44),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.84),
            },
        );

        elements.insert(
            "Ti".to_string(),
            Element {
                name: "Titanium".to_string(),
                atomic_number: 22,
                mass: Mass::new::<dalton>(47.867),
                valence_radius: Length::new::<angstrom>(1.32),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.76),
            },
        );

        elements.insert(
            "V".to_string(),
            Element {
                name: "Vanadium".to_string(),
                atomic_number: 23,
                mass: Mass::new::<dalton>(50.942),
                valence_radius: Length::new::<angstrom>(1.22),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.71),
            },
        );

        elements.insert(
            "Cr".to_string(),
            Element {
                name: "Chromium".to_string(),
                atomic_number: 24,
                mass: Mass::new::<dalton>(51.996),
                valence_radius: Length::new::<angstrom>(1.18),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.66),
            },
        );

        elements.insert(
            "Mn".to_string(),
            Element {
                name: "Manganese".to_string(),
                atomic_number: 25,
                mass: Mass::new::<dalton>(54.938),
                valence_radius: Length::new::<angstrom>(1.17),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.61),
            },
        );

        elements.insert(
            "Fe".to_string(),
            Element {
                name: "Iron".to_string(),
                atomic_number: 26,
                mass: Mass::new::<dalton>(55.845),
                valence_radius: Length::new::<angstrom>(1.32),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.17),
            },
        );

        elements.insert(
            "Co".to_string(),
            Element {
                name: "Cobalt".to_string(),
                atomic_number: 27,
                mass: Mass::new::<dalton>(58.933),
                valence_radius: Length::new::<angstrom>(1.26),
                valence_electrons: 9,
                atom_radius: Length::new::<angstrom>(1.16),
            },
        );

        elements.insert(
            "Ni".to_string(),
            Element {
                name: "Nickel".to_string(),
                atomic_number: 28,
                mass: Mass::new::<dalton>(58.693),
                valence_radius: Length::new::<angstrom>(1.24),
                valence_electrons: 10,
                atom_radius: Length::new::<angstrom>(1.15),
            },
        );

        elements.insert(
            "Cu".to_string(),
            Element {
                name: "Copper".to_string(),
                atomic_number: 29,
                mass: Mass::new::<dalton>(63.546),
                valence_radius: Length::new::<angstrom>(1.32),
                valence_electrons: 11,
                atom_radius: Length::new::<angstrom>(1.17),
            },
        );

        elements.insert(
            "Zn".to_string(),
            Element {
                name: "Zinc".to_string(),
                atomic_number: 30,
                mass: Mass::new::<dalton>(65.38),
                valence_radius: Length::new::<angstrom>(1.22),
                valence_electrons: 12,
                atom_radius: Length::new::<angstrom>(1.25),
            },
        );

        elements.insert(
            "Ga".to_string(),
            Element {
                name: "Gallium".to_string(),
                atomic_number: 31,
                mass: Mass::new::<dalton>(69.723),
                valence_radius: Length::new::<angstrom>(1.22),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.26),
            },
        );

        elements.insert(
            "Ge".to_string(),
            Element {
                name: "Germanium".to_string(),
                atomic_number: 32,
                mass: Mass::new::<dalton>(72.63),
                valence_radius: Length::new::<angstrom>(1.20),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.22),
            },
        );

        elements.insert(
            "As".to_string(),
            Element {
                name: "Arsenic".to_string(),
                atomic_number: 33,
                mass: Mass::new::<dalton>(74.922),
                valence_radius: Length::new::<angstrom>(1.19),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.19),
            },
        );

        elements.insert(
            "Se".to_string(),
            Element {
                name: "Selenium".to_string(),
                atomic_number: 34,
                mass: Mass::new::<dalton>(78.971),
                valence_radius: Length::new::<angstrom>(1.20),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.16),
            },
        );

        elements.insert(
            "Br".to_string(),
            Element {
                name: "Bromine".to_string(),
                atomic_number: 35,
                mass: Mass::new::<dalton>(79.904),
                valence_radius: Length::new::<angstrom>(1.20),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.14),
            },
        );

        elements.insert(
            "Kr".to_string(),
            Element {
                name: "Krypton".to_string(),
                atomic_number: 36,
                mass: Mass::new::<dalton>(83.798),
                valence_radius: Length::new::<angstrom>(1.16),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.10),
            },
        );

        // Add elements from Rb (37) to Og (118)
        elements.insert(
            "Rb".to_string(),
            Element {
                name: "Rubidium".to_string(),
                atomic_number: 37,
                mass: Mass::new::<dalton>(85.468),
                valence_radius: Length::new::<angstrom>(2.16),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(2.35),
            },
        );

        elements.insert(
            "Sr".to_string(),
            Element {
                name: "Strontium".to_string(),
                atomic_number: 38,
                mass: Mass::new::<dalton>(87.62),
                valence_radius: Length::new::<angstrom>(1.91),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(2.0),
            },
        );

        elements.insert(
            "Y".to_string(),
            Element {
                name: "Yttrium".to_string(),
                atomic_number: 39,
                mass: Mass::new::<dalton>(88.906),
                valence_radius: Length::new::<angstrom>(1.62),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.8),
            },
        );

        elements.insert(
            "Zr".to_string(),
            Element {
                name: "Zirconium".to_string(),
                atomic_number: 40,
                mass: Mass::new::<dalton>(91.224),
                valence_radius: Length::new::<angstrom>(1.45),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.6),
            },
        );

        elements.insert(
            "Nb".to_string(),
            Element {
                name: "Niobium".to_string(),
                atomic_number: 41,
                mass: Mass::new::<dalton>(92.906),
                valence_radius: Length::new::<angstrom>(1.34),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.45),
            },
        );

        elements.insert(
            "Mo".to_string(),
            Element {
                name: "Molybdenum".to_string(),
                atomic_number: 42,
                mass: Mass::new::<dalton>(95.95),
                valence_radius: Length::new::<angstrom>(1.30),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.4),
            },
        );

        elements.insert(
            "Tc".to_string(),
            Element {
                name: "Technetium".to_string(),
                atomic_number: 43,
                mass: Mass::new::<dalton>(98.0),
                valence_radius: Length::new::<angstrom>(1.27),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Ru".to_string(),
            Element {
                name: "Ruthenium".to_string(),
                atomic_number: 44,
                mass: Mass::new::<dalton>(101.07),
                valence_radius: Length::new::<angstrom>(1.25),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.3),
            },
        );

        elements.insert(
            "Rh".to_string(),
            Element {
                name: "Rhodium".to_string(),
                atomic_number: 45,
                mass: Mass::new::<dalton>(102.91),
                valence_radius: Length::new::<angstrom>(1.25),
                valence_electrons: 9,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Pd".to_string(),
            Element {
                name: "Palladium".to_string(),
                atomic_number: 46,
                mass: Mass::new::<dalton>(106.42),
                valence_radius: Length::new::<angstrom>(1.28),
                valence_electrons: 10,
                atom_radius: Length::new::<angstrom>(1.4),
            },
        );

        elements.insert(
            "Ag".to_string(),
            Element {
                name: "Silver".to_string(),
                atomic_number: 47,
                mass: Mass::new::<dalton>(107.87),
                valence_radius: Length::new::<angstrom>(1.34),
                valence_electrons: 11,
                atom_radius: Length::new::<angstrom>(1.6),
            },
        );

        elements.insert(
            "Cd".to_string(),
            Element {
                name: "Cadmium".to_string(),
                atomic_number: 48,
                mass: Mass::new::<dalton>(112.41),
                valence_radius: Length::new::<angstrom>(1.48),
                valence_electrons: 12,
                atom_radius: Length::new::<angstrom>(1.55),
            },
        );

        elements.insert(
            "In".to_string(),
            Element {
                name: "Indium".to_string(),
                atomic_number: 49,
                mass: Mass::new::<dalton>(114.82),
                valence_radius: Length::new::<angstrom>(1.44),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.55),
            },
        );

        elements.insert(
            "Sn".to_string(),
            Element {
                name: "Tin".to_string(),
                atomic_number: 50,
                mass: Mass::new::<dalton>(118.71),
                valence_radius: Length::new::<angstrom>(1.41),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.45),
            },
        );

        elements.insert(
            "Sb".to_string(),
            Element {
                name: "Antimony".to_string(),
                atomic_number: 51,
                mass: Mass::new::<dalton>(121.76),
                valence_radius: Length::new::<angstrom>(1.40),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.45),
            },
        );

        elements.insert(
            "Te".to_string(),
            Element {
                name: "Tellurium".to_string(),
                atomic_number: 52,
                mass: Mass::new::<dalton>(127.60),
                valence_radius: Length::new::<angstrom>(1.36),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.4),
            },
        );

        elements.insert(
            "I".to_string(),
            Element {
                name: "Iodine".to_string(),
                atomic_number: 53,
                mass: Mass::new::<dalton>(126.90),
                valence_radius: Length::new::<angstrom>(1.33),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.4),
            },
        );

        elements.insert(
            "Xe".to_string(),
            Element {
                name: "Xenon".to_string(),
                atomic_number: 54,
                mass: Mass::new::<dalton>(131.29),
                valence_radius: Length::new::<angstrom>(1.31),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.3),
            },
        );

        elements.insert(
            "Cs".to_string(),
            Element {
                name: "Cesium".to_string(),
                atomic_number: 55,
                mass: Mass::new::<dalton>(132.91),
                valence_radius: Length::new::<angstrom>(2.35),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(2.6),
            },
        );

        elements.insert(
            "Ba".to_string(),
            Element {
                name: "Barium".to_string(),
                atomic_number: 56,
                mass: Mass::new::<dalton>(137.33),
                valence_radius: Length::new::<angstrom>(1.98),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(2.15),
            },
        );

        elements.insert(
            "La".to_string(),
            Element {
                name: "Lanthanum".to_string(),
                atomic_number: 57,
                mass: Mass::new::<dalton>(138.91),
                valence_radius: Length::new::<angstrom>(1.69),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.95),
            },
        );

        elements.insert(
            "Ce".to_string(),
            Element {
                name: "Cerium".to_string(),
                atomic_number: 58,
                mass: Mass::new::<dalton>(140.116),
                valence_radius: Length::new::<angstrom>(1.65),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.85),
            },
        );

        elements.insert(
            "Pr".to_string(),
            Element {
                name: "Praseodymium".to_string(),
                atomic_number: 59,
                mass: Mass::new::<dalton>(140.908),
                valence_radius: Length::new::<angstrom>(1.65),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.85),
            },
        );

        elements.insert(
            "Nd".to_string(),
            Element {
                name: "Neodymium".to_string(),
                atomic_number: 60,
                mass: Mass::new::<dalton>(144.242),
                valence_radius: Length::new::<angstrom>(1.64),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.85),
            },
        );

        elements.insert(
            "Pm".to_string(),
            Element {
                name: "Promethium".to_string(),
                atomic_number: 61,
                mass: Mass::new::<dalton>(145.0),
                valence_radius: Length::new::<angstrom>(1.63),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.85),
            },
        );

        elements.insert(
            "Sm".to_string(),
            Element {
                name: "Samarium".to_string(),
                atomic_number: 62,
                mass: Mass::new::<dalton>(150.36),
                valence_radius: Length::new::<angstrom>(1.62),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.85),
            },
        );

        elements.insert(
            "Eu".to_string(),
            Element {
                name: "Europium".to_string(),
                atomic_number: 63,
                mass: Mass::new::<dalton>(151.964),
                valence_radius: Length::new::<angstrom>(1.85),
                valence_electrons: 9,
                atom_radius: Length::new::<angstrom>(1.85),
            },
        );

        elements.insert(
            "Gd".to_string(),
            Element {
                name: "Gadolinium".to_string(),
                atomic_number: 64,
                mass: Mass::new::<dalton>(157.25),
                valence_radius: Length::new::<angstrom>(1.61),
                valence_electrons: 10,
                atom_radius: Length::new::<angstrom>(1.80),
            },
        );

        elements.insert(
            "Tb".to_string(),
            Element {
                name: "Terbium".to_string(),
                atomic_number: 65,
                mass: Mass::new::<dalton>(158.925),
                valence_radius: Length::new::<angstrom>(1.59),
                valence_electrons: 11,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Dy".to_string(),
            Element {
                name: "Dysprosium".to_string(),
                atomic_number: 66,
                mass: Mass::new::<dalton>(162.500),
                valence_radius: Length::new::<angstrom>(1.59),
                valence_electrons: 12,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Ho".to_string(),
            Element {
                name: "Holmium".to_string(),
                atomic_number: 67,
                mass: Mass::new::<dalton>(164.930),
                valence_radius: Length::new::<angstrom>(1.58),
                valence_electrons: 13,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Er".to_string(),
            Element {
                name: "Erbium".to_string(),
                atomic_number: 68,
                mass: Mass::new::<dalton>(167.259),
                valence_radius: Length::new::<angstrom>(1.57),
                valence_electrons: 14,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Tm".to_string(),
            Element {
                name: "Thulium".to_string(),
                atomic_number: 69,
                mass: Mass::new::<dalton>(168.934),
                valence_radius: Length::new::<angstrom>(1.56),
                valence_electrons: 15,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Yb".to_string(),
            Element {
                name: "Ytterbium".to_string(),
                atomic_number: 70,
                mass: Mass::new::<dalton>(173.054),
                valence_radius: Length::new::<angstrom>(1.74),
                valence_electrons: 16,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Lu".to_string(),
            Element {
                name: "Lutetium".to_string(),
                atomic_number: 71,
                mass: Mass::new::<dalton>(174.967),
                valence_radius: Length::new::<angstrom>(1.56),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Hf".to_string(),
            Element {
                name: "Hafnium".to_string(),
                atomic_number: 72,
                mass: Mass::new::<dalton>(178.49),
                valence_radius: Length::new::<angstrom>(1.44),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.55),
            },
        );

        elements.insert(
            "Ta".to_string(),
            Element {
                name: "Tantalum".to_string(),
                atomic_number: 73,
                mass: Mass::new::<dalton>(180.948),
                valence_radius: Length::new::<angstrom>(1.34),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.45),
            },
        );

        elements.insert(
            "W".to_string(),
            Element {
                name: "Tungsten".to_string(),
                atomic_number: 74,
                mass: Mass::new::<dalton>(183.84),
                valence_radius: Length::new::<angstrom>(1.30),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Re".to_string(),
            Element {
                name: "Rhenium".to_string(),
                atomic_number: 75,
                mass: Mass::new::<dalton>(186.207),
                valence_radius: Length::new::<angstrom>(1.28),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Os".to_string(),
            Element {
                name: "Osmium".to_string(),
                atomic_number: 76,
                mass: Mass::new::<dalton>(190.23),
                valence_radius: Length::new::<angstrom>(1.26),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.30),
            },
        );

        elements.insert(
            "Ir".to_string(),
            Element {
                name: "Iridium".to_string(),
                atomic_number: 77,
                mass: Mass::new::<dalton>(192.217),
                valence_radius: Length::new::<angstrom>(1.27),
                valence_electrons: 9,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Pt".to_string(),
            Element {
                name: "Platinum".to_string(),
                atomic_number: 78,
                mass: Mass::new::<dalton>(195.084),
                valence_radius: Length::new::<angstrom>(1.30),
                valence_electrons: 10,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Au".to_string(),
            Element {
                name: "Gold".to_string(),
                atomic_number: 79,
                mass: Mass::new::<dalton>(196.967),
                valence_radius: Length::new::<angstrom>(1.34),
                valence_electrons: 11,
                atom_radius: Length::new::<angstrom>(1.35),
            },
        );

        elements.insert(
            "Hg".to_string(),
            Element {
                name: "Mercury".to_string(),
                atomic_number: 80,
                mass: Mass::new::<dalton>(200.592),
                valence_radius: Length::new::<angstrom>(1.49),
                valence_electrons: 12,
                atom_radius: Length::new::<angstrom>(1.50),
            },
        );

        elements.insert(
            "Tl".to_string(),
            Element {
                name: "Thallium".to_string(),
                atomic_number: 81,
                mass: Mass::new::<dalton>(204.38),
                valence_radius: Length::new::<angstrom>(1.48),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.90),
            },
        );

        elements.insert(
            "Pb".to_string(),
            Element {
                name: "Lead".to_string(),
                atomic_number: 82,
                mass: Mass::new::<dalton>(207.2),
                valence_radius: Length::new::<angstrom>(1.47),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.80),
            },
        );

        elements.insert(
            "Bi".to_string(),
            Element {
                name: "Bismuth".to_string(),
                atomic_number: 83,
                mass: Mass::new::<dalton>(208.980),
                valence_radius: Length::new::<angstrom>(1.46),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.60),
            },
        );

        elements.insert(
            "Po".to_string(),
            Element {
                name: "Polonium".to_string(),
                atomic_number: 84,
                mass: Mass::new::<dalton>(209.0),
                valence_radius: Length::new::<angstrom>(1.46),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.50),
            },
        );

        elements.insert(
            "At".to_string(),
            Element {
                name: "Astatine".to_string(),
                atomic_number: 85,
                mass: Mass::new::<dalton>(210.0),
                valence_radius: Length::new::<angstrom>(1.45),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.50),
            },
        );

        elements.insert(
            "Rn".to_string(),
            Element {
                name: "Radon".to_string(),
                atomic_number: 86,
                mass: Mass::new::<dalton>(222.0),
                valence_radius: Length::new::<angstrom>(1.43),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.50),
            },
        );

        elements.insert(
            "Fr".to_string(),
            Element {
                name: "Francium".to_string(),
                atomic_number: 87,
                mass: Mass::new::<dalton>(223.0),
                valence_radius: Length::new::<angstrom>(2.5),
                valence_electrons: 1,
                atom_radius: Length::new::<angstrom>(2.60),
            },
        );

        elements.insert(
            "Ra".to_string(),
            Element {
                name: "Radium".to_string(),
                atomic_number: 88,
                mass: Mass::new::<dalton>(226.0),
                valence_radius: Length::new::<angstrom>(2.1),
                valence_electrons: 2,
                atom_radius: Length::new::<angstrom>(2.15),
            },
        );

        elements.insert(
            "Ac".to_string(),
            Element {
                name: "Actinium".to_string(),
                atomic_number: 89,
                mass: Mass::new::<dalton>(227.0),
                valence_radius: Length::new::<angstrom>(1.95),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.95),
            },
        );

        elements.insert(
            "Th".to_string(),
            Element {
                name: "Thorium".to_string(),
                atomic_number: 90,
                mass: Mass::new::<dalton>(232.038),
                valence_radius: Length::new::<angstrom>(1.80),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.80),
            },
        );

        elements.insert(
            "Pa".to_string(),
            Element {
                name: "Protactinium".to_string(),
                atomic_number: 91,
                mass: Mass::new::<dalton>(231.036),
                valence_radius: Length::new::<angstrom>(1.80),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.80),
            },
        );

        elements.insert(
            "U".to_string(),
            Element {
                name: "Uranium".to_string(),
                atomic_number: 92,
                mass: Mass::new::<dalton>(238.029),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Np".to_string(),
            Element {
                name: "Neptunium".to_string(),
                atomic_number: 93,
                mass: Mass::new::<dalton>(237.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Pu".to_string(),
            Element {
                name: "Plutonium".to_string(),
                atomic_number: 94,
                mass: Mass::new::<dalton>(244.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Am".to_string(),
            Element {
                name: "Americium".to_string(),
                atomic_number: 95,
                mass: Mass::new::<dalton>(243.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 9,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Cm".to_string(),
            Element {
                name: "Curium".to_string(),
                atomic_number: 96,
                mass: Mass::new::<dalton>(247.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 10,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Bk".to_string(),
            Element {
                name: "Berkelium".to_string(),
                atomic_number: 97,
                mass: Mass::new::<dalton>(247.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 11,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Cf".to_string(),
            Element {
                name: "Californium".to_string(),
                atomic_number: 98,
                mass: Mass::new::<dalton>(251.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 12,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Es".to_string(),
            Element {
                name: "Einsteinium".to_string(),
                atomic_number: 99,
                mass: Mass::new::<dalton>(252.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 13,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Fm".to_string(),
            Element {
                name: "Fermium".to_string(),
                atomic_number: 100,
                mass: Mass::new::<dalton>(257.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 14,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Md".to_string(),
            Element {
                name: "Mendelevium".to_string(),
                atomic_number: 101,
                mass: Mass::new::<dalton>(258.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 15,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "No".to_string(),
            Element {
                name: "Nobelium".to_string(),
                atomic_number: 102,
                mass: Mass::new::<dalton>(259.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 16,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Lr".to_string(),
            Element {
                name: "Lawrencium".to_string(),
                atomic_number: 103,
                mass: Mass::new::<dalton>(262.0),
                valence_radius: Length::new::<angstrom>(1.75),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.75),
            },
        );

        elements.insert(
            "Rf".to_string(),
            Element {
                name: "Rutherfordium".to_string(),
                atomic_number: 104,
                mass: Mass::new::<dalton>(267.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Db".to_string(),
            Element {
                name: "Dubnium".to_string(),
                atomic_number: 105,
                mass: Mass::new::<dalton>(268.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Sg".to_string(),
            Element {
                name: "Seaborgium".to_string(),
                atomic_number: 106,
                mass: Mass::new::<dalton>(269.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Bh".to_string(),
            Element {
                name: "Bohrium".to_string(),
                atomic_number: 107,
                mass: Mass::new::<dalton>(270.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Hs".to_string(),
            Element {
                name: "Hassium".to_string(),
                atomic_number: 108,
                mass: Mass::new::<dalton>(277.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Mt".to_string(),
            Element {
                name: "Meitnerium".to_string(),
                atomic_number: 109,
                mass: Mass::new::<dalton>(278.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 9,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Ds".to_string(),
            Element {
                name: "Darmstadtium".to_string(),
                atomic_number: 110,
                mass: Mass::new::<dalton>(281.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 10,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Rg".to_string(),
            Element {
                name: "Roentgenium".to_string(),
                atomic_number: 111,
                mass: Mass::new::<dalton>(282.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 11,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Cn".to_string(),
            Element {
                name: "Copernicium".to_string(),
                atomic_number: 112,
                mass: Mass::new::<dalton>(285.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 12,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Nh".to_string(),
            Element {
                name: "Nihonium".to_string(),
                atomic_number: 113,
                mass: Mass::new::<dalton>(286.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 3,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Fl".to_string(),
            Element {
                name: "Flerovium".to_string(),
                atomic_number: 114,
                mass: Mass::new::<dalton>(289.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 4,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Mc".to_string(),
            Element {
                name: "Moscovium".to_string(),
                atomic_number: 115,
                mass: Mass::new::<dalton>(290.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 5,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Lv".to_string(),
            Element {
                name: "Livermorium".to_string(),
                atomic_number: 116,
                mass: Mass::new::<dalton>(293.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 6,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Ts".to_string(),
            Element {
                name: "Tennessine".to_string(),
                atomic_number: 117,
                mass: Mass::new::<dalton>(294.0),
                valence_radius: Length::new::<angstrom>(1.7),
                valence_electrons: 7,
                atom_radius: Length::new::<angstrom>(1.7),
            },
        );

        elements.insert(
            "Og".to_string(),
            Element {
                name: "Oganesson".to_string(),
                atomic_number: 118,
                mass: Mass::new::<dalton>(294.0),
                valence_radius: Length::new::<angstrom>(1.6),
                valence_electrons: 8,
                atom_radius: Length::new::<angstrom>(1.6),
            },
        );

        Self { elements }
    }
}

impl PeriodicTable {
    /// Look up an element by its chemical symbol (e.g. `"O"`, `"Fe"`).
    /// Returns a cloned [`Element`] if found, or `None` if the symbol is unknown.
    pub fn get_by_symbol(&self, symbol: &str) -> Option<Element> {
        self.elements.get(symbol).cloned()
    }

    /// Reverse-lookup: find the chemical symbol for a given atomic number.
    ///
    /// Performs a linear scan of the table (O(118)) which is acceptable for
    /// the small number of elements and the typical atom counts in structures.
    /// Returns `None` if no element with that atomic number is registered.
    pub fn symbol_of(&self, atomic_number: u16) -> Option<&str> {
        self.elements
            .iter()
            .find(|(_, el)| el.atomic_number == atomic_number)
            .map(|(sym, _)| sym.as_str())
    }
}
