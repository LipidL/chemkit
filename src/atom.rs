use super::coordinate::CoordinateSystem;
use super::periodic_table::Element;
use nalgebra::RealField;
use uom::Conversion;
use uom::si::{SI, Units};

#[derive(Debug, Clone)]
pub struct Atom<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub element: Element,
    pub position: CoordinateSystem<T>,
}
