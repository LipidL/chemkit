use nalgebra::RealField;
use uom::Conversion;
use uom::si::{
    SI, Units,
    angle::{Angle, radian},
    length::{Length, angstrom},
};

#[derive(Debug, Clone)]
pub struct CellParameters<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub a: Length<SI<T>, T>,
    pub b: Length<SI<T>, T>,
    pub c: Length<SI<T>, T>,
    pub alpha: Angle<SI<T>, T>,
    pub beta: Angle<SI<T>, T>,
    pub gamma: Angle<SI<T>, T>,
}

impl<T> CellParameters<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub fn new(
        a: Length<SI<T>, T>,
        b: Length<SI<T>, T>,
        c: Length<SI<T>, T>,
        alpha: Angle<SI<T>, T>,
        beta: Angle<SI<T>, T>,
        gamma: Angle<SI<T>, T>,
    ) -> Self {
        Self {
            a,
            b,
            c,
            alpha,
            beta,
            gamma,
        }
    }
}

impl<T> CellParameters<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
    angstrom: Conversion<T, T = T>,
    radian: Conversion<T, T = T>,
{
    /// Return a raw nalgebra vector whose components are the cell parameters in angstroms.
    pub fn to_matrix(&self) -> nalgebra::Matrix3<T> {
        let a = self.a.get::<angstrom>();
        let b = self.b.get::<angstrom>();
        let c = self.c.get::<angstrom>();
        let alpha = self.alpha.get::<radian>();
        let beta = self.beta.get::<radian>();
        let gamma = self.gamma.get::<radian>();
        nalgebra::Matrix3::new(
            a,
            b * gamma.cos(),
            c * beta.cos(),
            T::zero(),
            b * gamma.sin(),
            c * (beta.cos() - alpha.cos()) / gamma.sin(),
            T::zero(),
            T::zero(),
            c * alpha.sin() / gamma.sin(),
        )
    }
}
