use nalgebra::RealField;
use uom::Conversion;
use uom::si::{
    SI, Units,
    angle::{Angle, radian},
    length::Length,
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

    /// Return a raw nalgebra matrix whose columns are the cell vectors in the given length unit.
    pub fn to_matrix<LenUnit>(&self) -> nalgebra::Matrix3<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        let a = self.a.get::<LenUnit>();
        let b = self.b.get::<LenUnit>();
        let c = self.c.get::<LenUnit>();
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

    /// Construct a [`CellParameters`] from a raw lattice matrix whose columns are the cell vectors in the given length unit.
    pub fn from_matrix<LenUnit>(matrix: nalgebra::Matrix3<T>) -> Self
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        // cell vectors
        let va = matrix.column(0);
        let vb = matrix.column(1);
        let vc = matrix.column(2);

        // vector norms
        let va_norm = va.norm();
        let vb_norm = vb.norm();
        let vc_norm = vc.norm();

        // cell parameters
        let a = Length::new::<LenUnit>(va_norm);
        let b = Length::new::<LenUnit>(vb_norm);
        let c = Length::new::<LenUnit>(vc_norm);

        let alpha_cos = vb.dot(&vc) / (vb_norm * vc_norm);
        let beta_cos = va.dot(&vc) / (va_norm * vc_norm);
        let gamma_cos = va.dot(&vb) / (va_norm * vb_norm);
        let alpha = Angle::new::<radian>(alpha_cos.acos());
        let beta = Angle::new::<radian>(beta_cos.acos());
        let gamma = Angle::new::<radian>(gamma_cos.acos());

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
