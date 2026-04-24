use nalgebra::{Point3, RealField, Vector3};
use uom::Conversion;
use uom::si::{
    SI, Units,
    angle::{Angle, radian},
    length::Length,
};

/// Cartesian coordinates with compile-time unit checking.
#[derive(Debug, Clone)]
pub struct Cartesian<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub x: Length<SI<T>, T>,
    pub y: Length<SI<T>, T>,
    pub z: Length<SI<T>, T>,
}

impl<T> Cartesian<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    /// Create a [`Cartesian`] directly from three [`Length`] quantities.
    pub fn new(x: Length<SI<T>, T>, y: Length<SI<T>, T>, z: Length<SI<T>, T>) -> Self {
        Self { x, y, z }
    }

    /// Construct a [`Cartesian`] value from a raw nalgebra vector whose
    /// components are expressed in `LenUnit`.
    pub fn from_vector<LenUnit>(vector: Vector3<T>) -> Self
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
    {
        Self {
            x: Length::new::<LenUnit>(vector.x),
            y: Length::new::<LenUnit>(vector.y),
            z: Length::new::<LenUnit>(vector.z),
        }
    }

    /// Construct a [`Cartesian`] value from a raw nalgebra point whose
    /// components are expressed in `LenUnit`.
    pub fn from_point<LenUnit>(point: Point3<T>) -> Self
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
    {
        Self {
            x: Length::new::<LenUnit>(point.x),
            y: Length::new::<LenUnit>(point.y),
            z: Length::new::<LenUnit>(point.z),
        }
    }

    /// Return a raw nalgebra vector whose components are expressed in
    /// `LenUnit`.
    pub fn to_vector<LenUnit>(&self) -> Vector3<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
    {
        Vector3::new(
            self.x.get::<LenUnit>(),
            self.y.get::<LenUnit>(),
            self.z.get::<LenUnit>(),
        )
    }

    /// Return a raw nalgebra point whose components are expressed in
    /// `LenUnit`.
    pub fn to_point<LenUnit>(&self) -> Point3<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
    {
        Point3::new(
            self.x.get::<LenUnit>(),
            self.y.get::<LenUnit>(),
            self.z.get::<LenUnit>(),
        )
    }

    /// Convert to spherical coordinates.
    ///
    /// The radial component `r` of the returned [`Spherical`] is expressed in
    /// `LenUnit`; both angles are always stored in radians
    pub fn to_spherical<LenUnit>(&self) -> Spherical<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        let v = self.to_vector::<LenUnit>();
        Spherical {
            r: Length::new::<LenUnit>(v.norm()),
            // atan2(z, x) → azimuthal angle θ in the x-z plane
            theta: Angle::new::<radian>(v.z.atan2(v.x)),
            // asin(y / r) → elevation angle φ above the x-z plane
            phi: Angle::new::<radian>(v.y.asin()),
        }
    }
}

/// Spherical coordinates with compile-time unit checking.
///
/// Storing angles as [`Angle`] instead of a raw float prevents accidentally
/// passing a value in degrees where radians are expected.
#[derive(Debug, Clone)]
pub struct Spherical<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub r: Length<SI<T>, T>,
    pub theta: Angle<SI<T>, T>,
    pub phi: Angle<SI<T>, T>,
}

impl<T> Spherical<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    /// Create a [`Spherical`] directly from a [`Length`] and two [`Angle`]
    /// quantities.
    pub fn new(r: Length<SI<T>, T>, theta: Angle<SI<T>, T>, phi: Angle<SI<T>, T>) -> Self {
        Self { r, theta, phi }
    }

    /// Convert to Cartesian coordinates.
    pub fn to_cartesian<LenUnit>(&self) -> Cartesian<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        let r = self.r.get::<LenUnit>();
        let theta = self.theta.get::<radian>();
        let phi = self.phi.get::<radian>();

        Cartesian {
            x: Length::new::<LenUnit>(r * theta.cos() * phi.cos()),
            y: Length::new::<LenUnit>(r * phi.sin()),
            z: Length::new::<LenUnit>(r * theta.sin() * phi.cos()),
        }
    }

    /// Return a raw nalgebra vector by converting to Cartesian first.
    ///
    /// Components are expressed in `LenUnit`.
    pub fn to_vector<LenUnit>(&self) -> Vector3<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        self.to_cartesian::<LenUnit>().to_vector::<LenUnit>()
    }

    /// Return a raw nalgebra point by converting to Cartesian first.
    ///
    /// Components are expressed in `LenUnit`.
    pub fn to_point<LenUnit>(&self) -> Point3<T>
    where
        LenUnit: uom::si::length::Unit + Conversion<T, T = T>,
        radian: Conversion<T, T = T>,
    {
        self.to_cartesian::<LenUnit>().to_point::<LenUnit>()
    }
}

#[derive(Debug, Clone)]
pub enum CoordinateSystem<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    Cartesian(Cartesian<T>),
    Spherical(Spherical<T>),
}
