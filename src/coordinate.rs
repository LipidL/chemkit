use nalgebra::{Point3, RealField, Vector3};
use uom::Conversion;
use uom::si::{
    SI, Units,
    angle::{Angle, radian},
    length::{Length, angstrom},
};

// Design notes
//
// `Length<SI<T>, T>` and `Angle<SI<T>, T>` are both parameterised over the
// same `T`, so the choice of `f32` vs `f64` is preserved exactly as before.
//
// Two extra constraints are required on every impl block that calls
// `.get::<unit>()` or `Quantity::new::<unit>()`:
//
//   • `T: Conversion<T, T = T>`
//     The uom conversion-factor type must match the value type.
//     (True for f32 and f64; not true for integers.)
//
//   • `SI<T>: Units<T>`
//     Asserts that the full SI system is populated for storage type T.
//
//   • `angstrom: Conversion<T, T = T>`  (only where angstrom ops are used)
//   • `radian: Conversion<T, T = T>` (only where radian ops are used)
//
//     Rust cannot "unpack" these from `SI<T>: Units<T>` automatically, so
//     they must be stated explicitly in each impl block that calls the
//     corresponding unit method.  At call sites with concrete T = f32 or
//     T = f64 all bounds are trivially satisfied.
//

/// Cartesian coordinates with compile-time unit checking.
///
/// Fields are stored as [`Length<SI<T>, T>`], so the compiler rejects any
/// accidental unit mismatch (e.g. angstroms where metres are expected) at
/// compile time rather than at runtime.
///
/// # Example
/// ```rust,ignore
/// use uom::si::length::meter;
/// use uom::si::f64::Length;
///
/// // (Cartesian is re-exported from the chemkit crate root once the public
/// //  API surface is finalised.)
/// # use chemkit::coordinate::Cartesian;
/// let c = Cartesian::new(
///     Length::new::<meter>(1.0),
///     Length::new::<meter>(2.0),
///     Length::new::<meter>(3.0),
/// );
/// ```
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

/// Constructor — no unit-specific bounds beyond what the struct itself needs.
impl<T> Cartesian<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub fn new(x: Length<SI<T>, T>, y: Length<SI<T>, T>, z: Length<SI<T>, T>) -> Self {
        Self { x, y, z }
    }
}

/// Conversion methods — require the specific unit bounds used internally.
impl<T> Cartesian<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
    angstrom: Conversion<T, T = T>,
    radian: Conversion<T, T = T>,
{
    /// Construct a [`Cartesian`] value from a raw nalgebra vector whose components
    /// are expressed **in angstroms**.
    pub fn from_vector(vector: Vector3<T>) -> Self {
        Self {
            x: Length::new::<angstrom>(vector.x),
            y: Length::new::<angstrom>(vector.y),
            z: Length::new::<angstrom>(vector.z),
        }
    }

    /// Construct a [`Cartesian`] value from a raw nalgebra point whose components
    /// are expressed **in angstroms**.
    pub fn from_point(point: Point3<T>) -> Self {
        Self {
            x: Length::new::<angstrom>(point.x),
            y: Length::new::<angstrom>(point.y),
            z: Length::new::<angstrom>(point.z),
        }
    }

    /// Return a raw nalgebra vector whose components are the coordinate values
    /// expressed **in angstroms**
    pub fn to_vector(&self) -> Vector3<T> {
        Vector3::new(
            self.x.get::<angstrom>(),
            self.y.get::<angstrom>(),
            self.z.get::<angstrom>(),
        )
    }

    /// Return a raw nalgebra point whose components are expressed **in angstroms**.
    pub fn to_point(&self) -> Point3<T> {
        Point3::new(
            self.x.get::<angstrom>(),
            self.y.get::<angstrom>(),
            self.z.get::<angstrom>(),
        )
    }

    /// Convert to spherical coordinates.
    ///
    /// The returned [`Spherical`] value stores `r` as a [`Length`] and both
    /// angles as [`Angle`] quantities — all at the same precision `T`.
    /// The `r` component is expressed in angstroms.
    pub fn to_spherical(&self) -> Spherical<T> {
        let v = self.to_vector(); // raw Vector3<T> in metres
        Spherical {
            r: Length::new::<angstrom>(v.norm()),
            // atan2(z, x) → azimuthal angle θ in the x-z plane
            theta: Angle::new::<radian>(v.z.atan2(v.x)),
            // asin(y / r) → elevation angle φ above the x-z plane
            phi: Angle::new::<radian>(v.y.asin()),
        }
    }
}

/// Spherical coordinates with compile-time unit checking.
///
/// * `r`     — radial distance, stored as [`Length<SI<T>, T>`]
/// * `theta` — azimuthal angle (in the x-z plane), stored as [`Angle<SI<T>, T>`]
/// * `phi`   — elevation angle (from the x-z plane), stored as [`Angle<SI<T>, T>`]
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

/// Constructor.
impl<T> Spherical<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
{
    pub fn new(r: Length<SI<T>, T>, theta: Angle<SI<T>, T>, phi: Angle<SI<T>, T>) -> Self {
        Self { r, theta, phi }
    }
}

/// Conversion methods.
impl<T> Spherical<T>
where
    T: RealField + Copy + Conversion<T, T = T>,
    SI<T>: Units<T>,
    angstrom: Conversion<T, T = T>,
    radian: Conversion<T, T = T>,
{
    /// Convert to Cartesian coordinates.
    ///
    /// All intermediate arithmetic is performed on raw `T` values (metres /
    /// radians) and the result is re-wrapped in [`Length`] quantities.
    pub fn to_cartesian(&self) -> Cartesian<T> {
        let r = self.r.get::<angstrom>();
        let theta = self.theta.get::<radian>();
        let phi = self.phi.get::<radian>();

        Cartesian {
            x: Length::new::<angstrom>(r * theta.cos() * phi.cos()),
            y: Length::new::<angstrom>(r * phi.sin()),
            z: Length::new::<angstrom>(r * theta.sin() * phi.cos()),
        }
    }

    /// Return a raw nalgebra vector **in metres** by converting to Cartesian
    /// first.
    pub fn to_vector(&self) -> Vector3<T> {
        self.to_cartesian().to_vector()
    }

    /// Return a raw nalgebra point **in metres** by converting to Cartesian
    /// first.
    pub fn to_point(&self) -> Point3<T> {
        self.to_cartesian().to_point()
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
