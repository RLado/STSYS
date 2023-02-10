//! Data structures for STSYS
//! 

/// Defines an additive identity element for `Self`.
///
/// # Deriving
///
/// This trait can be automatically be derived using `#[deriving(Zero)]`
/// attribute. If you choose to use this, make sure that the laws outlined in
/// the documentation for `Zero::zero` still hold.
pub trait Zero: Sized + std::ops::Add<Self, Output = Self> {
    /// Returns the additive identity element of `Self`, `0`.
    ///
    /// # Laws
    ///
    /// ```{.text}
    /// a + 0 = a       ∀ a ∈ Self
    /// 0 + a = a       ∀ a ∈ Self
    /// ```
    ///
    fn zero() -> Self;

    /// Returns `true` if `self` is equal to the additive identity.
    fn is_zero(&self) -> bool;
}

macro_rules! zero_impl {
    ($t:ty, $v:expr) => {
        impl Zero for $t {
            #[inline]
            fn zero() -> $t {
                $v
            }
            #[inline]
            fn is_zero(&self) -> bool {
                *self == $v
            }
        }
    };
}

zero_impl!(usize, 0usize);
zero_impl!(u32, 0u32);
zero_impl!(u64, 0u64);
zero_impl!(isize, 0isize);
zero_impl!(i32, 0i32);
zero_impl!(i64, 0i64);
zero_impl!(f32, 0.0f32);
zero_impl!(f64, 0.0f64);

/// Trait defining numbers that can perform mathematical operations
///
pub trait Num:
    Zero
    + std::ops::Add
    + std::ops::AddAssign
    + std::ops::Sub
    + std::ops::SubAssign
    + std::ops::Mul
    + std::ops::MulAssign
    + Copy
    + std::fmt::Debug
{
}
impl Num for usize {}
impl Num for u32 {}
impl Num for u64 {}
impl Num for isize {}
impl Num for i32 {}
impl Num for i64 {}
impl Num for f32 {}
impl Num for f64 {}

/// 3D coord type
/// 
#[derive(Clone, Copy, Debug)]
pub struct Coord3D<T:Num>{
    /// X coordinate
    pub x: T,
    /// Y coordinate
    pub y: T,
    /// Z coordinate
    pub z: T,
}

impl<T: Num> Coord3D<T>{
    /// Create a new Vec3D
    /// 
    pub fn new(val: Option<Vec<T>>) -> Coord3D<T> {
        let r;

        if val.is_some() && val.as_ref().unwrap().len() == 3 {
            let val = val.unwrap();
            r = Coord3D{x: val[0], y: val[1], z: val[2]};
        }
        else if val.is_some(){
            panic!("Invalid input for Vec3D");
        }
        else {
            r = Coord3D{x: T::zero(), y: T::zero(), z: T::zero()};
        }

        return r;
    }
    
    /// Convert object to Vec
    /// 
    pub fn to_vec(&self) -> Vec<T>{
        return vec![self.x, self.y, self.z];
    }
}