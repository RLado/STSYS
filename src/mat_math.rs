//! Basic dense matrix and vector operations
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
            fn zero() -> $t { $v }
            #[inline]
            fn is_zero(&self) -> bool { *self == $v }
        }
    }
}

zero_impl!(usize, 0usize);
zero_impl!(u32,  0u32);
zero_impl!(u64,  0u64);
zero_impl!(isize, 0isize);
zero_impl!(i32, 0i32);
zero_impl!(i64, 0i64);
zero_impl!(f32, 0.0f32);
zero_impl!(f64, 0.0f64);


/// Trait defining numbers that can perform mathematical operations
/// 
pub trait Num: Zero+std::ops::Add+std::ops::AddAssign+std::ops::Mul+std::ops::MulAssign+Copy+std::fmt::Debug{}
impl Num for usize{}
impl Num for u32 {}
impl Num for u64 {}
impl Num for isize {}
impl Num for i32 {}
impl Num for i64 {}
impl Num for f32 {}
impl Num for f64 {}


/// Scalar plus dense matrix
/// 
pub fn scpmat<T: Num>(scalar: T, mat: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    let rows = mat.len();
    let columns = mat[0].len();

    let mut r = vec![vec![T::zero();columns];rows];
    for i in 0..mat.len(){
        for j in 0..mat[0].len(){
            r[i][j] = mat[i][j] + scalar;
        }
    }

    return r;
}

/// Scalar times dense matrix
/// 
pub fn scxmat<T: Num + std::ops::Mul<Output = T>>(scalar: T, mat: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    let rows = mat.len();
    let columns = mat[0].len();

    let mut r = vec![vec![T::zero();columns];rows];
    for i in 0..mat.len(){
        for j in 0..mat[0].len(){
            r[i][j] = mat[i][j] * scalar;
        }
    }

    return r;
}

/// Scalar plus dense vector
/// 
pub fn scpvec<T: Num>(scalar: T, v: &Vec<T>) -> Vec<T>{
    let mut r = vec![T::zero(); v.len()];
    for i in 0..v.len(){
        r[i] = v[i] + scalar;
    }

    return r;
}

/// Scalar times dense vector
/// 
pub fn scxvec<T: Num + std::ops::Mul<Output = T>>(scalar: T, v: &Vec<T>) -> Vec<T> {
    let mut r = vec![T::zero(); v.len()];
    for i in 0..v.len(){
        r[i] = v[i] * scalar;
    }

    return r;
}

/// Dense vector addition
///
pub fn add_vec<T: Num+std::ops::Add<Output = T>>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    let len = a.len();
    if len != b.len() {
        panic!("Vector dimensions do not match.");
    }

    let mut c = Vec::with_capacity(len);
    for i in 0..len {
        c[i] = a[i] + b[i];
    }
    return c;
}

//// Dense matrix add
///
pub fn add_mat <T: Num+std::ops::Add<Output = T>>(a: &Vec<Vec<T>>, b: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    // Check sizes
    let len_i = a.len();
    let len_j = a[0].len();
    if len_i != b.len() || len_j != b[0].len() {
        panic!("Matrix dimensions do not match.");
    }

    // Add
    let mut o = Vec::new();
    for i in 0..len_i{
        let mut row = Vec::new();
        for j in 0..len_j{
            row.push(a[i][j] + b[i][j]);
        }
        o.push(row);
    }

    return o;
}


//// Dense matrix multiplication
/// 
pub fn mul_mat<T: Num + std::ops::AddAssign<<T as std::ops::Mul>::Output>> (a: &Vec<Vec<T>>, b: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    let rows_a = a.len();
    let columns_a = a[0].len();
    let rows_b = b.len();
    let columns_b = b[0].len();

    if columns_a != rows_b{
        panic!("Matrix dimensions do not match.")
    }

    let mut r = vec![vec![T::zero();columns_b];rows_a];
    for i in 0..rows_a{
        for j in 0..columns_b{
            for k in 0..rows_b{
                r[i][j] += a[i][k]*b[k][j];
            }
        }
    }

    return r;
} 

/// Transpose a dense matrix
/// 
pub fn transpose<T: Num> (a: &Vec<Vec<T>>) -> Vec<Vec<T>>{
    let rows = a.len();
    let columns = a[0].len();

    let mut r = vec![vec![T::zero();rows];columns];
    for i in 0..rows{
        for j in 0..columns{
            r[j][i] = a[i][j];
        }
    }

    return r;
}