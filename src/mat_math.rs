//! Basic dense matrix and vector operations
//! 

/// Trait defining numbers that can perform mathematical operations
/// 
pub trait Num {}
impl Num for u32 {}
impl Num for u64 {}
impl Num for i32 {}
impl Num for i64 {}
impl Num for f32 {}
impl Num for f64 {}


/// Scalar plus dense matrix
/// 
/// `mat` is overwritten with the answer
/// 
pub fn scpmat<T: Num+std::ops::AddAssign+Copy>(scalar: T, mat: &mut Vec<Vec<T>>) {
    for i in 0..mat.len(){
        for j in 0..mat[0].len(){
            mat[i][j] += scalar;
        }
    }
}

/// Scalar times dense matrix
/// 
/// `mat` is overwritten with the answer
/// 
pub fn scxmat<T: Num+std::ops::MulAssign+Copy>(scalar: T, mat: &mut Vec<Vec<T>>) {
    for i in 0..mat.len(){
        for j in 0..mat[0].len(){
            mat[i][j] *= scalar;
        }
    }
}

/// Dense vector addition
///
pub fn add_vec<T: Num+std::ops::Add<Output = T>+Copy>(a: &Vec<T>, b: &Vec<T>) -> Vec<T> {
    let len = a.len();
    if len != b.len() {
        panic!("Vector sizes do not match.");
    }

    let mut c = Vec::with_capacity(len);
    for i in 0..len {
        c[i] = a[i] + b[i];
    }
    return c;
}

// TODO: Dense vector cross and dot product


//// Matrix add
///
pub fn add_mat <T: Num+std::ops::Add<Output = T>+Copy>(a: &Vec<Vec<T>>, b: &Vec<Vec<T>>) -> Vec<Vec<T>> {
    // Check sizes
    let len_i = a.len();
    let len_j = a[0].len();
    if len_i != b.len() || len_j != b[0].len() {
        panic!("Matrices sizes do not match.");
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


//// Matrix multiplication
/*
    for (int i = 0; i < R1; i++) {
        for (int j = 0; j < C2; j++) {
            rslt[i][j] = 0;
 
            for (int k = 0; k < R2; k++) {
                rslt[i][j] += mat1[i][k] * mat2[k][j];
            }
 
            printf("%d\t", rslt[i][j]);
        }
 
        printf("\n");
    }

*/