//! Module for calculating eigenvalues of a Sprs matrix. Uses nalgebra but it is
//! not memory or compute efficent. **Needs a rewrite.**
//! 

use nalgebra as na;
use nalgebra::Complex;
use rsparse::data::Sprs;

/// Calculate eigenvalues of a Sprs matrix
/// 
pub fn eigenval(mat: &Sprs) -> Vec<Complex<f64>> {
    // Flatten the matrix to convert into nalgebra format
    let flat = flatten(mat);

    let dm = na::DMatrix::from_vec(mat.m, mat.n, flat);
    let eig = dm.complex_eigenvalues();

    let mut r = Vec::with_capacity(mat.m);
    for e in &eig {
        r.push(*e);
    }

    return r;
}

/// Flatten a matrix
/// 
fn flatten (mat: &Sprs) -> Vec<f64> {
    let mat_dense = mat.to_dense();
    let mut flat = Vec::with_capacity(mat.m*mat.n);
    for i in 0..mat.m { // rows
        for j in 0..mat.n { // columns
            flat.push(mat_dense[i][j]);
        }
    }

    return flat;
}