//! Module for calculating eigenvalues of a Sprs matrix. Uses nalgebra but it is
//! not memory or compute efficent. **Needs a rewrite.**
//! 

use nalgebra as na;
use nalgebra::Complex;
use rsparse::data::{Sprs, Trpl};

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

/// Calculate sparse matrix QR decomposition
/// 
/// # Returns:
/// Q and R matrices, in that order.
/// 
pub fn qr_decomp(s: &Sprs) -> (Sprs, Sprs) {
    let sym = rsparse::sqr(&s, 2, true);
    let qr = rsparse::qr(&s,&sym);

    let mut eye = Trpl::new();
    for i in 0..s.m {
        eye.append(i, i, 1.);
    }
    let eye = eye.to_sprs();
    let mut q = eye.clone();
    // H_n = I - beta_n * v_n * v_n^T (v_n is the n-th column of V)
    // Q = H_1 * H_2 * ... * H_n
    for n in 0..qr.b.len() {
        let mut vn = Trpl::new();
        for i in qr.l.p[n]..qr.l.p[n+1]{
            vn.append(qr.l.i[i as usize], 0, qr.l.x[i as usize]);
        }
        let vn = vn.to_sprs();
        q = rsparse::multiply(&q, &rsparse::add(&eye, &scxmat(qr.b[n], &rsparse::multiply(&vn, &rsparse::transpose(&vn))), 1., -1.));
    }

    return (q, qr.u);
}

/// Scalar times sparse matrix
/// 
fn scxmat(s: f64, mat: &Sprs) -> Sprs {
    let mut r = mat.clone();
    for i in 0..mat.x.len() {
        r.x[i] = mat.x[i] * s;
    }
    return r;
}

/// Find the real eigenvalues of a `Sprs` matrix using the QR algorithm
/// 
/// iter ~ 100_000
/// 
/// Make this improvements please:
/// https://www.andreinc.net/2021/01/25/computing-eigenvalues-and-eigenvectors-using-qr-decomposition
/// 
pub fn eigen_qr(mat: &Sprs, iter: usize) -> Vec<f64> {
    let mut ak = mat.clone();
    // Construct an identity matrix
    let mut eye = Trpl::new();
    for i in 0..mat.m {
        eye.append(i, i, 1.);
    }
    let eye = eye.to_sprs();
    let mut qq = eye.clone();

    for _ in 0..iter {
        let (q, r) = qr_decomp(&ak);
        ak = rsparse::multiply(&r,&q);
        qq = rsparse::multiply(&qq, &q);
    }

    let mut lambda = vec![0.; ak.m];
    for i in 0..ak.n{
        for p in ak.p[i]..ak.p[i+1]{
            if ak.i[p as usize] == i {
                lambda[i] = ak.x[p as usize];
                break;
            }
        }
    }
    
    return lambda;
}