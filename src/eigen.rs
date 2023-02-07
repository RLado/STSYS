//! Module for calculating eigenvalues of a Sprs matrix. Uses nalgebra but it is
//! not memory or compute efficent. **Needs a rewrite.**
//! 

use nalgebra;
use rsparse::data::{Sprs, Trpl};

/// Calculate eigenvalues of a symetric matrix
/// 
pub fn eig_sym(mat: &Sprs) -> (Vec<f64>, Vec<Vec<f64>>) {
    // Flatten the matrix to convert into nalgebra format
    let flat = flatten(mat);

    let dm = nalgebra::DMatrix::from_vec(mat.m, mat.n, flat);
    let eig = nalgebra::linalg::SymmetricEigen::new(dm);

    let mut e_val = Vec::with_capacity(mat.m);
    let mut e_vec = Vec::with_capacity(mat.m);
    for e in eig.eigenvalues.iter() {
        e_val.push(*e);
    }
    for e in eig.eigenvectors.row_iter() {
        let mut v = Vec::with_capacity(mat.m);
        for i in 0..mat.m {
            v.push(e[i]);
        }
        e_vec.push(v);
    }

    return (e_val, e_vec);
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
        vn.m = qr.l.m;
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
/// iter ~ 10_000
/// 
pub fn eigen_qr(mat: &Sprs, iter: usize) -> Vec<f64> {
    let mut ak = mat.clone();
    // Construct an identity matrix
    let mut eye = Trpl::new();
    for i in 0..mat.m {
        eye.append(i, i, 1.);
    }
    let eye = eye.to_sprs();

    for _ in 0..iter {
        // s_k is the last element of the diagonal of A_k
        let mut sk = 0.;
        for i in ak.p[ak.n-1]..ak.p[ak.n]{
            if ak.i[i as usize] == ak.n-1 {
                sk = ak.x[i as usize];
                break;
            }
        }
        let smult = scxmat(sk, &eye);

        // QR decomposition of A_k - s_k * I
        let (q, r) = qr_decomp(&rsparse::add(&ak, &smult, 1., -1.));
        // Add smult back in
        ak = rsparse::add(&rsparse::multiply(&r,&q), &smult, 1., 1.);

        // round off small values to zero (avoid NaN)
        for i in 0..ak.x.len() {
            if ak.x[i].abs() < f64::EPSILON {
                ak.x[i] = 0.;
            }
        }
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