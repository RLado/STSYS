//! Module for Sprs matrix operations
//! 

use rsparse::data::{Sprs, Trpl};

use crate::mat_math;

/// Find the real eigenvalues of a square `Sprs` matrix using the QR algorithm
/// 
/// # Parameters:
/// mat: an n x n matrix
/// tol: tolerance for convergence
/// max_iter: maximum number of iterations (~200 is a good number)
/// 
/// # Returns:
/// (convergence, eigenvalues, eigenvectors)
/// 
/// Note that the eigenvectors may not be complete if eigenvalues are repeated.
/// 
pub fn eig(mat: &Sprs, tol: f64, max_iter: usize) -> (bool, Vec<f64>, Vec<Vec<f64>>) {
    let (mut ak, _) = hessenberg(&mat);
    let eye = eye(mat.m);

    let mut cb = false; // convergence break

    for _ in 0..max_iter {
        // s_k is a Rayleigh shift that starts on the last element of the diagonal of A_k and moves upwards
        let mut sk = 0.;
        let sd = sub_diag(&ak);
        for i in (0..sd.len()).rev() {
            if sd[i].abs() > tol {
                sk = diag(&ak)[i+1];
                break;
            }
            else if i == 0 {
                cb = true;
            }
        }
        if cb {
            break;
        }

        let smult = scxmat(sk, &eye);

        // QR decomposition of A_k - s_k * I
        let (q, r, p) = qr_decomp(&rsparse::add(&ak, &smult, 1., -1.));
        let q = rsparse::multiply(&p, &q);
        // Add smult back in
        ak = rsparse::add(&rsparse::multiply(&r,&q), &smult, 1., 1.);

        // round off small values to zero (avoid NaN)
        for i in 0..ak.x.len() {
            if ak.x[i].abs() < f64::EPSILON {
                ak.x[i] = 0.;
            }
        }
    }

    let lambda = diag(&ak);

    // calculate eigenvectors (A - λI) v = 0
    let mut evec = Vec::with_capacity(mat.m);
    for l in &lambda {
        let mut tol = f64::min(f64::EPSILON, tol);
        let mut v = vec![vec![]; mat.m];
        while v[0].len() == 0 && tol < 1e6 { // arbitrary 1e6 limit
            v = null(&rsparse::add(&mat, &scxmat(*l, &eye), 1., -1.), Some(tol));
            tol *= 10.;
        }
        if v[0].len() != 0{
            evec.push(mat_math::transpose(&v)[0].clone());
        }
    }

    return (cb, lambda, evec);
}

/// Calculate sparse matrix QR decomposition
/// 
/// # Returns:
/// Qx, R and P matrices, in that order. P is a permutation matrix used for 
/// obtaining Q (Q = P * Qx).
/// ```math
///     A = Q * R
///     A = P * Qx * R
/// ```
///
/// 
pub fn qr_decomp(s: &Sprs) -> (Sprs, Sprs, Sprs) {
    let sym = rsparse::sqr(&s, -1, true);
    let qr = rsparse::qr(&s,&sym);
    let eye = eye(s.n);
    let mut q = eye.clone();

    // H_n = I - beta_n * v_n * v_n^T (v_n is the n-th column of V)
    // Q = H_1 * H_2 * ... * H_n
    for n in 0..s.n {
        let mut vn = Trpl::new();
        vn.m = qr.l.m;

        for i in qr.l.p[n]..qr.l.p[n+1]{
            vn.append(qr.l.i[i as usize], 0, qr.l.x[i as usize]);
        }
        let vn = vn.to_sprs();

        q = rsparse::multiply(&q, &rsparse::add(&eye, &scxmat(qr.b[n], &rsparse::multiply(&vn, &rsparse::transpose(&vn))), 1., -1.));
    }

    // Calculate row permutation matrix P
    let mut p = Trpl::new();
    let pinv = sym.pinv.unwrap()[0..s.n].to_vec();
    for i in &pinv {
        let i = *i as usize;
        p.append(pinv[i] as usize, i, 1.);
    }
    let p = p.to_sprs();

    return (q, qr.u, rsparse::transpose(&p));
}

/// Scalar times sparse matrix
/// 
pub fn scxmat(s: f64, mat: &Sprs) -> Sprs {
    let mut r = mat.clone();
    for i in 0..mat.x.len() {
        r.x[i] = mat.x[i] * s;
    }
    return r;
}

/// Calculate the Hessemberg decomposition of a matrix A
/// 
/// # Parameters:
/// A: an n x n matrix
/// 
/// # Returns:
/// H: an n x n upper Hessemberg matrix
/// Q: an n x n orthogonal matrix
/// 
pub fn hessenberg(a: &Sprs) -> (Sprs, Sprs) {
    let mut h = a.clone();
    let mut q = eye(a.m);

    for k in 0..a.m-2 {
        let mut r = Trpl::new();
        r.m = a.m-(k+1);
        r.n = 1;
        let mut rc = 0; // row count
        for i in h.p[k]..h.p[k+1] {
            if h.i[i as usize] > k {
                r.append(rc, 0, h.x[i as usize]);
                rc += 1;
            }
        }
        let r = r.to_sprs();
        
        let mut u = Trpl::new();
        u.m = a.m-k;
        u.n = 1;
        let r1;
        if r.i[0] == 0 {
            r1 = r.x[0];
        } else {
            r1 = 0.;
        }
        
        u.append(0, 0, -(sign(r1)*not_zero(r1)+is_zero(r1))*rsparse::multiply(&rsparse::transpose(&r), &r).x[0].sqrt());
        let u = u.to_sprs();

        let mut v = rsparse::add(&r, &u, 1., -1.);
        let tvd = rsparse::multiply(&rsparse::transpose(&v), &v).x[0].sqrt();
        for i in 0..v.x.len() {
            v.x[i] = v.x[i] / tvd;
        }

        // Write W matrix by blocks
        let mut w = Trpl::new();
        w.m = a.m;
        w.n = a.n;
        // First top block
        for i in 0..(k+1) {
            for j in 0..w.n {
                if i == j {
                    w.append(i, j, 1.);
                }
            }
        }
        // Calculate v*v^T
        let vv = rsparse::multiply(&v, &rsparse::transpose(&v));
        let vv_dense = vv.to_dense();

        // Second bottom block
        for i in 0..w.m-(k+1){
            for j in 0..w.n-(k+1) {
                if i == j {
                    w.append(i+(k+1), j+(k+1), 1.-2.*vv_dense[i][j]);
                }
                else {
                    w.append(i+(k+1), j+(k+1), -2.*vv_dense[i][j]);
                }
            }
        }
        let w = w.to_sprs();

        let wt = &rsparse::transpose(&w);
        // Calculate H_n
        h = rsparse::multiply(&rsparse::multiply(&w, &h), &wt);
        // Calculate Q_n
        q = rsparse::multiply(&q, &wt);
    }

    return (h, q);
}

/// Find the null space of a square `Sprs` matrix using QR decomposition
/// 
pub fn null(a: &Sprs, tol: Option<f64>) -> Vec<Vec<f64>> {
    let t;
    if tol.is_none() {
        t = f64::EPSILON;
    }
    else {
        t = tol.unwrap();
    }
    let (q, r, p) = qr_decomp(&rsparse::transpose(&a)); // A^T = QR
    let q = rsparse::multiply(&p, &q);
    let rnk = a.n - diag(&r).iter().filter(|&&x| x < t).count();
    let qd = q.to_dense();
    let mut null = vec![vec![0.; q.n-rnk]; q.m];
    for i in 0..q.m {
        for j in rnk..q.n {
            null[i][j-rnk] = qd[i][j];
        }
    }
    return null;
}

/// Create a identity `Sprs` matrix
/// 
pub fn eye(n: usize) -> Sprs {
    // Construct an identity matrix
    let mut eye = Trpl::new();
    for i in 0..n {
        eye.append(i, i, 1.);
    }
    return eye.to_sprs();
}

/// Copy the diagonal of a matrix to a vector
/// 
pub fn diag(mat: &Sprs) -> Vec<f64> {
    let mut diag = vec![0.; mat.m];
    for i in 0..mat.n{
        for p in mat.p[i]..mat.p[i+1]{
            if mat.i[p as usize] == i {
                diag[i] = mat.x[p as usize];
                break;
            }
        }
    }

    return diag;
}

/// Copy the sub-diagonal of a matrix to a vector
/// 
pub fn sub_diag(mat: &Sprs) -> Vec<f64> {
    let mut sdiag = vec![0.; mat.m-1];
    for i in 0..mat.n{
        for p in mat.p[i]..mat.p[i+1]{
            if mat.i[p as usize] == i + 1 {
                sdiag[i] = mat.x[p as usize];
                break;
            }
        }
    }

    return sdiag;
}

/// Compute the 2-norm of a `Sprs` matrix
/// 
pub fn norm2(v: &Sprs) -> f64 {
    let mut norm = 0.;
    for i in 0..v.x.len() {
        norm += v.x[i]*v.x[i];
    }
    return norm.sqrt();
}

/// Return the sign of a given f64
/// 
fn sign(x: f64) -> f64 {
    if x > 0. {
        return 1.;
    } else if x < 0. {
        return -1.;
    } else {
        return 0.;
    }
}

/// Is not zero f64?
/// 
/// if x is less than epsilon, return 0. else 1.
/// 
fn not_zero(x: f64) -> f64 {
    if x.abs() < f64::EPSILON {
        return 0.;
    } else {
        return 1.;
    }
}

/// Is zero f64?
/// 
/// if x is less than epsilon, return 0. else 1.
/// 
fn is_zero(x: f64) -> f64 {
    if x.abs() < f64::EPSILON {
        return 1.;
    } else {
        return 0.;
    }
}


/// Remove a column from a sparse matrix
/// 
/// # Example:
/// ```
/// use stsys_lib::sprs_ops::sprs_remove_column;
/// use rsparse::data::Sprs;
/// 
/// let a = Sprs{
///     nzmax: 4,
///     m: 3,
///     n: 3,
///     p: vec![0, 1, 2, 4],
///     i: vec![0, 1, 2, 0],
///     x: vec![1., 2., 3., 4.],
/// };
/// 
/// let b = sprs_remove_column(&a, 1);
/// 
/// assert_eq!(b.nzmax, 3);
/// assert_eq!(b.m,3);
/// assert_eq!(b.n, 2);
/// assert_eq!(b.to_dense(), vec![vec![1.0, 4.0], vec![0.0, 0.0], vec![0.0, 3.0]]);
/// ```
/// 
pub fn sprs_remove_column(sprs: &Sprs, col: usize) -> Sprs {
    let mut temp = Trpl{
        m: sprs.m,
        n: sprs.n-1,
        p: Vec::new(),
        i: Vec::new(),
        x: Vec::new(),
    };

    for i in 0..sprs.p.len()-1{
        for j in sprs.p[i] as usize..sprs.p[i+1] as usize {
            if i < col {
                temp.append(sprs.i[j], i, sprs.x[j]);
            } else if i > col {
                temp.append(sprs.i[j], i-1, sprs.x[j]);
            }
        }
    }
    
    return temp.to_sprs();
}

/// Remove a row from a sparse matrix
/// 
/// # Example:
/// 
/// ```
/// use rsparse;
/// use rsparse::data::Sprs;
/// use stsys_lib::sprs_ops::sprs_remove_row;
/// 
/// let a = vec![
///     vec![0.951851, 0.980789, 0.538168, 0.597793, 0.729354],
///     vec![0.427680, 0.511328, 0.794301, 0.969392, 0.702270],
///     vec![0.294124, 0.453990, 0.932289, 0.842932, 0.803577],
///     vec![0.045583, 0.318977, 0.735981, 0.090698, 0.312947],
///     vec![0.285703, 0.371392, 0.758594, 0.961243, 0.282974],
/// ];
/// let a_n1 = vec![ // a with no row 1
///     vec![0.951851, 0.980789, 0.538168, 0.597793, 0.729354],
///     vec![0.294124, 0.453990, 0.932289, 0.842932, 0.803577],
///     vec![0.045583, 0.318977, 0.735981, 0.090698, 0.312947],
///     vec![0.285703, 0.371392, 0.758594, 0.961243, 0.282974],
/// ];
/// 
/// let mut a_s = Sprs::new();
/// a_s.from_vec(&a);
/// 
/// let a_n1_s = sprs_remove_row(&a_s, 1);
/// 
/// assert_eq!(a_n1_s.nzmax, 20);
/// assert_eq!(a_n1_s.m, 4);
/// assert_eq!(a_n1_s.n, 5);
/// assert_eq!(
///     &a_n1_s.to_dense(),
///     &a_n1
/// );
/// ```
/// 
pub fn sprs_remove_row(sprs: &Sprs, row: usize) -> Sprs {
    return rsparse::transpose(&sprs_remove_column(&rsparse::transpose(&sprs), row));
}

/// Invert square matrix
/// 
pub fn inv_sq(mat: &Sprs) -> Sprs {
    let mut temp = Trpl{
        m: mat.m,
        n: mat.n,
        p: Vec::new(),
        i: Vec::new(),
        x: Vec::new(),
    };
    let mut col;

    for i in 0..mat.n {
        // Generate identity column
        col = vec![0.; mat.n];
        col[i] = 1.;

        // Solve for column
        rsparse::lusol(&mat, &mut col, 1, 1e-12);
        
        // Append column to temp
        for j in 0..mat.n {
            temp.append(j, i, col[j]);
        }
    }
    
    return temp.to_sprs();
}

/// Invert positive definite matrix
/// 
pub fn inv_pos_def(mat: &Sprs) -> Sprs {
    let mut temp = Trpl{
        m: mat.m,
        n: mat.n,
        p: Vec::new(),
        i: Vec::new(),
        x: Vec::new(),
    };
    let mut col;

    for i in 0..mat.n {
        // Generate identity column
        col = vec![0.; mat.n];
        col[i] = 1.;

        // Solve for column
        rsparse::cholsol(&mat, &mut col, 0);
        
        // Append column to temp
        for j in 0..mat.n {
            temp.append(j, i, col[j]);
        }
    }
    
    return temp.to_sprs();
}