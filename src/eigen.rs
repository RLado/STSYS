//! Module for calculating eigenvalues of a Sprs matrix
//! 

use rsparse::data::{Sprs, Trpl};

/// Find the real eigenvalues of a `Sprs` matrix using the QR algorithm
/// 
/// max_iter ~ 1_000
/// 
pub fn eig(mat: &Sprs, tol: f64, max_iter: usize) -> (Vec<f64>, Vec<Vec<f64>>) {
    let (mut ak, _) = hessenberg(&mat);
    let eye = eye(mat.m);

    for _ in 0..max_iter {
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

        // check for convergence
        let mut conv = true;
        for i in 0..ak.n-1 {
            for p in ak.p[i]..ak.p[i+1]{ // loop over columns
                if ak.i[p as usize] > i {
                    if ak.x[p as usize].abs() > tol {
                        conv = false;
                        break;
                    }
                }
            }
        }
        if conv {
            break;
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

    // calculate eigenvectors (A - λI) v = 0
    let mut evec = Vec::with_capacity(mat.m);
    let mut b;
    for l in &lambda {
        let a = rsparse::add(&mat, &scxmat(*l, &eye), 1., -1.);
        b = vec![0.; mat.m];
        rsparse::lusol(&a, &mut b, 1, 1e-12);
        evec.push(b);
    }

    
    return (lambda, evec);
}

/// Calculate sparse matrix QR decomposition
/// 
/// # Returns:
/// Q and R matrices, in that order.
/// 
pub fn qr_decomp(s: &Sprs) -> (Sprs, Sprs) {
    let sym = rsparse::sqr(&s, 2, true);
    let qr = rsparse::qr(&s,&sym);
    let eye = eye(s.m);
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