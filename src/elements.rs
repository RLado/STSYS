//! This module contains all the elements stiffness matrices and its helper 
//! functions.
//! 
//! Elements are generally named after the number of DoF they have. 

use std::f64::consts::PI;
use crate::data_structs::Coord3D;

use crate::{mat_math::{scxvec, scxmat, add_vec, mul_mat, transpose, cross_prod_3d, dot_prod, add_mat}};

/// Calculate the distance between two `f64` points in 3D space [a -> b]
/// 
/// # Example:
/// ```
/// use stsys_lib::data_structs::Coord3D;
/// use stsys_lib::elements::dist_3df;
///
/// 
/// let a = Coord3D::new(Some(vec![0.,0.,0.]));
/// let b = Coord3D::new(Some(vec![1.,2.,3.]));
/// 
/// assert!(dist_3df(&a, &b) == 14f64.powf(0.5));
/// ```
/// 
pub fn dist_3df(a: &Coord3D<f64>, b: &Coord3D<f64>) -> f64 {
    return f64::sqrt((b.x-a.x).powi(2)+(b.y-a.y).powi(2)+(b.z-a.z).powi(2));
}

/// Vector norm
/// 
/// # Example:
/// ```
/// use stsys_lib::data_structs::Coord3D;
/// use stsys_lib::elements::norm;
/// 
/// let a = vec![1.,2.,3.];
/// 
/// assert!(norm(&a) == 14f64.powf(0.5));
/// ```
/// 
pub fn norm(v: &Vec<f64>) -> f64 {
    let mut n = 0.;
    for i in v {
        n += i.powi(2);
    }
    return n.sqrt();
}

/// Align a vector with a direction
/// # Parameters:
/// - v: vector
/// - d: direction
/// 
pub fn align_vector_3df(v: &Coord3D<f64>, d: &Coord3D<f64>) -> Vec<Vec<f64>> {
    // Normalize v and d
    let v = scxvec(norm(&v.to_vec()).powi(-1), &v.to_vec());
    let d = scxvec(norm(&d.to_vec()).powi(-1), &d.to_vec());

    let w = cross_prod_3d(&v, &d);
    let c = dot_prod(&v, &d);

    let h = 1. / (1. + c);

    let vmat = vec![
        vec![0., -w[2], w[1]],
        vec![w[2], 0., -w[0]],
        vec![-w[1], w[0], 0.],
    ];

    let eye = vec![
        vec![1., 0., 0.],
        vec![0., 1., 0.],
        vec![0., 0., 1.],
    ];

    let vmat_sq = mul_mat(&vmat, &vmat);

    let r = add_mat(&add_mat(&eye, &vmat), &scxmat(h, &vmat_sq));

    return r;
}

/// 3D Timoshenko beam element
/// 3D beam element with 12 DoF and 2 nodes. (reactions X,Y,Z and Mx,My,Mz for 
/// both nodes)
#[derive(Clone, Debug)]
pub struct Beam12{
    /// Nodes
    nodes: [Coord3D<f64>;2], 
    /// Element rotation with respect of local X axis (deg)
    x_rot: f64,
    /// Elastic modulus
    e: Coord3D<f64>,
    /// Shear modulus
    g: Coord3D<f64>,
    /// Moment of inertia {Ix,Iy,Iz}
    i: Coord3D<f64>,
    /// Torsional constant of cross-section \[m^4\]
    j: f64,
    /// Area {Ax,Ay,Az}
    a: Coord3D<f64>,
    /// Stiffness matrix (local axis)
    pub stff_l: Vec<Vec<f64>>,
    /// Stiffness matrix (global axis)
    pub stff_g: Vec<Vec<f64>>,
    /// Rotation matrix
    pub rot: Vec<Vec<f64>>,
}

impl Beam12{
    /// Creates a new empty element (without stiffness matrices)
    /// 
    /// # Parameters:
    /// - nodes: Nodes defining the element's span. eg. {{0,0,0},{1,0,0}}
    /// - x_rot: Element rotation with respect of local X axis (deg)
    /// - e: Elastic modulus {Ex, Ey, Ez}
    /// - g: Shear modulus {Gx, Gy, Gz}
    /// - i: Moment of inertia {Ix,Iy,Iz}
    /// - j: Torsional constant of cross-section
    /// - a: Area {Ax,Ay,Az}
    /// 
    pub fn new(nodes: [Coord3D<f64>;2], x_rot: f64, e: Coord3D<f64>, g: Coord3D<f64>, i: Coord3D<f64>, j: f64, a: Coord3D<f64>)-> Beam12{
        // create an instance of Beam12
        let elem = Beam12{
            nodes: nodes,
            x_rot: x_rot,
            e: e,
            g: g,
            i: i,
            j: j,
            a: a,
            stff_l: Vec::new(),
            stff_g: Vec::new(),
            rot: Vec::new(),
        };

        return elem;
    }

    /// Calculates the element's intrinsic properties
    /// 
    /// ⚠ If you use this function more than once on a `Beam12` be mindful of the
    /// rotation applied to `e`,`g`, and `i`.
    /// 
    /// # Parameters:
    /// - nodes: Nodes defining the element's span. eg. {{0,0,0},{1,0,0}}
    /// - x_rot: Element rotation with respect of local X axis (deg)
    /// - e: Elastic modulus {Ex, Ey, Ez}
    /// - g: Shear modulus {Gx, Gy, Gz}
    /// - i: Moment of inertia {Ix,Iy,Iz}. Note: Ix = Torsional constant of cross-section (j)
    /// - a: Area {Ax,Ay,Az}
    /// 
    pub fn gen_stff(&mut self) {
        // calculate elements properties
        let l = dist_3df(&self.nodes[0], &self.nodes[1]); //elem_len
        // let a = l/2.;

        // populate the rotation matrix
        let p2 = vec![l, 0., 0.]; // p1 = [0., 0., 0.]

        // align p3 to specified beam orientation
        let x_rot_mat = vec![
            vec![1., 0., 0.],
            vec![0., f64::cos(self.x_rot/180.*PI), -f64::sin(self.x_rot/180.*PI)],
            vec![0., f64::sin(self.x_rot/180.*PI), f64::cos(self.x_rot/180.*PI)],
        ];

        // rotate properties to align with beam orientation
        let e_temp = &transpose(&mul_mat(&x_rot_mat, &transpose(&vec![self.e.to_vec()])))[0];
        let g_temp = &transpose(&mul_mat(&x_rot_mat, &transpose(&vec![self.g.to_vec()])))[0];
        let i_temp = &transpose(&mul_mat(&x_rot_mat, &transpose(&vec![self.i.to_vec()])))[0];
        self.e = Coord3D::new(Some(e_temp.to_vec()));
        self.g = Coord3D::new(Some(g_temp.to_vec()));
        self.i = Coord3D::new(Some(i_temp.to_vec()));

        // create a transformation matrix T3 that aligns P with nodes
        let t3 = align_vector_3df(&Coord3D::new(Some(p2)), &Coord3D::new(Some(add_vec(&self.nodes[1].to_vec(), &self.nodes[0].to_vec(), 1., -1.))));

        // assemble the 12x12 rotation matrix
        //  [ T3  0   0   0 ]
        //  [ 0   T3  0   0 ]
        //  [ 0   0   T3  0 ]
        //  [ 0   0   0   T3]
        self.rot = vec![vec![0.;12];12];
        for d in (0..12).step_by(3){
            for i in 0..3{
                for j in 0..3{
                    self.rot[d+i][d+j] = t3[i][j];
                }
            }
        }

        // populate the stiffness matrix
        let phi_y = 12.*(self.e.z*self.i.z/(self.g.y*self.a.y*l.powi(2)));
        let phi_z = 12.*(self.e.y*self.i.y/(self.g.z*self.a.z*l.powi(2)));

        let kz = scxmat(
            self.e.z*self.i.z/((1. + phi_y)*l.powi(2)), 
            &vec![
                vec![12., 6.*l, -12., 6.*l],
                vec![0., (4.+phi_y)*l.powi(2), -6.*l, (2.-phi_y)*l.powi(2)],
                vec![0., 0., 12., -6.*l],
                vec![0., 0., 0., (4.+phi_y)*l.powi(2)]
            ]
        );
        let ky = scxmat(
            self.e.y*self.i.y/((1. + phi_z)*l.powi(2)), 
            &vec![
                vec![12., -6.*l, -12., -6.*l],
                vec![0., (4.+phi_z)*l.powi(2), 6.*l, (2.-phi_z)*l.powi(2)],
                vec![0., 0., 12., 6.*l],
                vec![0., 0., 0., (4.+phi_z)*l.powi(2)]
            ]
        );

        let eal = self.e.x*self.a.x/l;
        let gjl = self.g.x*self.j/l;
        self.stff_l = vec![
            vec![eal, 0., 0., 0., 0., 0., -eal, 0., 0., 0., 0., 0.], // 1
            vec![0., kz[0][0], 0., 0., 0., kz[0][1], 0., kz[0][2], 0., 0., 0., kz[0][3]], // 2
            vec![0., 0., ky[0][0], 0., ky[0][1], 0., 0., 0., ky[0][2], 0., ky[0][3], 0.], // 3
            vec![0., 0., 0., gjl, 0., 0., 0., 0., 0., -gjl, 0., 0.], // 4
            vec![0., 0., ky[0][1], 0., ky[1][1], 0., 0., 0., ky[1][2], 0., ky[1][3], 0.], // 5
            vec![0., kz[0][1], 0., 0., 0., kz[1][1], 0., kz[1][2], 0., 0., 0., kz[1][3]], // 6
            vec![-eal, 0., 0., 0., 0., 0., eal, 0., 0., 0., 0., 0.], // 7
            vec![0., kz[0][2], 0., 0., 0., kz[1][2], 0., kz[2][2], 0., 0., 0., kz[2][3]], // 8
            vec![0., 0., ky[0][2], 0., ky[1][2], 0., 0., 0., ky[2][2], 0., ky[2][3], 0.], // 9
            vec![0., 0., 0., -gjl, 0., 0., 0., 0., 0., gjl, 0., 0.], // 10
            vec![0., 0., ky[0][3], 0., ky[1][3], 0., 0., 0., ky[2][3], 0., ky[3][3], 0.], // 11
            vec![0., kz[0][3], 0., 0., 0., kz[1][3], 0., kz[2][3], 0., 0., 0., kz[3][3]], // 12
        ];

        // transform local stiffness matrix into global stiffness matrix
        self.stff_g = mul_mat(&mul_mat(&transpose(&self.rot), &self.stff_l), &self.rot);
    }
}