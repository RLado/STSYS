//! This module contains all the elements stiffness matrices and its helper 
//! functions.
//! 
//! Elements are generally named after the number of DoF they have. 

use crate::{utils::Coord3D, mat_math::{scxvec, scxmat, add_vec}};

/// Calculate the distance between two `f64` points in 3D space [a -> b]
/// 
/// # Example:
/// ```
/// use stsys_lib::utils::Coord3D;
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
    return f64::powf((b.x-a.x).powi(2)+(b.y-a.y).powi(2)+(b.z-a.z).powi(2), 0.5);
}

/// 3D Timoshenko beam element
/// 3D beam element with 12 DoF and 2 nodes. (reactions X,Y,Z and Mx,My,Mz for 
/// both nodes)
#[derive(Clone, Debug)]
pub struct Beam12{
    /// Nodes
    nodes: [Coord3D<f64>;2], 
    /// Elastic modulus
    e: f64,
    /// Shear modulus
    g: f64,
    /// Moment of inertia {None,Iy,Iz}
    i: Coord3D<f64>,
    /// Torsional constant of cross-section [m^4]
    j: f64,
    /// Area {Ax,Ay,Az}
    a: Coord3D<f64>,
    /// Stiffness matrix
    pub stff: Vec<Vec<f64>>,
    /// Rotation matrix
    pub rot: Vec<Vec<f64>>,
}

impl Beam12{
    /// Create a new element
    /// 
    /// # Arguments:
    /// - nodes: Nodes defining the element's span. eg. [[0,0,0],[1,0,0]]
    /// - e: Elastic modulus
    /// - g: Shear modulus
    /// - i: Moment of inertia {Ix,Iy,Iz}
    /// - a: Area {Ax,Ay,Az}
    /// 
    pub fn new(nodes: [Coord3D<f64>;2], e: f64, g: f64, i: Coord3D<f64>, j: f64, a: Coord3D<f64>) -> Beam12 {
        // create an instance of Beam12
        let mut elem = Beam12{
            nodes: nodes,
            e: e,
            g: g,
            i: i,
            j: j,
            a: a,
            stff: Vec::new(),
            rot: Vec::new(),
        };

        // calculate elements properties
        let l = dist_3df(&nodes[0], &nodes[1]); //elem_len
        // let a = l/2.;
        
        // populate the stiffness matrix
        let phi_y = 12.*(elem.e*elem.i.y/(elem.g*elem.a.y*l.powi(2)));
        let phi_z = 12.*(elem.e*elem.i.z/(elem.g*elem.a.z*l.powi(2)));

        let kz = scxmat(
            elem.e*elem.i.z/((1. + phi_y)*l.powi(2)), 
            &vec![
                vec![12., 6.*l, -12., 6.*l],
                vec![0., (4.+phi_y)*l.powi(2), -6.*l, (2.-phi_y)*l.powi(2)],
                vec![0., 0., 12., -6.*l],
                vec![0., 0., 0., (4.+phi_y)*l.powi(2)]
            ]
        );
        let ky = scxmat(
            elem.e*elem.i.y/((1. + phi_z)*l.powi(2)), 
            &vec![
                vec![12., -6.*l, -12., -6.*l],
                vec![0., (4.+phi_z)*l.powi(2), 6.*l, (2.-phi_z)*l.powi(2)],
                vec![0., 0., 12., 6.*l],
                vec![0., 0., 0., (4.+phi_z)*l.powi(2)]
            ]
        );

        let eal = elem.e*elem.a.x/l;
        let gjl = elem.g*elem.j/l;
        elem.stff = vec![
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

        // populate the rotation matrix
        // define 3rd node (must be in the xy plane, but not on the x-axis)
        // -- make a vector that goes from node 1 to node 2 and divide it by 2
        let u = scxvec(0.5, &add_vec(&elem.nodes[1].to_vec(), &elem.nodes[0].to_vec(), 1., -1.));
        // -- add unit y-axis component to the vector
        let v = add_vec(&u, &vec![0.,1.,0.], 1., 1.);
        // -- apply the resulting vector to node 1
        let n3 = Coord3D::new(Some(add_vec(&elem.nodes[0].to_vec(), &v, 1., 1.)));
        let n2 = &elem.nodes[1];
        let n1 = &elem.nodes[0];

        // calculate values for T3
        let x21 = n2.x-n1.x;
        let y21 = n2.y-n1.y;
        let z21 = n2.z-n1.z;
        let x31 = n3.x-n1.x;
        let y31 = n3.y-n1.y;
        let z31 = n3.z-n1.z;

        let a123 = f64::powf((y21*z31-y31*z21).powi(2)+(z21*x31-z31*x21).powi(2)+(x21*y31-x31*y21).powi(2), 0.5);

        let lx = x21/l;
        let mx = y21/l;
        let nx = z21/l;

        let lz = 1./(2.*a123)*(y21*z31-y31*z21);
        let mz = 1./(2.*a123)*(z21*x31-z31*x21);
        let nz = 1./(2.*a123)*(x21*y31-x31*y21);

        let ly = mz*nx-nz*mx;
        let my = nz*lx-lz*nx;
        let ny = lz*mx-mz*lx;

        // the 12x12 rotation matrix is made up of 4 copies of:
        let rot_t3 = vec![
            vec![lx, mx, nx],
            vec![ly, my, ny],
            vec![lz, mz, nz],
        ];

        // assemble the 12x12 rotation matrix


        // return the element
        return elem;
    }
}