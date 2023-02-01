//! This module contains all the elements stiffness matrices and its helper 
//! functions.
//! 
//! Elements are generally named after the number of DoF they have. 

use crate::{utils::Coord3D, mat_math::scxmat};

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

/// 3D beam element with 12 DoF and 2 nodes. (reactions X,Y,Z and Mx,My,Mz for 
/// both nodes)
#[derive(Clone, Debug)]
pub struct Beam12{
    /// Nodes
    pub nodes: [Coord3D<f64>;2], 
    /// Elastic modulus
    pub e: f64,
    /// Shear modulus
    pub g: f64,
    /// Moment of inertia {None,Iy,Iz}
    pub i: Coord3D<f64>,
    /// Torsional constant of cross-section [m^4]
    pub j: f64,
    /// Area {Ax,Ay,Az}
    pub a: Coord3D<f64>,
    /// Stiffness matrix
    pub stff: Vec<Vec<f64>>,
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
        };

        // calculate elements properties
        let elem_len = dist_3df(&nodes[0], &nodes[1]);
        
        // populate the stiffness matrix
        let phi_y = 12.*(elem.e*elem.i.y/(elem.g*elem.a.y*elem_len.powi(2)));
        let phi_z = 12.*(elem.e*elem.i.z/(elem.g*elem.a.z*elem_len.powi(2)));

        let kz = scxmat(
            elem.e*elem.i.z/((1. + phi_y)*elem_len.powi(2)), 
            &vec![
                vec![12., 6.*elem_len, -12., 6.*elem_len],
                vec![0., (4.+phi_y)*elem_len.powi(2), -6.*elem_len, (2.-phi_y)*elem_len.powi(2)],
                vec![0., 0., 12., -6.*elem_len],
                vec![0., 0., 0., (4.+phi_y)*elem_len.powi(2)]
            ]
        );
        let ky = scxmat(
            elem.e*elem.i.y/((1. + phi_z)*elem_len.powi(2)), 
            &vec![
                vec![12., -6.*elem_len, -12., -6.*elem_len],
                vec![0., (4.+phi_z)*elem_len.powi(2), 6.*elem_len, (2.-phi_z)*elem_len.powi(2)],
                vec![0., 0., 12., 6.*elem_len],
                vec![0., 0., 0., (4.+phi_z)*elem_len.powi(2)]
            ]
        );

        let eal = elem.e*elem.a.x/elem_len;
        let gkl = elem.g*elem.j/elem_len;
        elem.stff = vec![
            vec![eal, 0., 0., 0., 0., 0., -eal, 0., 0., 0., 0., 0.], // 1
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 2
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 3
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 4
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 5
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 6
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 7
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 8
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 9
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 10
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 11
            vec![0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.], // 12
        ];


        // return the element
        return elem;
    }
}