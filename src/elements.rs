//! This module contains all the elements stiffness matrices and its helper 
//! functions.
//! 
//! Elements are generally named after the number of DoF they have. 

use super::mat_math;

/// 3D beam element with 12 DoF and 2 nodes. (reactions X,Y,Z and Mx,My,Mz for 
/// both nodes)
#[derive(Clone, Debug)]
pub struct Beam12{
    /// Nodes
    nodes: [[f64;3];2], 
    /// Elastic modulus
    e: f64,
    /// Shear modulus
    g: f64,
    /// Moment of inertia {Ix,Iy,Iz}
    i: [f64;3],
    /// Area {Ax,Ay,Az}
    a: [f64;3],
    /// Stiffness matrix
    stff: Vec<Vec<f64>>,
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
    pub fn new(nodes: [[f64;3];2], e: f64, g: f64, i: [f64;3], a: [f64;3]) -> Beam12 {
        // create an instance of Beam12
        let mut elem = Beam12{
            nodes: nodes,
            e: e,
            g: g,
            i: i,
            a: a,
            stff: vec![vec![0.;12];12],
        };
        
        // populate the stiffness matrix


        // return the element
        return elem;
    }
}