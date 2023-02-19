//! Structural FEM analysis software
//! 

use rsparse;
use rsparse::data::{Sprs, Trpl};

pub mod data_structs;
pub mod elements;
pub mod mat_math;
pub mod utils;
pub mod sprs_ops;
pub mod parser;
pub mod writer;

/// Assembles an element's stiffness/mass matrix into a given global stiffness/mass 
/// matrix
/// 
/// Note that this function expects a local matrix already transformed into the 
/// global axis coordinate system.
/// 
pub fn beam12_assemble_elem(global_mat: &Sprs, mat_g: &Vec<Vec<f64>>, node_n: [usize;2]) -> Sprs {
    // Check the positions are valid in the global stiffness matrix
    let max_p = std::cmp::max(node_n[0], node_n[1]);
    if (max_p+1)*6 > global_mat.m {
        panic!("Atleast one of the given nodes is outside the bounds of the matrix.");
    }

    // Proceed to assemble one element into the matrix
    let mut temp = Trpl{
        m: global_mat.m,
        n: global_mat.n,
        p: Vec::with_capacity(global_mat.nzmax),
        i: Vec::with_capacity(global_mat.nzmax),
        x: Vec::with_capacity(global_mat.nzmax),
    };

    // Element stiffness matrix quadrant's definition
    //  n#1 n#2
    // --------
    // |Q1||Q2| n#1
    // ---++---
    // |Q3||Q4| n#2
    // --------

    for i in 0..6{
        for j in 0..6{
            // Q1
            temp.append(node_n[0]*6+i, node_n[0]*6+j, mat_g[i][j]);
            // Q2
            temp.append(node_n[1]*6+i, node_n[0]*6+j, mat_g[i+6][j]);
            // Q3
            temp.append(node_n[0]*6+i, node_n[1]*6+j, mat_g[i][j+6]);
            // Q4
            temp.append(node_n[1]*6+i, node_n[1]*6+j, mat_g[i+6][j+6]);
        }
    }

    let mut temp_s = Sprs::new();
    temp_s.from_trpl(&temp);
    return rsparse::add(global_mat, &temp_s, 1., 1.);
}

/// Generates a global stiffness matrix given elements and conectivity
/// 
/// {F} = |K|{d}
/// 
pub fn beam12_gen_stiffness(elements: &Vec<elements::Beam12>, conectivity: &Vec<[usize;2]>) -> Sprs {
    let mut mnl = Vec::with_capacity(conectivity.len());
    for c in conectivity{
        mnl.push(std::cmp::max(c[0], c[1]));
    }
    let max_node_num = mnl.iter().max();
    let mut global_stff = Sprs::zeros((*max_node_num.unwrap()+1)*6, (*max_node_num.unwrap()+1)*6, 0);
    
    for i in 0..conectivity.len(){
        let mut elem = elements[i].clone();
        elem.gen_stff();
        global_stff = beam12_assemble_elem(&global_stff, &elem.stff_g, conectivity[i]);
    }

    return global_stff;
}

/// Generates a global mass matrix given elements and conectivity
/// 
pub fn beam12_gen_mass(elements: &Vec<elements::Beam12>, conectivity: &Vec<[usize;2]>) -> Sprs {
    let mut mnl = Vec::with_capacity(conectivity.len());
    for c in conectivity{
        mnl.push(std::cmp::max(c[0], c[1]));
    }
    let max_node_num = mnl.iter().max();
    let mut global_mass = Sprs::zeros((*max_node_num.unwrap()+1)*6, (*max_node_num.unwrap()+1)*6, 0);
    
    for i in 0..conectivity.len(){
        let mut elem = elements[i].clone();
        elem.gen_stff();
        global_mass = beam12_assemble_elem(&global_mass, &elem.mass_g, conectivity[i]);
    }

    return global_mass;
}

/// Statically solve a structure given constraints and the stiffness matrix (elastic)
/// 
/// {F} = |K|{d}
/// Unknowns must be marked using `None`
/// 
/// # Returns:
/// - Vec<f64>: f_vec without unknowns
/// - Vec<f64>: d_vec without unknowns
/// 
pub fn solve_static (f_vec: Vec<Option<f64>>, global_stff: &Sprs, d_vec: Vec<Option<f64>>) -> (Vec<f64>, Vec<f64>) {
    let len_f_vec = f_vec.len();
    let len_d_vec = d_vec.len();

    // Check inputs
    if !(len_f_vec == len_d_vec && len_f_vec == global_stff.n  && len_f_vec == global_stff.m) {
        panic!("Dimensions of input arguments do not agree.");
    }

    // List unknowns
    let mut unk_i_f = Vec::new();
    let mut kn_f = Vec::new();
    for i in 0..len_f_vec {
        if f_vec[i].is_none() == d_vec[i].is_none() {
            panic!("Can not solve with the given constraints, please check your inputs.");
        }
        if f_vec[i].is_none(){
            unk_i_f.push(i);
        }
        else {
            kn_f.push(f_vec[i].unwrap());
        }
    }

    // Remove all rows and columns corresponding to the unknowns from the stiffness matrix
    let mut red_global_stff = global_stff.clone(); // reduced matrix
    for i in unk_i_f.into_iter().rev() {
        red_global_stff = sprs_ops::sprs_remove_column(&red_global_stff, i);
        red_global_stff = sprs_ops::sprs_remove_row(&red_global_stff, i);
    }

    // Solve for the d_vec unknowns
    rsparse::cholsol(&red_global_stff, &mut kn_f, 0);
    // rsparse::lusol(&red_global_stff, &mut kn_f, 1, 1e-12); // No need to use this, it is slower.
    
    // Populate full d_vec
    let mut d_vec_sol = Vec::new();
    let mut j = 0;
    for i in 0..len_d_vec{
        if d_vec[i].is_none(){
            d_vec_sol.push(kn_f[j]);
            j += 1;
        }
        else{
            d_vec_sol.push(d_vec[i].unwrap());
        }
    }

    // Multiply the full d_vec by global_stff to get the full f_vec
    let f_vec_sol = rsparse::gaxpy(&global_stff, &d_vec_sol, &vec![0.;len_f_vec]);

    // Return the solution
    return (f_vec_sol, d_vec_sol);
}

/// Performs a free body modal analysis (WORK IN PROGRESS)
/// 
/// Returns all the real sqrt(eigenvalues) & eigenvectors of the given mass
/// / stiffness matrix.
/// 
pub fn modal_free(global_mass: &Sprs, global_stff: &Sprs) -> (Vec<f64>, Vec<Vec<f64>>) {
    let m_s_mat = rsparse::multiply(&sprs_ops::inv_sq(&global_mass), &global_stff);
    let (_cb, mut freq, modes) = sprs_ops::eig(&m_s_mat, f64::EPSILON, 1_000);

    // λ = ω²
    for i in 0..freq.len(){
        freq[i] = freq[i].sqrt();
    }

    return (freq, modes);
}