//! Structural FEM analysis software
//! 

use rsparse;
use rsparse::data::{Sprs, Trpl};

pub mod elements;
pub mod mat_math;
pub mod utils;

pub fn beam12_assemble_elem(global_stff: &Sprs, stff_g: &Sprs, pos: [f64;2]) -> Sprs {
    let mut temp = Trpl{
        m: global_stff.m,
        n: global_stff.n,
        p: Vec::with_capacity(global_stff.nzmax),
        i: Vec::with_capacity(global_stff.nzmax),
        x: Vec::with_capacity(global_stff.nzmax),
    };

    let mut temp_s = Sprs::new();
    temp_s.from_triplet(&temp);
    return rsparse::add(global_stff, &temp_s, 1., 1.);
}