//! Support functions and structs
//! 

use std::collections::HashMap;
use rsparse;
use rsparse::data::{Sprs, Trpl};
use cfg_mgr::CfgData;
use crate::mat_math::Num;

/// 3D coord type
/// 
#[derive(Clone, Copy, Debug)]
pub struct Coord3D<T:Num>{
    /// X coordinate
    pub x: T,
    /// Y coordinate
    pub y: T,
    /// Z coordinate
    pub z: T,
}

impl<T: Num> Coord3D<T>{
    /// Create a new Vec3D
    /// 
    pub fn new(val: Option<Vec<T>>) -> Coord3D<T> {
        let r;

        if val.is_some() && val.as_ref().unwrap().len() == 3 {
            let val = val.unwrap();
            r = Coord3D{x: val[0], y: val[1], z: val[2]};
        }
        else if val.is_some(){
            panic!("Invalid input for Vec3D");
        }
        else {
            r = Coord3D{x: T::zero(), y: T::zero(), z: T::zero()};
        }

        return r;
    }
    
    /// Convert object to Vec
    /// 
    pub fn to_vec(&self) -> Vec<T>{
        return vec![self.x, self.y, self.z];
    }
}

/// Remove a column from a sparse matrix
/// 
/// # Example:
/// ```
/// use stsys_lib::utils::sprs_remove_column;
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
/// use stsys_lib::utils::sprs_remove_row;
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

/// Prints a matrix
/// 
pub fn print_matrix<T: Num>(mat: &Vec<Vec<T>>) {
    for row in mat {
        println!("{:?}", row);
    }
}

/// Print all the keys of the configuration HashMap
///
pub fn print_cfg(cfg: &HashMap<String, CfgData>) {
    println!("Configuration file read:\n");
    println!("-------");

    for key in cfg.keys() {
        print!("{}: ", key);

        // Print all numerical values (if any) for a particular key
        for i in 0..cfg.get(key).unwrap().numeric.len() {
            print!("{}, ", cfg.get(key).unwrap().numeric[i]);
        }

        // Print the string data of a key (if any) (separate using ;)
        println!(";{}", cfg.get(key).unwrap().string);
    }

    println!("-------\n");
}