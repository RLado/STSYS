//! Support functions and structs
//! 

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
}

/// Prints a matrix
/// 
pub fn print_matrix<T: Num>(vec: &Vec<Vec<T>>) {
    for row in vec {
        println!("{:?}", row);
    }
}