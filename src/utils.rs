//! Support functions
//! 

use std::collections::HashMap;
use cfg_mgr::CfgData;

use crate::data_structs::Num;

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