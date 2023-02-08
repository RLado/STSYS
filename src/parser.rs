//! Configuration parsers to create solvable models
//! 

use std::collections::HashMap;
use cfg_mgr;
use cfg_mgr::CfgData;


/// General configuration parser
/// 
/// # Parameters:
/// path: Path to model file
/// 
/// # Returns:
/// - Elements
/// - Connections
/// - f_vec (if needed)
/// - d_vec (if needed)
/// 
pub fn cfg_parser(path: &str) -> (Vec<elements::Beam12>, Vec<[usize;2]>, Option<Vec<Option<f64>>>, Option<Vec<Option<f64>>>) {

}

/// Model geometry configuration parser
/// 
/// # Parameters:
/// path: Path to model file
/// 
pub fn geometry_parser(path: &str) -> (Vec<elements::Beam12>, Vec<[usize;2]>) {

}

/// f_vec / d_vec configuration parser
/// 
/// # Parameters:
/// path: Path to model file
/// 
/// # Returns:
/// - f_vec
/// - d_vec
/// 
pub fn cfg_parser(path: &str) -> (Vec<Option<f64>>, Vec<Option<f64>>) {

}