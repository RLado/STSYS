//! Configuration parsers to create solvable models
//! 

use std::collections::HashMap;
use cfg_mgr;
use cfg_mgr::CfgData;

use crate::elements;
use crate::utils;

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
    let config = cfg_mgr::load(&path);

}

/// Model geometry configuration parser
/// 
/// # Parameters:
/// path: Path to model file
/// 
pub fn geometry_parser(path: &str) -> (Vec<elements::Beam12>, Vec<[usize;2]>) {
    let config = cfg_mgr::load(&path).unwrap();
    let mut elems = Vec::new();
    let mut connections = Vec::new();

    // gather the components defined in the model string
    let mut pc = 0; // position count

    // define memory buffers
    let mut pos_buff = [0, 0];
    let mut kpt_c_buff = [vec![0.;3]; 2];
    let mut e_buff;
    let mut g_buff;
    let mut i_buff;
    let mut j_buff;
    let mut a_buff;
    let mut x_rot_buff;

    for e in config["model"].string.trim().split(',') {
        e = e.trim();

        if pc != 2 { // read keypoints
            let kpt = config[e].numeric;
            // read connection unit
            pos_buff[pc] = kpt[0] as usize;
            // store keypoint coordinates
            kpt_c_buff[pc] = kpt[1..=3].to_vec();
        }
        else { // read material settings
            let m = config[e].string.trim().split(',').collect();
            
            // E
            e_buff = config[m[0].trim()].numeric;
            // G
            g_buff = config[m[1].trim()].numeric;
            // I
            i_buff = config[m[2].trim()].numeric;
            // J
            j_buff = config[m[3].trim()].numeric[0];
            // A
            a_buff = config[m[4].trim()].numeric;
            // x_Rot
            x_rot_buff = config[m[5].trim()].numeric[0];
        }

        if pc == 2 {
            // flush buffers
            connections.push(pos_buff);
            elems.push(
                elements::Beam12::new(
                    [utils::Coord3D::new(Some(kpt_c_buff[0]), utils::Coord3D::new(Some(kpt_c_buff[1])))], // mm
                    x_rot_buff,
                    utils::Coord3D::new(Some(e_buff)), // N/mm²
                    utils::Coord3D::new(Some(g_buff)), // N/mm²
                    utils::Coord3D::new(Some(i_buff)),
                    j_buff,
                    utils::Coord3D::new(Some(a_buff)), //mm²
                )
            );
            // reset pc
            pc = 0;
        }
        else {
            pc += 1;
        }
    }

    return (elems, connections);
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
    let config = cfg_mgr::load(&path);
}