//! Configuration parsers to create solvable models
//! 

use cfg_mgr;

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
    let config = cfg_mgr::load(&path).unwrap();
    
    if config["analysis_type"].string.trim() == "static" {
        let (elems, connections) = geometry_parser(&path);
        let (f_vec, d_vec) = constraints_parser(&path);
        return (elems, connections, Some(f_vec), Some(d_vec));
    }
    else if config["analysis_type"].string.trim() == "modal" {
        let (elems, connections) = geometry_parser(&path);
        return (elems, connections, None, None);
    }
    else{
        panic!("Invalid analysis type");
    }

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
    let mut kpt1_c_buff = vec![0.;3];
    let mut kpt2_c_buff = vec![0.;3];
    let mut e_buff = vec![0.;3];
    let mut g_buff = vec![0.;3];
    let mut i_buff = vec![0.;3];
    let mut j_buff = 0.;
    let mut a_buff = vec![0.;3];
    let mut x_rot_buff = 0.;

    for mut e in config["model"].string.trim().split(',') {
        e = e.trim();

        if pc != 2 { // read keypoints
            let kpt = &config[e].numeric;
            // read connection unit
            pos_buff[pc] = kpt[0] as usize;
            // store keypoint coordinates
            if pc == 0 { kpt1_c_buff = kpt[1..=3].to_vec(); }
            else { kpt2_c_buff = kpt[1..=3].to_vec(); }
        }
        else { // read material settings
            let m: Vec<&str> = config[e].string.trim().split(',').collect();
            
            // E
            e_buff = config[m[0].trim()].numeric.clone();
            // G
            g_buff = config[m[1].trim()].numeric.clone();
            // I
            i_buff = config[m[2].trim()].numeric.clone();
            // J
            j_buff = config[m[3].trim()].numeric.clone()[0];
            // A
            a_buff = config[m[4].trim()].numeric.clone();
            // x_Rot
            x_rot_buff = config[m[5].trim()].numeric.clone()[0];
        }

        if pc == 2 {
            // flush buffers
            connections.push(pos_buff);
            elems.push(
                elements::Beam12::new(
                    [utils::Coord3D::new(Some(kpt1_c_buff.clone())), utils::Coord3D::new(Some(kpt2_c_buff.clone()))], // mm
                    x_rot_buff,
                    utils::Coord3D::new(Some(e_buff.clone())), // N/mm²
                    utils::Coord3D::new(Some(g_buff.clone())), // N/mm²
                    utils::Coord3D::new(Some(i_buff.clone())),
                    j_buff,
                    utils::Coord3D::new(Some(a_buff.clone())), //mm²
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
pub fn constraints_parser(path: &str) -> (Vec<Option<f64>>, Vec<Option<f64>>) {
    let config = cfg_mgr::load(&path).unwrap();
    let mut f_vec = Vec::new();
    let mut d_vec = Vec::new();

    // gather f_vec
    for mut i in config["f_vec"].string.trim().split(',') {
        i = i.trim();
        for j in &config[i].numeric {
            if j.is_nan() { f_vec.push(None); }
            else { f_vec.push(Some(*j)); }
        }
    }

    // gather d_vec
    for mut i in config["d_vec"].string.trim().split(',') {
        i = i.trim();
        for j in &config[i].numeric {
            if j.is_nan() { d_vec.push(None); }
            else { d_vec.push(Some(*j)); }
        }
    }

    return (f_vec, d_vec);
}