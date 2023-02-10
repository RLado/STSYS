//! File writing utilities.
//! 

use std::fs::File;
use std::io::Error;
use std::io::Write;

/// Write the results of a static analysis to a file:
///
/// # Output format:
/// `node#, Fx, Fy, Fz, Mx, My, Mz, Dx, Dy, Dz, Rx, Ry, Rz`
///
pub fn write_results_static(
    out_file: &str,
    f_vec: &Vec<f64>,
    d_vec: &Vec<f64>
) -> Result<(), Error> {
    let mut f = File::create(out_file)?;
    write!(f, "node#, \tFx, \tFy, \tFz, \tMx, \tMy, \tMz, \tDx, \tDy, \tDz, \tRx, \tRy, \tRz\n")?;

    for i in 0..f_vec.len()/6 {
        write!(f, "{}, ", i)?;
        for j in 0..6 {
            write!(f, "{}, ", f_vec[i*6+j])?;
        }
        for j in 0..6 {
            write!(f, "{}, ", d_vec[i*6+j])?;
        }
        write!(f, "\n")?;
    }

    return Ok(());
}

/// Write the results of a static analysis to a file:
///
/// # Output format:
/// `freq, mode`
/// 
pub fn write_results_modal(
    out_file: &str,
    freq: &Vec<f64>,
    modes: &Vec<Vec<f64>>
) -> Result<(), Error> {
    let mut f = File::create(out_file)?;
    write!(f, "freq \t-- \tmode\n")?;

    for i in 0..freq.len() {
        write!(f, "{}\t-- \t", freq[i])?;
        for j in 0..modes[i].len() {
            write!(f, "{}, ", modes[i][j])?;
        }
        write!(f, "\n")?;
    }

    return Ok(());
}
