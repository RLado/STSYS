use stsys_lib::{elements, utils};
mod test_utils;

#[test]
fn beam12_1() {
    // square bar 1m long 50mmx50mm
    let mut elem = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![0.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        utils::Coord3D::new(Some(vec![2500., 50.*1000., 50.*1000.])), //mm²
    );
    elem.gen_stff();

    // 12x12 identity matrix
    let mut eye12 = vec![vec![0.;12];12];
    for i in 0..12 {
        eye12[i][i] = 1.;
    }

    utils::print_matrix(&elem.rot);
    println!("---");
    utils::print_matrix(&eye12);
    assert_eq!(&elem.rot, &eye12);

    println!("---");

    utils::print_matrix(&elem.stff_l);
    println!("---");
    utils::print_matrix(&elem.stff_g);
    assert_eq!(&elem.stff_l, &elem.stff_g);
}