use stsys_lib;
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

#[test]
fn beam12_cantilever_4e() {
    // build model
    // square bar 1m long 50mmx50mm
    let elem0 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![0.,0.,0.])), utils::Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem1 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![250.,0.,0.])), utils::Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem2 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![500.,0.,0.])), utils::Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem3 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![750.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let connections = vec![[0,1], [1,2], [2,3], [3,4]];

    // Set constraints
    let f_vec = vec![
        None, None, None, None, None, None,
        Some(0.), Some(0.), Some(0.), Some(0.), Some(0.), Some(0.),
        Some(0.), Some(0.), Some(0.), Some(0.), Some(0.), Some(0.),
        Some(0.), Some(0.), Some(0.), Some(0.), Some(0.), Some(0.),
        Some(0.), Some(500.), Some(0.), Some(0.), Some(0.), Some(0.),
    ];
    let d_vec = vec![
        Some(0.), Some(0.), Some(0.), Some(0.), Some(0.), Some(0.),
        None, None, None, None, None, None,
        None, None, None, None, None, None,
        None, None, None, None, None, None,
        None, None, None, None, None, None,
    ];

    // Assemble global stiffness matrix
    let stff = stsys_lib::beam12_gen_stiffness(&vec![elem0, elem1, elem2, elem3], &connections);

    // Solve the model
    let f_vec_sol;
    let d_vec_sol;
    (f_vec_sol, d_vec_sol) = stsys_lib::solve(f_vec, &stff, d_vec);

    // Check results
    dbg!(&f_vec_sol, &d_vec_sol);
    assert_eq!(f_vec_sol, vec![0.;30]);
    assert_eq!(d_vec_sol, vec![0.;30]);
}