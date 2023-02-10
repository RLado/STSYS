use stsys_lib::{elements, utils, parser};
mod test_utils;

#[test]
fn beam12_1() {
    // square bar 1m long 50mmx50mm
    let mut elem = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![0.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
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
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem1 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![250.,0.,0.])), utils::Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem2 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![500.,0.,0.])), utils::Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem3 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![750.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
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
    (f_vec_sol, d_vec_sol) = stsys_lib::solve_static(f_vec, &stff, d_vec);

    // Check results (NEED VERIFICATION!)
    dbg!(&f_vec_sol, &d_vec_sol);
    test_utils::assert_eq_f_vec(
        &f_vec_sol, 
        &vec![0.0, -500.0, 0.0, 0.0, 0.0, -499999.99999999406, 0.0, 3.637978807091713e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 9.313225746154785e-10, 0.0, 7.275957614183426e-12, 0.0, 0.0, 0.0, 3.259629011154175e-9, 0.0, 499.9999999999927, 0.0, 0.0, 0.0, 1.3969838619232178e-9],
        1e-6
    );
    test_utils::assert_eq_f_vec(
        &d_vec_sol, 
        &vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0005243161881617512, 0.0, 0.0, 0.0, 4.000002560001587e-6, 0.0, 0.001905775782038127, 0.0, 0.0, 0.0, 6.857147245716998e-6, 0.0, 0.0038586643130575846, 0.0, 0.0, 0.0, 8.571434057146244e-6, 0.0, 0.006097267312648582, 0.0, 0.0, 0.0, 9.142862994289323e-6],
        1e-6
    );
}

#[test]
fn beam12_cantilever_4e_load() {
    // build model
    // square bar 1m long 50mmx50mm
    let elem0 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![0.,0.,0.])), utils::Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem1 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![250.,0.,0.])), utils::Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem2 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![500.,0.,0.])), utils::Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem3 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![750.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
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

    // assemble global stiffness matrix
    let stff = stsys_lib::beam12_gen_stiffness(&vec![elem0, elem1, elem2, elem3], &connections);

    // Solve the model
    let f_vec_sol;
    let d_vec_sol;
    (f_vec_sol, d_vec_sol) = stsys_lib::solve_static(f_vec, &stff, d_vec);

    // load model --------------------------------------------------------------
    let (elements_l, connections_l, f_vec_l, d_vec_l) = parser::cfg_parser("tests/test_models/model_1.cfg");

    // assemble global stiffness matrix
    let stff_l = stsys_lib::beam12_gen_stiffness(&elements_l, &connections_l);

    // solve the model
    let f_vec_sol_l;
    let d_vec_sol_l;
    (f_vec_sol_l, d_vec_sol_l) = stsys_lib::solve_static(f_vec_l.unwrap(), &stff_l, d_vec_l.unwrap());

    // check results -----------------------------------------------------------
    dbg!(&f_vec_sol, &d_vec_sol);
    test_utils::assert_eq_f_vec(
        &f_vec_sol, 
        &f_vec_sol_l,
        1e-6
    );
    test_utils::assert_eq_f_vec(
        &d_vec_sol, 
        &d_vec_sol_l,
        1e-6
    );
}

#[test]
fn beam12_cantilever_4e_modal() {
    // build model
    // square bar 1m long 50mmx50mm
    let elem0 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![0.,0.,0.])), utils::Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem1 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![250.,0.,0.])), utils::Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem2 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![500.,0.,0.])), utils::Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem3 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![750.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let connections = vec![[0,1], [1,2], [2,3], [3,4]];

    // assemble global stiffness matrix
    let stff = stsys_lib::beam12_gen_stiffness(&vec![elem0, elem1, elem2, elem3], &connections);

    // solve the model
    let (mut freq, modes) = stsys_lib::modal_free(&stff);

    // check results (NEED VERIFICATION!)
    dbg!(&freq);
    dbg!(&modes);
    let mut gt = vec![
        -2.65851622368382e-05,
        -2.04337832783118e-05,
        -5.26292332509211e-07,
        6.90674711719440e-06,
        2.30665474477819e-05,
        5.48687222883322e-05,
        802128.623638435,
        1592027.79885425,
        1592027.79889512,
        2902128.62362539,
        5497871.37638022,
        7597871.37638004,
        13054606.2815968,
        13054606.2816378,
        48719051.6758052,
        48719051.6758665,
        106048247.346345,
        383686163.346342,
        726865500.653718,
        1004503416.65370,
        313637831054.392,
        313637831054.393,
        343906807930.745,
        343906807930.745,
        647765746177.055,
        647765746177.055,
        948147235456.339,
        948147235456.340,
        1191980957851.76,
        1191980957851.76
    ];

    freq.sort_by(|a, b| a.partial_cmp(b).unwrap());
    gt.sort_by(|a, b| a.partial_cmp(b).unwrap());
    test_utils::assert_eq_f_vec(
        &freq, 
        &gt,
        1e-1
    );
}

#[test]
fn beam12_cantilever_4e_modal_load() {
    // build model
    // square bar 1m long 50mmx50mm
    let elem0 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![0.,0.,0.])), utils::Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem1 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![250.,0.,0.])), utils::Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem2 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![500.,0.,0.])), utils::Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let elem3 = elements::Beam12::new(
        [utils::Coord3D::new(Some(vec![750.,0.,0.])), utils::Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        utils::Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        utils::Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        utils::Coord3D::new(Some(vec![0., 5.20833e5, 5.20833e5])),
        8.78601e5,
        utils::Coord3D::new(Some(vec![2500., 50.*250., 50.*250.])), //mm²
    );

    let connections = vec![[0,1], [1,2], [2,3], [3,4]];

    // assemble global stiffness matrix
    let stff = stsys_lib::beam12_gen_stiffness(&vec![elem0, elem1, elem2, elem3], &connections);

    // solve the model
    let (freq, modes) = stsys_lib::modal_free(&stff);

    // load model --------------------------------------------------------------
    let (elements_l, connections_l, _, _) = parser::cfg_parser("tests/test_models/model_2.cfg");

    // assemble global stiffness matrix
    let stff_l = stsys_lib::beam12_gen_stiffness(&elements_l, &connections_l);

    // solve the model
    let (freq_l, modes_l) = stsys_lib::modal_free(&stff_l);

    // check results -----------------------------------------------------------
    test_utils::assert_eq_f_vec(
        &freq, 
        &freq_l,
        1e-12
    );
    test_utils::assert_eq_f2d_vec(
        &modes, 
        &modes_l,
        1e-12
    );
}