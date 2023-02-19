use stsys_lib::{elements, utils, parser, data_structs::Coord3D};
mod test_utils;

#[test]
fn beam12_1() {
    // square bar 1m long 50mmx50mm
    let mut elem = elements::Beam12::new(
        [Coord3D::new(Some(vec![0.,0.,0.])), Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );
    elem.gen_stff();

    // 12x12 identity matrix
    let mut eye12 = vec![vec![0.;12];12];
    for i in 0..12 {
        eye12[i][i] = 1.;
    }

    utils::print_matrix(&elem.stff_l);
    println!("---");
    utils::print_matrix(&elem.stff_g);
    test_utils::assert_eq_f2d_vec(&elem.stff_l, &elem.stff_g, 1e-12);
}

#[test]
fn beam12_cantilever_4e() {
    // build model
    // square bar 1m long 50mmx50mm
    let elem0 = elements::Beam12::new(
        [Coord3D::new(Some(vec![0.,0.,0.])), Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem1 = elements::Beam12::new(
        [Coord3D::new(Some(vec![250.,0.,0.])), Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem2 = elements::Beam12::new(
        [Coord3D::new(Some(vec![500.,0.,0.])), Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem3 = elements::Beam12::new(
        [Coord3D::new(Some(vec![750.,0.,0.])), Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
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

    // Check results
    dbg!(&f_vec_sol, &d_vec_sol);
    test_utils::assert_eq_f_vec(
        &f_vec_sol, 
        &vec![0.0, -500.0, 0.0, 0.0, 0.0, -500000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 500.0, 0.0, 0.0, 0.0, 0.0],
        1e-2
    );
    test_utils::assert_eq_f_vec(
        &d_vec_sol, 
        &vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.12693, 0.0, 0.0, 0.0, 0.0, 0.0, 0.46814, 0.0, 0.0, 0.0, 0.0, 0.0, 0.95221, 0.0, 0.0, 0.0, 0.0, 0.0, 1.5077, 0.0, 0.0, 0.0, 0.0],
        1e-1
    );
}

#[test]
fn beam12_cantilever_4e_load() {
    // build model
    // square bar 1m long 50mmx50mm
    let elem0 = elements::Beam12::new(
        [Coord3D::new(Some(vec![0.,0.,0.])), Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem1 = elements::Beam12::new(
        [Coord3D::new(Some(vec![250.,0.,0.])), Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem2 = elements::Beam12::new(
        [Coord3D::new(Some(vec![500.,0.,0.])), Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem3 = elements::Beam12::new(
        [Coord3D::new(Some(vec![750.,0.,0.])), Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
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
    let (_, elements_l, connections_l, f_vec_l, d_vec_l) = parser::cfg_parser("tests/test_models/model_1.cfg");

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
        [Coord3D::new(Some(vec![0.,0.,0.])), Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem1 = elements::Beam12::new(
        [Coord3D::new(Some(vec![250.,0.,0.])), Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem2 = elements::Beam12::new(
        [Coord3D::new(Some(vec![500.,0.,0.])), Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem3 = elements::Beam12::new(
        [Coord3D::new(Some(vec![750.,0.,0.])), Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let connections = vec![[0,1], [1,2], [2,3], [3,4]];

    // assemble global mass matrix
    let mass = stsys_lib::beam12_gen_mass(&vec![elem0.clone(), elem1.clone(), elem2.clone(), elem3.clone()], &connections);

    // assemble global stiffness matrix
    let stff = stsys_lib::beam12_gen_stiffness(&vec![elem0, elem1, elem2, elem3], &connections);

    // solve the model
    let (freq, modes) = stsys_lib::modal_free(&mass, &stff);

    // check results (NEED VERIFICATION!)
    dbg!(&freq);
    dbg!(&modes);
    let gt = vec![
        0.07442440855768964,
        0.06740994408305259,
        0.06221630230794882,
        0.03913285068762656,
        0.03250401024068036,
        0.040480329397483424,
        f64::NAN,
        f64::NAN,
        0.02213857785751406,
        0.016820494225038767,
        0.026100805790126202,
        0.025378660979669985,
        0.02468266132164015,
        0.021978565122426894,
        0.012545874093499155,
        f64::NAN,
        0.01666962310110084,
        0.01479406979964629,
        0.008732585400490088,
        0.004562605265055274,
        0.002068655051093652,
        0.009097382995165675,
        0.0046328019247477556,
        0.0016720750486304153,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,    
    ];

    //freq.sort_by(|a, b| a.partial_cmp(b).unwrap());
    //gt.sort_by(|a, b| a.partial_cmp(b).unwrap());
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
        [Coord3D::new(Some(vec![0.,0.,0.])), Coord3D::new(Some(vec![250.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem1 = elements::Beam12::new(
        [Coord3D::new(Some(vec![250.,0.,0.])), Coord3D::new(Some(vec![500.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem2 = elements::Beam12::new(
        [Coord3D::new(Some(vec![500.,0.,0.])), Coord3D::new(Some(vec![750.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let elem3 = elements::Beam12::new(
        [Coord3D::new(Some(vec![750.,0.,0.])), Coord3D::new(Some(vec![1000.,0.,0.]))], // mm
        0.,
        Coord3D::new(Some(vec![210000., 210000., 210000.])), // N/mm²
        Coord3D::new(Some(vec![79000. ,79000., 79000.])), // N/mm²
        Coord3D::new(Some(vec![8.78601e5, 5.20833e5, 5.20833e5])),
        2500., //mm²
        7850., // kg/m³
    );

    let connections = vec![[0,1], [1,2], [2,3], [3,4]];

    // assemble global mass matrix
    let mass = stsys_lib::beam12_gen_mass(&vec![elem0.clone(), elem1.clone(), elem2.clone(), elem3.clone()], &connections);

    // assemble global stiffness matrix
    let stff = stsys_lib::beam12_gen_stiffness(&vec![elem0, elem1, elem2, elem3], &connections);

    // solve the model
    let (freq, modes) = stsys_lib::modal_free(&mass, &stff);

    // load model --------------------------------------------------------------
    let (_, elements_l, connections_l, _, _) = parser::cfg_parser("tests/test_models/model_2.cfg");

    // assemble global mass matrix
    let mass_l = stsys_lib::beam12_gen_mass(&elements_l, &connections_l);

    // assemble global stiffness matrix
    let stff_l = stsys_lib::beam12_gen_stiffness(&elements_l, &connections_l);

    // solve the model
    let (freq_l, modes_l) = stsys_lib::modal_free(&mass_l, &stff_l);

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