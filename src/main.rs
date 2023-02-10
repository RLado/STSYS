use stsys_lib;

fn main() {
    // Argparser ---------------------------------------------------------------
    let args: Vec<String> = std::env::args().collect();
    let mut config_path: Option<&str> = None;
    let mut result_path: Option<&str> = None;
    //dbg!(&args);

    // Ensure there is an even number of arguments (flag/value pairs)
    assert!(
        (args.len() - 1) % 2 == 0,
        "Uneven number of arguments. Please check your input."
    );
    // Ensure the required arguments are present
    assert!(
        args.contains(&"-c".to_string()),
        "No model configuration file provided. Please provide a configuration file using the '-c' flag."
    );
    assert!(
        args.contains(&"-r".to_string()),
        "No result file provided. Please provide a result file using the '-r' flag."
    );

    // Borrow values into variables (for readability)
    for i in (1..args.len()).step_by(2) {
        match &args[i] as &str {
            "-c" => config_path = Some(&args[i + 1]),
            "-r" => result_path = Some(&args[i + 1]),
            &_ => panic!("{} flag is not recognized", args[i]),
        };
    }
    // -------------------------------------------------------------------------

    // Parse configuration file
    let (analysis_type, elems, connections, f_vec, d_vec) = stsys_lib::parser::cfg_parser(config_path.unwrap());

    // Calculate requested results and write to the output file
    if analysis_type == "static" {
        // Assemble global stiffness matrix
        let stff = stsys_lib::beam12_gen_stiffness(&elems, &connections);

        // Solve static analysis
        let (f, d) = stsys_lib::solve_static(f_vec.unwrap(), &stff, d_vec.unwrap());

        // Write results to file
        _ = stsys_lib::writer::write_results_static(&result_path.unwrap(), &f, &d);
    }
    else if analysis_type == "modal" {
        // Assemble global stiffness matrix
        let stff = stsys_lib::beam12_gen_stiffness(&elems, &connections);

        // Solve modal analysis
        let (freq, modes) = stsys_lib::modal_free(&stff);

        // Write results to file
        _ = stsys_lib::writer::write_results_modal(&result_path.unwrap(), &freq, &modes);
    }
    else {
        panic!("Invalid analysis type");
    }
    
}
