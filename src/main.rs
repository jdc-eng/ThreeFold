use nalgebra::{DMatrix, RowDVector, Vector6, MatrixXx6, OVector};
// use nalgebra::base::dimension::{Const, Dyn};
// use nalgebra::Matrix1x6;
mod dynamics;
mod trajectory;
mod ode_solvers;



fn main() {
    // let initial_state = RowDVector::from_vec(vec![8.2574814261319895E-1, 2.8339816127588502E-28, 7.9884023525316067E-2, -9.3900450555366247E-16, 1.9356640516861148E-1, -9.8582632717838736E-15]);
    // let mut trajectory = trajectory::Trajectory::new(initial_state);
    // let dynamics = dynamics::Dynamics::new(0.01215058560962404);


    // let t_start = 0.0;
    // let t_end = 2.7764143751594395;
    // let num_steps = 1000;
    // let x0 = trajectory.states[0].clone(); // fix this to convert RowDVector to Vector6

    // ode_solvers::rk4_simulate(t_start, t_end, num_steps, x0, &dynamics, &mut trajectory);

    // for state in trajectory.states.iter() {
    //     println!("{:?}", state);
    // }

    println!("Hello, world!");
    let vec = vec![8.2574814261319895E-1, 2.8339816127588502E-28, 7.9884023525316067E-2, -9.3900450555366247E-16, 1.9356640516861148E-1, -9.8582632717838736E-15];
    let mut state = MatrixXx6::from_vec(vec);
    // let state = Matrix3::new(0, 1, 2, 
    //                     3, 4, 5, 
    //                     6, 7, 8);
    let delta = Vector6::from_vec(vec![1e-6, 0.0, 0.0, 0.0, 0.0, 0.0]);
    
    println!("Initial State: {}", &state);
    
    let new_state = &state.push(2.0);

    // Print out the modified state

    // println!("Delta: {}", &delta);
    println!("Modified State: {}", &new_state);



}
