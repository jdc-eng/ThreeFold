// use nalgebra::{DMatrix, RowDVector, Vector6};


// pub struct Trajectory {

//     // Optional name. Should implement something to auto-assign from state library.
//     pub name: std::option::Option::None<String>,

//     // Vector of the spacecraft state along the trajectory. make this type DMatrix
//     pub states: Vec<RowDVector<f64>>,

// }


// impl Trajectory {

//     pub fn new(initial_state: RowDVector<f64>) -> Self {
//         Self {
//             name: Option<String,
//             states: vec![initial_state],
//         }
//     }

//     pub fn append_state(&mut self, new_state: RowDVector<f64>) {
//         self.states.push(new_state);
//     }
// }



// // keeps info about the simulation parameters. User will create this
// pub struct Simulation {
//     pub t_start: f64,
//     pub t_end: f64,
//     pub num_steps: usize,
//     pub trajectory: Trajectory,

//     // should also have fields for integration method, which is of type Fn(f64, Vector6<f64>) -> Vector6<f64>. also should inherit Dynamics struct
//     // or just make it of type Propogator, which is a struct that holds integration method
//     pub integration_method: fn(f64, Vector6<f64>) -> Vector6<f64>,

// }


// impl Simulation {
//     pub fn new(t_start: f64, t_end: f64, num_steps: usize, trajectory: Trajectory) -> Self {
//         Self {
//             t_start,
//             t_end,
//             num_steps,
//             trajectory,
//         }
//     }

//     pub fn default() -> Self {
//         Self {
//             t_start: 0.0,
//             t_end: 1.0,
//             num_steps: 100,
//             trajectory: Trajectory::new(RowDVector::from_vec(vec![0.0; 6])),
//         }
//     }
// }


// // maybe another struct that holds the propogation info; like integration method, step size, etc.

