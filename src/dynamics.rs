// use nalgebra::{DMatrix, RowDVector};
// use nalgebra::Vector6;



// pub struct Dynamics {
//     // Gravitational parameter of the central body
//     pub mu: f64 ,

// }


// impl Dynamics {
//     pub fn new(mu: f64) -> Self {
//         Self { mu: 0.01215058560962404 }
//     }

//     // Function to compute the derivatives of the state vector
//     pub fn compute_derivatives(&self, _t: f64, state: Vector6<f64>) -> Vector6<f64> {
//         let x = state[0];
//         let y = state[1];
//         let z = state[2];
//         let vx = state[3];
//         let vy = state[4];
//         let vz = state[5];
        
//         let m1 = 1.0 - self.mu; // Mass of primary body
//         let m2 = self.mu;     // Mass of secondary body

//         let rb1 = x + self.mu; // Position relative to primary
//         let rb2 = x - (1.0 - self.mu); // Position relative to secondary

//         let r1 = (rb1 * rb1 + y * y + z * z).sqrt();
//         let r2 = (rb2 * rb2 + y * y + z * z).sqrt();


//         let mut derivatives = Vector6::zeros();

//         // Position derivatives (velocity)
//         derivatives.fixed_rows_mut::<3>(0).copy_from(&state.fixed_rows::<3>(3));

//         // Velocity derivatives (acceleration due to gravity)
//         derivatives[3] = x + 2.0 * vy - m1 * (rb1) / (r1 * r1 * r1) - m2 * (rb2) / (r2 * r2 * r2);
//         derivatives[4] = y - 2.0 * vx - m1 * y / (r1 * r1 * r1) - m2 * y / (r2 * r2 * r2);
//         derivatives[5] = -m1 * z / (r1 * r1 * r1) - m2 * z / (r2 * r2 * r2);

//         derivatives
//     }
// }

