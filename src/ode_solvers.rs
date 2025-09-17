// use nalgebra::Vector6;
// use crate::trajectory::Trajectory;

// pub fn rk4<F>(dt:f64, t:f64, x:Vector6<f64>, dynamics: &F ) -> Vector6<f64>
// where F: Fn(f64, Vector6<f64>) -> Vector6<f64>,
// {
//     let t_half: f64 = t + dt / 2.0;
//     let t_next: f64 = t + dt;

//     let k1: Vector6<f64> = dt * dynamics(t, x);
//     let k2: Vector6<f64> = dt * dynamics(t_half, x + 0.5* k1);
//     let k3: Vector6<f64> = dt * dynamics(t_half, x + 0.5* k2);
//     let k4: Vector6<f64> = dt * dynamics(t_next, x + k3);

//     let x_delta: Vector6<f64> = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
//     x + x_delta
// }

// pub fn rk4_simulate<F>(t_start:f64, t_end:f64, num_steps:u32 , x0:Vector6<f64>, dynamics: &F, trajectory: &mut Trajectory) -> Vector6<f64>
// where F: Fn(f64, Vector6<f64>) -> Vector6<f64>,
// {
//     let dt: f64 = (t_end - t_start ) / (num_steps as f64);
//     let mut x = x0;

//     for i_step in 0..num_steps {
//         let alpha = (i_step as f64) / (num_steps as f64);
//         let t = t_start + alpha * (t_end - t_start);

//         x = rk4(dt, t, x, dynamics );
//         trajectory.states.push(x);
//     }

// }

