use nalgebra::Vector2;


pub struct Trajectory {

    // Optional name. Should implement something to auto-assign from state library.
    pub name: Option<String>,

    // Vector of the spacecraft state along the trajectory
    pub states: Vector2<f64>,
}

// Need impl for Trajectory::new()

impl Trajectory {

    pub fn new() -> Self {
        Self {
            name: None,
            states: Vec::new(),
        }
    }
}