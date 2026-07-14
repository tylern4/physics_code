use crate::physics::constants::*;
use crate::reader::event::Event;

/// Delta-t timing for particle identification
#[derive(Debug, Clone)]
pub struct DeltaT {
    pub vertex: f32,
    pub elec_array: Vec<f32>,
    pub proton_array: Vec<f32>,
    pub pion_array: Vec<f32>,
    pub kaon_array: Vec<f32>,
}

impl DeltaT {
    /// Create a new DeltaT from an event
    pub fn new(event: &Event) -> Self {
        let vertex = vertex_time(event.sc_t(0), event.sc_r(0), 1.0);

        let elec_array = delta_t_vec(event, MASS_E);
        let proton_array = delta_t_vec(event, MASS_P);
        let pion_array = delta_t_vec(event, MASS_PIP);
        let kaon_array = delta_t_vec(event, MASS_KP);

        Self {
            vertex,
            elec_array,
            proton_array,
            pion_array,
            kaon_array,
        }
    }

    pub fn get_dt_e(&self, part: usize) -> f32 {
        *self.elec_array.get(part).unwrap_or(&0.0)
    }

    pub fn get_dt_p(&self, part: usize) -> f32 {
        *self.proton_array.get(part).unwrap_or(&0.0)
    }

    pub fn get_dt_pi(&self, part: usize) -> f32 {
        *self.pion_array.get(part).unwrap_or(&0.0)
    }

    pub fn get_dt_k(&self, part: usize) -> f32 {
        *self.kaon_array.get(part).unwrap_or(&0.0)
    }

    pub fn get_vertex(&self) -> f32 {
        self.vertex
    }
}

/// Calculate vertex time
pub fn vertex_time(sc_time: f32, sc_pathlength: f32, relativistic_beta: f32) -> f32 {
    sc_time - sc_pathlength / (relativistic_beta * SOL)
}

/// Calculate delta-t for a given mass
pub fn delta_t(mass: f32, momentum: f32, sc_t: f32, sc_r: f32, vertex: f32) -> f32 {
    let cut_beta = 1.0 / (1.0 + (mass / momentum) * (mass / momentum)).sqrt();
    vertex - vertex_time(sc_t, sc_r, cut_beta)
}

/// Calculate delta-t array for a given mass
fn delta_t_vec(event: &Event, mass: f32) -> Vec<f32> {
    let gpart = event.gpart as usize;
    let mut dt_array = vec![0.0; gpart];

    for i in 0..gpart {
        let sct = event.sc_t(i);
        let scr = event.sc_r(i);
        let mom = event.p(i);
        let vertex = vertex_time(event.sc_t(0), event.sc_r(0), 1.0);

        dt_array[i] = delta_t(mass, mom, sct, scr, vertex);
    }

    dt_array
}
