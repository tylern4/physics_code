use crate::physics::constants::*;

/// Virtual photon flux calculation
#[derive(Debug, Clone)]
pub struct PhotonFlux {
    beam_energy: f32,
    target_mass: f32,
    beam_momentum: f32,
    nu: f32,
    scattered_energy: f32,
    scattered_momentum: f32,
    w: f32,
    q2: f32,
    flux: f32,
}

impl PhotonFlux {
    pub fn new(w: f32, q2: f32, beam_energy: f32) -> Self {
        let target_mass = MASS_P;
        let beam_momentum = (beam_energy * beam_energy - MASS_E * MASS_E).sqrt();
        let nu = ((w * w + q2) / target_mass - target_mass) / 2.0;
        let scattered_energy = beam_energy - nu;
        let scattered_momentum = (scattered_energy * scattered_energy - MASS_E * MASS_E).sqrt();

        let mut pf = Self {
            beam_energy,
            target_mass,
            beam_momentum,
            nu,
            scattered_energy,
            scattered_momentum,
            w,
            q2,
            flux: 0.0,
        };

        pf.flux = pf.calculate_flux();
        pf
    }

    fn theta_calc(&self) -> f32 {
        ((self.beam_energy * self.scattered_energy - self.q2 / 2.0 - MASS_E * MASS_E)
            / (self.beam_momentum * self.scattered_momentum))
            .acos()
    }

    fn epsilon_calc(&self) -> f32 {
        let theta = self.theta_calc();
        let tan_half_theta = (theta / 2.0).tan();
        let factor = 1.0 + 2.0 * (1.0 + (self.nu * self.nu) / self.q2) * tan_half_theta * tan_half_theta;
        1.0 / factor
    }

    fn calculate_flux(&self) -> f32 {
        let epsilon = self.epsilon_calc();
        FS_ALPHA / (4.0 * PI * self.q2)
            * self.w
            / (self.beam_energy * self.beam_energy * self.target_mass * self.target_mass)
            * (self.w * self.w - self.target_mass * self.target_mass)
            / (1.0 - epsilon)
    }

    pub fn get_flux(&self) -> f32 {
        self.flux
    }
}
