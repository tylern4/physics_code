use crate::physics::constants::*;
use crate::physics::four_momentum::FourMomentum;
use crate::reader::event::Event;

/// Momentum correction using polynomial coefficients
pub struct MomCorr;

impl MomCorr {
    /// Theta correction factor for electron
    fn theta_correction_factor(phi_e: f32, theta_e: f32, sec: usize) -> f32 {
        if sec == 0 || sec > 6 {
            return 0.0;
        }
        if sec == 2 || sec == 4 {
            return 0.0;
        }

        let s = sec - 1;
        let phi4 = phi_e.powi(4);
        let phi3 = phi_e.powi(3);
        let phi2 = phi_e.powi(2);
        let theta2 = theta_e.powi(2);

        let a = (MOM_CORR_ELECTRON[s][0][0] * theta2 as f64
            + MOM_CORR_ELECTRON[s][0][1] * theta_e as f64
            + MOM_CORR_ELECTRON[s][0][2]) * phi4 as f64;
        let b = (MOM_CORR_ELECTRON[s][1][0] * theta2 as f64
            + MOM_CORR_ELECTRON[s][1][1] * theta_e as f64
            + MOM_CORR_ELECTRON[s][1][2]) * phi3 as f64;
        let c = (MOM_CORR_ELECTRON[s][2][0] * theta2 as f64
            + MOM_CORR_ELECTRON[s][2][1] * theta_e as f64
            + MOM_CORR_ELECTRON[s][2][2]) * phi2 as f64;
        let d = (MOM_CORR_ELECTRON[s][3][0] * theta2 as f64
            + MOM_CORR_ELECTRON[s][3][1] * theta_e as f64
            + MOM_CORR_ELECTRON[s][3][2]) * phi_e as f64;
        let e = MOM_CORR_ELECTRON[s][4][0] * theta2 as f64
            + MOM_CORR_ELECTRON[s][4][1] * theta_e as f64
            + MOM_CORR_ELECTRON[s][4][2];

        (a + b + c + d + e) as f32
    }

    /// Momentum correction factor for electron
    fn p_correction_factor(phi_e: f32, theta_e: f32, sec: usize) -> f32 {
        if sec == 0 || sec > 6 {
            return 1.0;
        }
        if sec == 2 {
            return 1.0;
        }

        let s = sec - 1;
        let phi4 = phi_e.powi(4);
        let phi3 = phi_e.powi(3);
        let phi2 = phi_e.powi(2);
        let theta2 = theta_e.powi(2);

        let a = (MOM_CORR_ELECTRON_P[s][0][0] * theta2 as f64
            + MOM_CORR_ELECTRON_P[s][0][1] * theta_e as f64
            + MOM_CORR_ELECTRON_P[s][0][2]) * phi4 as f64;
        let b = (MOM_CORR_ELECTRON_P[s][1][0] * theta2 as f64
            + MOM_CORR_ELECTRON_P[s][1][1] * theta_e as f64
            + MOM_CORR_ELECTRON_P[s][1][2]) * phi3 as f64;
        let c = (MOM_CORR_ELECTRON_P[s][2][0] * theta2 as f64
            + MOM_CORR_ELECTRON_P[s][2][1] * theta_e as f64
            + MOM_CORR_ELECTRON_P[s][2][2]) * phi2 as f64;
        let d = (MOM_CORR_ELECTRON_P[s][3][0] * theta2 as f64
            + MOM_CORR_ELECTRON_P[s][3][1] * theta_e as f64
            + MOM_CORR_ELECTRON_P[s][3][2]) * phi_e as f64;
        let e = MOM_CORR_ELECTRON_P[s][4][0] * theta2 as f64
            + MOM_CORR_ELECTRON_P[s][4][1] * theta_e as f64
            + MOM_CORR_ELECTRON_P[s][4][2];

        (a + b + c + d + e) as f32
    }

    /// Corrected electron 4-vector
    pub fn corrected_electron(event: &Event) -> FourMomentum {
        let p = event.p(0);
        let sec = event.dc_sect(0) as usize;

        let temp = FourMomentum::new(event.px(0), event.py(0), event.pz(0), MASS_E);

        let theta = temp.theta();
        let phi = temp.phi();

        let theta_corr = theta + Self::theta_correction_factor(phi, theta, sec);
        let p_corr = p; // * Self::p_correction_factor(phi, theta, sec);

        let px = p_corr * theta_corr.sin() * phi.cos();
        let py = p_corr * theta_corr.sin() * phi.sin();
        let pz = p_corr * theta_corr.cos();

        FourMomentum::new(px, py, pz, MASS_E)
    }

    /// Corrected 4-vector for any particle
    pub fn corrected_vector(px: f32, py: f32, pz: f32, particle_type: i32) -> FourMomentum {
        FourMomentum::new(px, py, pz, get_mass(particle_type))
    }
}
