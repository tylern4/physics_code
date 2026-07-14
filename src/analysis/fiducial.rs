use crate::physics::constants::*;

/// Electron fiducial cut (phi vs theta, sector-dependent)
pub fn electron_fiducial(theta: f32, phi_cent: f32) -> bool {
    let c = [0.04, 0.03, 15.0];
    let y = c[0] * phi_cent * phi_cent + c[1] * phi_cent + c[2];
    theta >= y
}

/// Hadron fiducial cut using Arjun's parameterization
pub fn hadron_fiducial(
    theta: f32,
    theta_rad: f32,
    phi_c: f32,
    sector: i32,
    p: f32,
) -> bool {
    if sector == 0 {
        return false;
    }

    // Theta min cuts per sector
    if theta_rad < 0.174533 {
        return false;
    }
    if sector == 3 && theta_rad < 0.314159 {
        return false;
    }

    // Top line cut
    if p > 0.25 {
        let x = p - 0.2;
        let s = (x.powf(0.02) * 210.0 - 100.0) * (-0.5 * x).exp() - 10.0;
        if theta >= s {
            return false;
        }
    }

    // Min momentum cut
    if p < 0.18 {
        return false;
    }

    // Sector-specific cuts
    let x = match sector {
        1 => p,
        2 => p + 0.1,
        3 => p - 0.01,
        4 => p - 0.01,
        5 => p - 0.01,
        6 => p + 0.03,
        _ => p,
    };
    let offset = match sector {
        5 => -4.0,
        6 => 15.0,
        _ => 0.0,
    };
    let s = (300.0 * x * x * x - 300.0 * x * x + 497.462 * x + 15.0 + offset) * (-2.15 * x).exp() + 14.3;
    if theta >= s {
        return false;
    }

    // Theta min
    let theta_min_par = [2.5, 4.0, 0.001, 4.8];
    let ans = theta_min_par[0] + (theta_min_par[1] / (p + theta_min_par[2])) + theta_min_par[3] * p;
    if ans <= theta_rad * DEG2RAD {
        return false;
    }

    // Fiducial phi bounds
    let sec_idx = (sector - 1) as usize;
    if sec_idx >= 6 {
        return false;
    }
    let phi_min = -(A0MH[sec_idx] * (1.0 - (-A1MH[sec_idx] * (theta - A2MH[sec_idx])).exp()) - A3MH[sec_idx]);
    let phi_max = A0XH[sec_idx] * (1.0 - (-A1XH[sec_idx] * (theta - A2XH[sec_idx])).exp()) + A3XH[sec_idx];

    phi_c >= phi_min && phi_c <= phi_max
}
