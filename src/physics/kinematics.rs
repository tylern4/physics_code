use crate::physics::constants::*;
use crate::physics::four_momentum::FourMomentum;

/// Calculate Q^2 from beam and scattered electron 4-vectors
/// q^mu^2 = (e_mu - e_mu')^2 = -Q^2
pub fn q2_calc(e_beam: &FourMomentum, e_scattered: &FourMomentum) -> f32 {
    let q = *e_beam - *e_scattered;
    -q.m2()
}

/// Calculate W from beam and scattered electron 4-vectors
/// W^2 = (gamma + P)^2 = M_p^2 - Q^2 + 2*M_p*nu
pub fn w_calc(e_beam: &FourMomentum, e_scattered: &FourMomentum) -> f32 {
    let q = *e_beam - *e_scattered;
    w_calc_from_gamma(&q)
}

/// Calculate W from the photon 4-vector
pub fn w_calc_from_gamma(gamma: &FourMomentum) -> f32 {
    let target = FourMomentum::at_rest(MASS_P);
    (target + *gamma).m()
}

/// Calculate Q^2 from the photon 4-vector
pub fn q2_calc_from_gamma(gamma: &FourMomentum) -> f32 {
    -gamma.m2()
}

/// Calculate Bjorken x
pub fn xb_calc(q2: f32, e_prime: f32, beam_energy: f32) -> f32 {
    let nu = beam_energy - e_prime;
    if nu <= 0.0 {
        return 0.0;
    }
    q2 / (2.0 * MASS_P * nu)
}

/// Calculate Bjorken x from the photon 4-vector
pub fn xb_calc_from_gamma(gamma: &FourMomentum) -> f32 {
    let q2 = q2_calc_from_gamma(gamma);
    let target = FourMomentum::at_rest(MASS_P);
    let dot = gamma.dot4(&target);
    if dot <= 0.0 {
        return 0.0;
    }
    q2 / (2.0 * dot)
}

/// Calculate theta from cz direction cosine (returns degrees)
pub fn theta_calc(cz: f32) -> f32 {
    cz.acos() * RAD2DEG
}

/// Calculate theta in radians from cz direction cosine
pub fn theta_calc_rad(cz: f32) -> f32 {
    cz.acos()
}

/// Calculate theta in radians from cz direction cosine
pub fn theta_rad(cz: f32) -> f32 {
    cz.acos()
}

/// Calculate phi from cx, cy direction cosines (returns degrees)
/// NOTE: matches C++ convention atan2(cosx, cosy) - NOT standard atan2(y,x)
pub fn phi_calc(cx: f32, cy: f32) -> f32 {
    cx.atan2(cy) * RAD2DEG
}

/// Calculate phi in radians from cx, cy direction cosines
pub fn phi_calc_rad(cx: f32, cy: f32) -> f32 {
    cx.atan2(cy)
}

/// Calculate phi in radians from cx, cy direction cosines
pub fn phi_rad(cx: f32, cy: f32) -> f32 {
    cx.atan2(cy)
}

/// Calculate center phi (degrees) - phi shifted by 30 degrees, wrapped to [0, 360]
pub fn center_phi_calc(cx: f32, cy: f32) -> f32 {
    let mut phi0 = cx.atan2(cy) * RAD2DEG + 30.0;
    if phi0 < 0.0 {
        phi0 += 360.0;
    }
    if phi0 > 360.0 {
        phi0 -= 360.0;
    }
    phi0
}

/// Calculate center phi in radians
pub fn center_phi_calc_rad(cx: f32, cy: f32) -> f32 {
    inv_tan(cy, cx)
}

/// Get CLAS sector number (1-6) from phi in degrees [-180, 180]
pub fn get_sector(phi: f32) -> i32 {
    if phi >= 60.0 && phi < 120.0 {
        1
    } else if phi >= 0.0 && phi < 60.0 {
        2
    } else if phi >= -60.0 && phi < 0.0 {
        3
    } else if phi >= -120.0 && phi < -60.0 {
        4
    } else if phi >= -180.0 && phi < -120.0 {
        5
    } else if phi >= 120.0 && phi < 180.0 {
        6
    } else {
        0
    }
}

/// Inverse tangent that returns angle in [0, 2*pi]
pub fn inv_tan(y: f32, x: f32) -> f32 {
    if x > 0.0 && y > 0.0 {
        (y / x).atan() // 1st Quad
    } else if x < 0.0 && y > 0.0 {
        (y / x).atan() + PI // 2nd Quad
    } else if x < 0.0 && y < 0.0 {
        (y / x).atan() + PI // 3rd Quad
    } else if x > 0.0 && y < 0.0 {
        (y / x).atan() + 2.0 * PI // 4th Quad
    } else if x == 0.0 && y > 0.0 {
        PI / 2.0
    } else if x == 0.0 && y < 0.0 {
        3.0 * PI / 2.0
    } else {
        f32::NAN
    }
}

/// Phi of a boosted vector (in [0, 2*pi])
pub fn phi_boosted(vec: &FourMomentum) -> f32 {
    inv_tan(vec.py, vec.px)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sector() {
        assert_eq!(get_sector(90.0), 1);
        assert_eq!(get_sector(30.0), 2);
        assert_eq!(get_sector(-30.0), 3);
        assert_eq!(get_sector(-90.0), 4);
        assert_eq!(get_sector(-150.0), 5);
        assert_eq!(get_sector(150.0), 6);
    }

    #[test]
    fn test_q2_elastic() {
        // For elastic scattering: W = M_p
        let beam = FourMomentum::new(0.0, 0.0, 4.817, MASS_E);
        // Scattered electron at some angle
        let scattered = FourMomentum::new(1.0, 0.0, 3.5, MASS_E);
        let q2 = q2_calc(&beam, &scattered);
        assert!(q2 > 0.0);
    }
}
