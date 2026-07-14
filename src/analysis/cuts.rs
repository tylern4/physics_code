use crate::physics::constants::*;
use crate::physics::kinematics;
use crate::reader::event::Event;

use super::delta_t::DeltaT;

/// Fit function for sampling fraction cuts: par[0] + par[1]*x + par[2]*x^5
fn ec_fit_func(x: f64, par: &[f64; 3]) -> f64 {
    par[0] + par[1] * x + par[2] * x * x * x * x * x * x
}

/// Delta-t polynomial: sum(params[i] * p^i)
fn dt_poly4(params: &[f64; 5], p: f32) -> f32 {
    let p = p as f64;
    (params[0] * p * p * p * p
        + params[1] * p * p * p
        + params[2] * p * p
        + params[3] * p
        + params[4]) as f32
}

/// Log-pol2 function for delta-t cuts
fn log_pol2(params: &[f64; 5], p: f32) -> f32 {
    let p = p as f64;
    (params[0] * (params[1] * p).ln() + params[2] * p * p + params[3] * p + params[4]) as f32
}

/// Cherenkov fiducial cut
fn fid_chern(x: f32, y: f32) -> bool {
    let p0 = 48.0;
    let p1 = 1.75;
    x > p0 - p1 * y && x > p0 + p1 * y
}

/// Cherenkov fiducial cut for MC
fn fid_chern_mc(x: f32, y: f32) -> bool {
    let p0 = 48.0;
    let p1 = 1.9;
    x > p0 - p1 * y && x > p0 + p1 * y
}

/// Pion sector cuts
fn pip_sec_cut(sector: i32, p: f32, theta: f32) -> bool {
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
        1 => 0.0,
        2 => 0.0,
        3 => 0.0,
        4 => 0.0,
        5 => -4.0,
        6 => 15.0,
        _ => 0.0,
    };
    let s = (300.0 * x * x * x - 300.0 * x * x + 497.462 * x + 15.0 + offset) * (-2.15 * x).exp() + 14.3;
    theta < s
}

/// Theta minimum cut
fn theta_min(p: f32, theta_rad: f32) -> bool {
    let par = [2.5, 4.0, 0.001, 4.8];
    let ans = par[0] + (par[1] / (p + par[2])) + par[3] * p;
    ans > theta_rad * DEG2RAD
}

/// Cuts struct for particle identification
#[derive(Debug, Clone)]
pub struct Cuts<'a> {
    event: &'a Event,
    pub dt: DeltaT,
    theta: f32,
    phi: f32,
    phi_cent: f32,
    sec: i32,
}

impl<'a> Cuts<'a> {
    pub fn new(event: &'a Event) -> Self {
        let dt = DeltaT::new(event);
        let theta = kinematics::theta_calc(event.cz(0));
        let phi = kinematics::phi_calc(event.cx(0), event.cy(0));
        let sec = event.dc_sect(0);

        let phi_cent = match sec {
            1 => phi - 90.0,
            2 => phi - 30.0,
            3 => phi + 30.0,
            4 => phi + 90.0,
            5 => phi + 150.0,
            6 => phi - 150.0,
            _ => phi,
        };

        Self {
            event,
            dt,
            theta,
            phi,
            phi_cent,
            sec,
        }
    }

    /// Check basic bank requirements
    pub fn check_banks(&self) -> bool {
        let gpart = self.event.gpart;
        if gpart <= 0 || gpart >= 5 {
            return false;
        }
        self.event.q(0) == NEGATIVE
            && self.event.ec(0) > 0
            && self.event.cc(0) > 0
            && self.event.stat(0) > 0
            && self.event.sc(0) > 0
            && self.event.dc(0) > 0
            && self.event.dc_stat(0) > 0
    }

    /// Base electron identification
    pub fn is_electron(&self) -> bool {
        let mut elec = true;
        elec &= self.check_banks();
        elec &= self.event.nphe(0) > 3;
        elec &= self.event.ec_ei(0) >= 0.05;
        elec &= self.event.p(0) > MIN_P_CUT;
        elec
    }

    /// Cherenkov fiducial cut
    pub fn fid_chern_cut(&self) -> bool {
        let a = -0.000785;
        let b = 0.0;
        let c = -0.00168;
        let d = 1.0;

        let p0_x = self.event.dc_xsc(0);
        let p0_y = self.event.dc_ysc(0);
        let p0_z = self.event.dc_zsc(0);
        let n_x = self.event.dc_cxsc(0);
        let n_y = self.event.dc_cysc(0);
        let n_z = self.event.dc_czsc(0);

        let numer = a * p0_x + b * p0_y + c * p0_z + d;
        let denom = a * n_x + b * n_y + c * n_z;

        let t = (numer / denom).abs();
        let t_x = n_x * t;
        let t_y = n_y * t;
        let t_z = n_z * t;

        let final_x = p0_x + t_x;
        let final_y = p0_y + t_y;
        let final_z = p0_z + t_z;

        let mag = (final_x * final_x + final_y * final_y + final_z * final_z).sqrt();
        let cc_theta = (final_z / mag).acos();
        let cc_phi = final_y.atan2(final_x);

        let cc_x = self.event.cc_r(0) * cc_theta.sin() * cc_phi.cos();
        let cc_y = self.event.cc_r(0) * cc_theta.sin() * cc_phi.sin();

        fid_chern(cc_x, cc_y)
    }

    /// Hadron fiducial phi
    pub fn hadron_fid_phi(&self, part: usize) -> f32 {
        let phi = kinematics::phi_calc(self.event.cx(part), self.event.cy(part));
        let sec = self.event.dc_sect(part);
        match sec {
            1 => phi - 90.0,
            2 => phi - 30.0,
            3 => phi + 30.0,
            4 => phi + 90.0,
            5 => phi + 150.0,
            6 => phi - 150.0,
            _ => f32::NAN,
        }
    }

    /// Hadron fiducial phi min
    pub fn hadron_fid_phi_min(&self, theta: f32, sector: usize) -> f32 {
        if sector >= 6 {
            return f32::NAN;
        }
        -(A0MH[sector] * (1.0 - (-A1MH[sector] * (theta - A2MH[sector])).exp()) - A3MH[sector])
    }

    /// Hadron fiducial phi max
    pub fn hadron_fid_phi_max(&self, theta: f32, sector: usize) -> f32 {
        if sector >= 6 {
            return f32::NAN;
        }
        A0XH[sector] * (1.0 - (-A1XH[sector] * (theta - A2XH[sector])).exp()) + A3XH[sector]
    }

    /// Hadron fiducial cut (Arjun's version)
    pub fn hadron_fid_arjun(&self, part: usize) -> bool {
        let theta = kinematics::theta_calc(self.event.cz(part));
        let theta_rad = kinematics::theta_calc_rad(self.event.cz(part));
        let pip_p = self.event.p(part);
        let phi_c = self.hadron_fid_phi(part);
        let sector = self.event.dc_sect(part);

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
        if pip_p > 0.25 {
            let x = pip_p - 0.2;
            let s = (x.powf(0.02) * 210.0 - 100.0) * (-0.5 * x).exp() - 10.0;
            if theta >= s {
                return false;
            }
        }

        // Min momentum cut
        if pip_p < 0.18 {
            return false;
        }

        // Sector-specific cuts
        if !pip_sec_cut(sector, pip_p, theta) {
            return false;
        }

        // Theta min
        if !theta_min(pip_p, theta_rad) {
            return false;
        }

        // Fiducial phi bounds
        let phi_min = self.hadron_fid_phi_min(theta, (sector - 1) as usize);
        let phi_max = self.hadron_fid_phi_max(theta, (sector - 1) as usize);

        phi_c >= phi_min && phi_c <= phi_max
    }

    /// Pi+ identification
    pub fn pip(&self, part: usize) -> bool {
        self.event.q(part) == POSITIVE
            && self.hadron_fid_arjun(part)
            && self.dt_pip_cut(part)
    }

    /// Pi- identification
    pub fn pim(&self, part: usize) -> bool {
        self.event.q(part) == NEGATIVE
            && self.hadron_fid_arjun(part)
            && self.dt_pip_cut(part)
    }

    /// Proton identification
    pub fn prot(&self, part: usize) -> bool {
        self.event.q(part) == POSITIVE
            && self.hadron_fid_arjun(part)
            && self.dt_p_cut(part)
    }

    /// Delta-t cut for pi+
    pub fn dt_pip_cut(&self, part: usize) -> bool {
        let dt = self.dt.get_dt_pi(part);
        let sec = self.event.sc_sect(part);
        if sec == 0 {
            return false;
        }
        let p = self.event.p(part);

        dt <= log_pol2(&DT_PIP_TOP, p) && dt >= log_pol2(&DT_PIP_BOTTOM, p)
    }

    /// Delta-t cut for proton
    pub fn dt_p_cut(&self, part: usize) -> bool {
        let dt = self.dt.get_dt_p(part);
        let sec = self.event.dc_sect(part) - 1;
        if sec < 0 {
            return false;
        }
        let p = self.event.p(part);

        dt <= log_pol2(&DT_P_TOP, p) && dt >= log_pol2(&DT_P_BOT, p)
    }

    /// Delta-t cut for kaon
    pub fn dt_k_cut(&self, part: usize) -> bool {
        let dt = self.dt.get_dt_k(part);
        let sec = self.event.dc_sect(part) - 1;
        if sec < 0 {
            return false;
        }
        let p = self.event.p(part);

        dt <= dt_poly4(&DT_PIP_CONST_TOP, p)
            && dt >= dt_poly4(&DT_PIP_CONST_BOTTOM, p)
            && p < 3.0
    }

    /// Bad SC paddle cut
    pub fn bad_sc_cut(&self, part: usize) -> bool {
        !bad_sc_paddles(self.event.sc_sect(part), self.event.sc_pd(part))
    }

    /// Electron fiducial cut
    pub fn elec_fid_cut(&self) -> bool {
        let c = [0.04, 0.03, 15.0];
        let y = c[0] * self.phi_cent * self.phi_cent + c[1] * self.phi_cent + c[2];
        self.theta >= y
    }
}

/// E1D-specific cuts
pub struct E1dCuts<'a> {
    base: Cuts<'a>,
}

impl<'a> E1dCuts<'a> {
    pub fn new(event: &'a Event) -> Self {
        Self { base: Cuts::new(event) }
    }

    pub fn is_electron(&self) -> bool {
        let mut elec = true;
        elec &= self.base.is_electron();
        elec &= self.beam_cut();

        // Sampling fraction cut
        let sf = self.base.event.etot(0) / self.base.event.p(0);
        elec &= self.sf_cut(sf, self.base.event.p(0));
        elec &= self.base.bad_sc_cut(0);

        if !elec {
            return elec;
        }

        elec &= self.base.fid_chern_cut();

        let sec = self.base.event.dc_sect(0);
        let t = kinematics::theta_calc_rad(self.base.event.cz(0));

        // Sector-specific hand cuts
        if sec == 5 {
            elec &= !(t >= 0.00006 * self.base.phi_cent * self.base.phi_cent + 0.0005 * self.base.phi_cent + 0.59
                && t <= 0.00005 * self.base.phi_cent * self.base.phi_cent + 0.0006 * self.base.phi_cent + 0.65);

            elec &= !(t >= 0.00008 * self.base.phi_cent * self.base.phi_cent + 0.0005 * self.base.phi_cent + 0.456
                && t <= 0.00008 * self.base.phi_cent * self.base.phi_cent + 0.0005 * self.base.phi_cent + 0.48);

            elec &= !(t >= 0.00005 * self.base.phi_cent * self.base.phi_cent + 0.36
                && t <= 0.00005 * self.base.phi_cent * self.base.phi_cent + 0.375);
        }

        if sec == 3 {
            elec &= !(t >= -0.00002 * self.base.phi_cent * self.base.phi_cent + 0.0003 * self.base.phi_cent + 0.38
                && t <= 0.000045 * self.base.phi_cent * self.base.phi_cent + 0.395);
        }

        if sec == 1 {
            elec &= !(t >= -0.000045 * self.base.phi_cent * self.base.phi_cent + 0.345
                && t <= 0.000045 * self.base.phi_cent * self.base.phi_cent + 0.355);
        }

        elec
    }

    pub fn beam_cut(&self) -> bool {
        self.base.event.dc_vx(0) > 0.2
            && self.base.event.dc_vx(0) < 0.4
            && self.base.event.dc_vy(0) > -0.1
            && self.base.event.dc_vy(0) < 0.16
            && self.base.event.dc_vz(0) > -5.0
            && self.base.event.dc_vz(0) < 5.0
    }

    pub fn sf_top_fit(&self, p: f64) -> f64 {
        let par = [0.4128860178888431, -0.014187665284655775, 1.6610245247603538e-05];
        ec_fit_func(p, &par)
    }

    pub fn sf_bot_fit(&self, p: f64) -> f64 {
        let par = [0.08883889850547949, 0.04546759865840269, -4.28522184014871e-05];
        ec_fit_func(p, &par)
    }

    pub fn sf_cut(&self, sf: f32, p: f32) -> bool {
        let sf = sf as f64;
        let p = p as f64;
        sf > self.sf_bot_fit(p) && sf < self.sf_top_fit(p)
    }

    pub fn pip(&self, part: usize) -> bool {
        self.base.event.q(part) == POSITIVE
            && self.base.hadron_fid_arjun(part)
            && self.dt_pip_cut(part)
    }

    pub fn dt_pip_cut(&self, part: usize) -> bool {
        self.base.dt_pip_cut(part) && self.base.bad_sc_cut(part)
    }

    pub fn fid_chern_cut(&self) -> bool {
        let a = -0.000785;
        let b = 0.0;
        let c = -0.00168;
        let d = 1.0;

        let p0_x = self.base.event.dc_xsc(0);
        let p0_y = self.base.event.dc_ysc(0);
        let p0_z = self.base.event.dc_zsc(0);
        let n_x = self.base.event.dc_cxsc(0);
        let n_y = self.base.event.dc_cysc(0);
        let n_z = self.base.event.dc_czsc(0);

        let numer = a * p0_x + b * p0_y + c * p0_z + d;
        let denom = a * n_x + b * n_y + c * n_z;

        let t = (numer / denom).abs();
        let final_x = p0_x + n_x * t;
        let final_y = p0_y + n_y * t;
        let final_z = p0_z + n_z * t;

        let mag = (final_x * final_x + final_y * final_y + final_z * final_z).sqrt();
        let cc_theta = (final_z / mag).acos();
        let cc_phi = final_y.atan2(final_x);

        let cc_x = self.base.event.cc_r(0) * cc_theta.sin() * cc_phi.cos();
        let cc_y = self.base.event.cc_r(0) * cc_theta.sin() * cc_phi.sin();

        fid_chern(cc_x, cc_y)
    }
}

/// E1F-specific cuts
pub struct E1fCuts<'a> {
    base: Cuts<'a>,
}

impl<'a> E1fCuts<'a> {
    pub fn new(event: &'a Event) -> Self {
        Self { base: Cuts::new(event) }
    }

    pub fn is_electron(&self) -> bool {
        let mut elec = true;
        elec &= self.base.is_electron();
        elec &= self.beam_cut();

        if !elec {
            return elec;
        }

        elec &= self.base.fid_chern_cut();
        elec
    }

    pub fn beam_cut(&self) -> bool {
        self.base.event.dc_vx(0).abs() < 0.3
            && self.base.event.dc_vy(0).abs() < 0.4
    }
}

/// E16-specific cuts
pub struct E16Cuts<'a> {
    base: Cuts<'a>,
}

impl<'a> E16Cuts<'a> {
    pub fn new(event: &'a Event) -> Self {
        Self { base: Cuts::new(event) }
    }

    pub fn is_electron(&self) -> bool {
        let mut elec = true;
        elec &= self.base.is_electron();
        elec &= self.beam_cut();

        if !elec {
            return elec;
        }

        elec &= self.base.fid_chern_cut();
        elec
    }

    pub fn beam_cut(&self) -> bool {
        // No beam cut for E16
        true
    }
}
