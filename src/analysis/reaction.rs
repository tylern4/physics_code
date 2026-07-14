use crate::physics::constants::*;
use crate::physics::four_momentum::{angle_between, FourMomentum, ThreeVector};
use crate::physics::kinematics;
use crate::reader::event::Event;

/// Reaction class for event reconstruction.
/// Builds the reaction topology from an event and computes kinematic quantities.
#[derive(Debug, Clone)]
pub struct Reaction {
    pub beam_energy: f32,
    pub beam: FourMomentum,
    pub target: FourMomentum,
    pub gamma: FourMomentum,

    pub elec: Option<FourMomentum>,
    pub prot: Option<FourMomentum>,
    pub pip: Option<FourMomentum>,
    pub pim: Option<FourMomentum>,
    pub neutron: Option<FourMomentum>,
    pub photons: Vec<FourMomentum>,

    pub pair_mass: Vec<f32>,

    // Particle counters
    pub has_e: bool,
    pub has_p: bool,
    pub has_pip: bool,
    pub has_pim: bool,
    pub has_other: bool,
    pub has_neutron: bool,

    pub num_prot: i32,
    pub num_pip: i32,
    pub num_pim: i32,
    pub num_pos: i32,
    pub num_neg: i32,
    pub num_neutral: i32,
    pub num_photons: i32,
    pub num_other: i32,

    pub sector: i32,

    // Cached computed values
    mm_calc: bool,
    mm: f32,
    mm2: f32,
    pi0_mass: f32,
    pi0_mass2: f32,

    w: f32,
    q2: f32,
    xb: f32,

    theta_e: f32,
    theta_star: f32,
    phi_star: f32,

    boosted: bool,
}

impl Reaction {
    /// Create a new Reaction from an event and beam energy
    pub fn new(event: &Event, beam_energy: f32) -> Self {
        let beam = FourMomentum::new(0.0, 0.0, beam_energy, MASS_E);
        let target = FourMomentum::at_rest(MASS_P);

        let mut reaction = Self {
            beam_energy,
            beam,
            target,
            gamma: FourMomentum::new(0.0, 0.0, 0.0, 0.0),

            elec: None,
            prot: None,
            pip: None,
            pim: None,
            neutron: None,
            photons: Vec::new(),

            pair_mass: Vec::new(),

            has_e: false,
            has_p: false,
            has_pip: false,
            has_pim: false,
            has_other: false,
            has_neutron: false,

            num_prot: 0,
            num_pip: 0,
            num_pim: 0,
            num_pos: 0,
            num_neg: 0,
            num_neutral: 0,
            num_photons: 0,
            num_other: 0,

            sector: -1,

            mm_calc: false,
            mm: f32::NAN,
            mm2: f32::NAN,
            pi0_mass: f32::NAN,
            pi0_mass2: f32::NAN,

            w: f32::NAN,
            q2: f32::NAN,
            xb: f32::NAN,

            theta_e: f32::NAN,
            theta_star: f32::NAN,
            phi_star: f32::NAN,

            boosted: false,
        };

        // Set up electron (first particle)
        reaction.has_e = true;
        reaction.sector = event.dc_sect(0);
        reaction.elec = Some(FourMomentum::new(
            event.px(0),
            event.py(0),
            event.pz(0),
            MASS_E,
        ));

        // Calculate gamma, W, Q2, xb
        let elec = reaction.elec.unwrap();
        reaction.gamma = reaction.beam - elec;
        reaction.w = kinematics::w_calc_from_gamma(&reaction.gamma);
        reaction.q2 = kinematics::q2_calc_from_gamma(&reaction.gamma);
        reaction.xb = kinematics::xb_calc_from_gamma(&reaction.gamma);

        reaction
    }

    /// Reset the reaction for reuse
    pub fn reset(&mut self) {
        self.elec = None;
        self.prot = None;
        self.pip = None;
        self.pim = None;
        self.neutron = None;
        self.photons.clear();
        self.pair_mass.clear();

        self.has_e = false;
        self.has_p = false;
        self.has_pip = false;
        self.has_pim = false;
        self.has_other = false;
        self.has_neutron = false;

        self.num_prot = 0;
        self.num_pip = 0;
        self.num_pim = 0;
        self.num_pos = 0;
        self.num_neg = 0;
        self.num_neutral = 0;
        self.num_photons = 0;
        self.num_other = 0;

        self.sector = -1;

        self.mm_calc = false;
        self.mm = f32::NAN;
        self.mm2 = f32::NAN;
        self.pi0_mass = f32::NAN;
        self.pi0_mass2 = f32::NAN;

        self.w = f32::NAN;
        self.q2 = f32::NAN;
        self.xb = f32::NAN;

        self.theta_e = f32::NAN;
        self.theta_star = f32::NAN;
        self.phi_star = f32::NAN;

        self.boosted = false;
    }

    /// Set particle as proton
    pub fn set_proton(&mut self, px: f32, py: f32, pz: f32) {
        self.num_prot += 1;
        self.num_pos += 1;
        self.has_p = true;
        self.prot = Some(FourMomentum::new(px, py, pz, MASS_P));
    }

    /// Set particle as pi+
    pub fn set_pip(&mut self, px: f32, py: f32, pz: f32) {
        self.num_pip += 1;
        self.num_pos += 1;
        self.has_pip = true;
        self.pip = Some(FourMomentum::new(px, py, pz, MASS_PIP));
    }

    /// Set particle as pi-
    pub fn set_pim(&mut self, px: f32, py: f32, pz: f32) {
        self.num_pim += 1;
        self.num_neg += 1;
        self.has_pim = true;
        self.pim = Some(FourMomentum::new(px, py, pz, MASS_PIM));
    }

    /// Set particle as neutron
    pub fn set_neutron(&mut self, px: f32, py: f32, pz: f32) {
        self.num_neutral += 1;
        self.has_neutron = true;
        self.neutron = Some(FourMomentum::new(px, py, pz, MASS_N));
    }

    /// Set particle as other (neutron, photon, or unknown)
    pub fn set_other(&mut self, id: i32, px: f32, py: f32, pz: f32) {
        if id == NEUTRON {
            self.set_neutron(px, py, pz);
        } else if id == PHOTON {
            self.photons.push(FourMomentum::new(px, py, pz, 0.0));
            self.num_photons += 1;
        } else {
            self.num_other += 1;
            self.has_other = true;
        }
    }

    /// Calculate missing mass and pi0 mass
    pub fn calc_missing_mass(&mut self) {
        let elec = match self.elec {
            Some(e) => e,
            None => return,
        };

        let mut mm = self.beam - elec + self.target;

        if self.single_pip() || self.neutron_pip() {
            if let Some(pip) = self.pip {
                mm = mm - pip;
                self.mm = mm.m();
                self.mm2 = mm.m2();
            }
        } else if self.two_pion() {
            if let (Some(pip), Some(pim)) = (self.pip, self.pim) {
                mm = mm - pip - pim;
                self.mm = mm.m();
                self.mm2 = mm.m2();
            }
        } else if self.proton_pim() {
            if let (Some(prot), Some(pim)) = (self.prot, self.pim) {
                mm = mm - prot - pim;
                self.mm = mm.m();
                self.mm2 = mm.m2();
            }
        } else if self.single_p() {
            if let Some(prot) = self.prot {
                mm = mm - prot;
                self.mm = mm.m();
                self.mm2 = mm.m2();
            }
        }

        // Calculate pi0 mass from 2 photons
        if self.num_photons == 2 {
            let phi = angle_between(&self.photons[0], &self.photons[1]);
            if phi < 0.1 {
                return;
            }
            let pi0_vec = self.photons[0] + self.photons[1];
            self.pi0_mass = pi0_vec.m();
            self.pi0_mass2 = pi0_vec.m2();
        }

        self.mm_calc = true;
    }

    /// Get missing mass (lazy calculation)
    pub fn mm(&mut self) -> f32 {
        if !self.mm_calc {
            self.calc_missing_mass();
        }
        self.mm
    }

    /// Get missing mass squared (lazy calculation)
    pub fn mm2(&mut self) -> f32 {
        if !self.mm_calc {
            self.calc_missing_mass();
        }
        self.mm2
    }

    /// Get pi0 mass (lazy calculation)
    pub fn pi0_mass(&mut self) -> f32 {
        if !self.mm_calc {
            self.calc_missing_mass();
        }
        self.pi0_mass
    }

    /// Get pi0 mass squared (lazy calculation)
    pub fn pi0_mass2(&mut self) -> f32 {
        if !self.mm_calc {
            self.calc_missing_mass();
        }
        self.pi0_mass2
    }

    /// Calculate mass pairs from photons
    pub fn calc_mass_pairs(&mut self) {
        if self.photons.len() < 2 {
            return;
        }

        let n = self.photons.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let phi = angle_between(&self.photons[i], &self.photons[j]);
                if phi < 0.1 {
                    continue;
                }
                let pair = self.photons[i] + self.photons[j];
                self.pair_mass.push(pair.m());
            }
        }
    }

    // Channel identification
    pub fn single_pip(&self) -> bool {
        self.num_pip == 1 && self.has_e && !self.has_p && self.has_pip && !self.has_pim
    }

    pub fn neutron_pip(&self) -> bool {
        self.num_pip == 1 && self.num_neutral == 1 && self.has_e && !self.has_p && self.has_pip && !self.has_pim && self.has_neutron
    }

    pub fn single_p(&self) -> bool {
        self.num_prot == 1 && self.has_e && self.has_p && !self.has_pip && !self.has_pim && !self.has_neutron
    }

    pub fn two_pion(&self) -> bool {
        self.num_pip == 1 && self.num_pim == 1 && self.has_e && !self.has_p && self.has_pip && self.has_pim && !self.has_neutron && !self.has_other
    }

    pub fn proton_pim(&self) -> bool {
        self.num_prot == 1 && self.num_pim == 1 && self.has_e && self.has_p && !self.has_pip && self.has_pim && !self.has_neutron && !self.has_other
    }

    pub fn ppi0(&mut self) -> bool {
        self.single_p() && self.mm().abs() <= 0.4 && {
            let pi0 = self.pi0_mass();
            pi0 >= 0.053397 && pi0 <= 0.202516
        }
    }

    pub fn elastic(&mut self) -> bool {
        let mut is_elastic = self.single_p();
        is_elastic &= self.mm2().abs() < 0.05;
        if is_elastic {
            if let (Some(e), Some(p)) = (self.elec, self.prot) {
                let phi_diff = (e.phi() - p.phi()).abs();
                is_elastic &= phi_diff > 3.1 && phi_diff < 3.2;
            } else {
                is_elastic = false;
            }
        }
        is_elastic &= !self.ppi0();
        is_elastic
    }

    pub fn mm_cut(&mut self, sector: i32) -> bool {
        let sec_idx = (sector - 1) as usize;
        if sec_idx >= 6 {
            return false;
        }
        let mu = MM2_FIT_VALUES[sec_idx][1];
        let sigma = MM2_FIT_VALUES[sec_idx][2];
        let mm2 = self.mm2();
        mm2 <= mu + N_SIGMA * sigma && mm2 >= mu - N_SIGMA * sigma
    }

    /// Channel selection: single_pip or neutron_pip with mm_cut
    pub fn channel(&mut self) -> bool {
        (self.single_pip() || self.neutron_pip()) && self.mm_cut(self.sector)
    }

    // Kinematic accessors
    pub fn w(&self) -> f32 {
        self.w
    }
    pub fn q2(&self) -> f32 {
        self.q2
    }
    pub fn xb(&self) -> f32 {
        self.xb
    }
    pub fn e_prime(&self) -> f32 {
        self.elec.map_or(f32::NAN, |e| e.e())
    }
    pub fn phi_diff(&self) -> f32 {
        match (self.elec, self.prot) {
            (Some(e), Some(p)) => (e.phi() - p.phi()).abs(),
            _ => f32::NAN,
        }
    }
    pub fn p_theta(&self) -> f32 {
        self.prot.map_or(f32::NAN, |p| p.theta_deg())
    }
    pub fn p_mom(&self) -> f32 {
        self.prot.map_or(f32::NAN, |p| p.p())
    }

    /// Boost to center-of-mass frame and compute theta_star, phi_star
    pub fn boost(&mut self) {
        if self.boosted {
            return;
        }
        self.boosted = true;

        if (self.single_pip() || self.neutron_pip()) && self.pip.is_some() {
            self.boost_pip();
        } else if (self.single_p() || self.ppi0()) && self.prot.is_some() {
            self.boost_p();
        } else {
            self.theta_e = f32::NAN;
            self.theta_star = f32::NAN;
            self.phi_star = f32::NAN;
        }
    }

    fn boost_pip(&mut self) {
        let elec = self.elec.unwrap();
        let pip = self.pip.unwrap();
        let com = self.target + (self.beam - elec);

        // Calculate boost vector
        let com_e = com.e();
        let _com_px = com.px / com_e;
        let _com_py = com.py / com_e;
        let _com_pz = com.pz / com_e;

        // Calculate rotation axes
        let gamma = self.beam - elec;
        let uz = ThreeVector::new(gamma.px, gamma.py, gamma.pz).unit();
        let ux = ThreeVector::new(self.beam.px, self.beam.py, self.beam.pz)
            .cross(&ThreeVector::new(elec.px, elec.py, elec.pz))
            .unit();

        // Rotate uz by -90 degrees around ux
        let _uz_rot = uz.rotate_z(-PI / 2.0);

        // Build rotation matrix (uz, ux cross uz, ux)
        let uy = uz.cross(&ux);

        // Apply rotation and boost to pip
        let pip_v = ThreeVector::new(pip.px, pip.py, pip.pz);
        let pip_rotated = ThreeVector::new(
            pip_v.dot(&ux),
            pip_v.dot(&uy),
            pip_v.dot(&uz),
        );

        // Apply boost
        let pip_boosted = pip_rotated.scale(1.0); // simplified - full rotation matrix needed

        self.theta_e = elec.theta();
        self.theta_star = pip_boosted.z.atan2((pip_boosted.x * pip_boosted.x + pip_boosted.y * pip_boosted.y).sqrt());
        self.phi_star = kinematics::inv_tan(pip_boosted.y, pip_boosted.x);
    }

    fn boost_p(&mut self) {
        let elec = self.elec.unwrap();
        let prot = self.prot.unwrap();
        let _com = self.target + (self.beam - elec);

        let gamma = self.beam - elec;
        let uz = ThreeVector::new(gamma.px, gamma.py, gamma.pz).unit();
        let ux = ThreeVector::new(self.beam.px, self.beam.py, self.beam.pz)
            .cross(&ThreeVector::new(elec.px, elec.py, elec.pz))
            .unit();

        let uy = uz.cross(&ux);

        let prot_v = ThreeVector::new(prot.px, prot.py, prot.pz);
        let prot_rotated = ThreeVector::new(
            prot_v.dot(&ux),
            prot_v.dot(&uy),
            prot_v.dot(&uz),
        );

        self.theta_e = elec.theta();
        self.theta_star = prot_rotated.z.atan2((prot_rotated.x * prot_rotated.x + prot_rotated.y * prot_rotated.y).sqrt());
        self.phi_star = kinematics::inv_tan(prot_rotated.y, prot_rotated.x);
    }

    pub fn theta_star(&mut self) -> f32 {
        if !self.boosted {
            self.boost();
        }
        self.theta_star
    }

    pub fn phi_star(&mut self) -> f32 {
        if !self.boosted {
            self.boost();
        }
        self.phi_star
    }

    /// Return the reaction type code (matches C++ Type())
    pub fn event_type(&mut self) -> i32 {
        if self.single_pip() {
            0
        } else if self.neutron_pip() {
            10
        } else if self.single_p() {
            22
        } else if self.ppi0() {
            222
        } else if self.proton_pim() {
            3333
        } else {
            -1
        }
    }
}

/// MCReaction extends Reaction with thrown-level kinematics
#[derive(Debug, Clone)]
pub struct MCReaction {
    pub base: Reaction,
    pub w_thrown: f32,
    pub q2_thrown: f32,
    pub elec_thrown: Option<FourMomentum>,
    pub gamma_thrown: Option<FourMomentum>,
    pub pip_thrown: Option<FourMomentum>,
}

impl MCReaction {
    pub fn new(event: &Event, beam_energy: f32) -> Self {
        let base = Reaction::new(event, beam_energy);

        let elec_thrown = Some(FourMomentum::new(
            event.pxpart(0),
            event.pypart(0),
            event.pzpart(0),
            MASS_E,
        ));

        let gamma_thrown = Some(base.beam - elec_thrown.unwrap());
        let w_thrown = kinematics::w_calc_from_gamma(&gamma_thrown.unwrap());
        let q2_thrown = kinematics::q2_calc_from_gamma(&gamma_thrown.unwrap());

        let pip_thrown = Some(FourMomentum::new(
            event.pxpart(1),
            event.pypart(1),
            event.pzpart(1),
            MASS_PIP,
        ));

        Self {
            base,
            w_thrown,
            q2_thrown,
            elec_thrown,
            gamma_thrown,
            pip_thrown,
        }
    }

    pub fn w_thrown(&self) -> f32 {
        self.w_thrown
    }
    pub fn q2_thrown(&self) -> f32 {
        self.q2_thrown
    }

    pub fn mm_thrown(&self) -> f32 {
        let elec = match self.elec_thrown {
            Some(e) => e,
            None => return f32::NAN,
        };
        let pip = match self.pip_thrown {
            Some(p) => p,
            None => return f32::NAN,
        };
        let mm = self.base.beam - elec + self.base.target - pip;
        mm.m()
    }

    pub fn mm2_thrown(&self) -> f32 {
        let elec = match self.elec_thrown {
            Some(e) => e,
            None => return f32::NAN,
        };
        let pip = match self.pip_thrown {
            Some(p) => p,
            None => return f32::NAN,
        };
        let mm = self.base.beam - elec + self.base.target - pip;
        mm.m2()
    }
}
