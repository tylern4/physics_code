use serde::{Deserialize, Serialize};

use crate::physics::constants::MAX_PARTS;

/// Event data matching the CLAS12 h10 TTree schema.
/// Each field corresponds to a branch in the ROOT TTree.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Event {
    // Event-level scalars
    pub npart: i32,
    pub evstat: i32,
    pub intt: i32,
    pub evntid: i32,
    pub evtype: i32,
    pub evntclas: i32,
    pub evthel: i32,
    pub evntclas2: i32,
    pub q_l: f32,
    pub t_l: f32,
    pub tr_time: f32,
    pub rf_time1: f32,
    pub rf_time2: f32,
    pub gpart: i32,

    // Per-particle arrays (indexed by particle number)
    pub id: Vec<i32>,
    pub stat: Vec<i32>,
    pub dc: Vec<i32>,
    pub cc: Vec<i32>,
    pub sc: Vec<i32>,
    pub ec: Vec<i32>,
    pub lec: Vec<i32>,
    pub ccst: Vec<i32>,
    pub p: Vec<f32>,
    pub q: Vec<i32>,
    pub b: Vec<f32>,
    pub cx: Vec<f32>,
    pub cy: Vec<f32>,
    pub cz: Vec<f32>,
    pub vx: Vec<f32>,
    pub vy: Vec<f32>,
    pub vz: Vec<f32>,

    // DC bank
    pub dc_part: i32,
    pub dc_sect: Vec<i32>,
    pub dc_trk: Vec<i32>,
    pub dc_stat: Vec<i32>,
    pub dc_vx: Vec<f32>,
    pub dc_vy: Vec<f32>,
    pub dc_vz: Vec<f32>,
    pub dc_vr: Vec<f32>,
    pub dc_xsc: Vec<f32>,
    pub dc_ysc: Vec<f32>,
    pub dc_zsc: Vec<f32>,
    pub dc_cxsc: Vec<f32>,
    pub dc_cysc: Vec<f32>,
    pub dc_czsc: Vec<f32>,
    pub dc_c2: Vec<f32>,

    // EC bank
    pub ec_part: i32,
    pub ec_stat: Vec<i32>,
    pub ec_sect: Vec<i32>,
    pub ec_whol: Vec<i32>,
    pub ec_inst: Vec<i32>,
    pub ec_oust: Vec<i32>,
    pub etot: Vec<f32>,
    pub ec_ei: Vec<f32>,
    pub ec_eo: Vec<f32>,
    pub ec_t: Vec<f32>,
    pub ec_r: Vec<f32>,
    pub ech_x: Vec<f32>,
    pub ech_y: Vec<f32>,
    pub ech_z: Vec<f32>,
    pub ec_m2: Vec<f32>,
    pub ec_m3: Vec<f32>,
    pub ec_m4: Vec<f32>,
    pub ec_c2: Vec<f32>,

    // SC bank
    pub sc_part: i32,
    pub sc_sect: Vec<i32>,
    pub sc_hit: Vec<i32>,
    pub sc_pd: Vec<i32>,
    pub sc_stat: Vec<i32>,
    pub edep: Vec<f32>,
    pub sc_t: Vec<f32>,
    pub sc_r: Vec<f32>,
    pub sc_c2: Vec<f32>,

    // CC bank
    pub cc_part: i32,
    pub cc_sect: Vec<i32>,
    pub cc_hit: Vec<i32>,
    pub cc_segm: Vec<i32>,
    pub nphe: Vec<i32>,
    pub cc_t: Vec<f32>,
    pub cc_r: Vec<f32>,
    pub cc_c2: Vec<f32>,

    // MC thrown-level branches
    pub nprt: i32,
    pub pidpart: Vec<i32>,
    pub xpart: Vec<f32>,
    pub ypart: Vec<f32>,
    pub zpart: Vec<f32>,
    pub epart: Vec<f32>,
    pub pxpart: Vec<f32>,
    pub pypart: Vec<f32>,
    pub pzpart: Vec<f32>,
    pub qpart: Vec<f32>,
    pub flagspart: Vec<i32>,
}

impl Event {
    /// Create a new empty event with pre-allocated arrays
    pub fn new() -> Self {
        let mut e = Self::default();
        e.id = vec![0; MAX_PARTS];
        e.stat = vec![0; MAX_PARTS];
        e.dc = vec![0; MAX_PARTS];
        e.cc = vec![0; MAX_PARTS];
        e.sc = vec![0; MAX_PARTS];
        e.ec = vec![0; MAX_PARTS];
        e.lec = vec![0; MAX_PARTS];
        e.ccst = vec![0; MAX_PARTS];
        e.p = vec![0.0; MAX_PARTS];
        e.q = vec![0; MAX_PARTS];
        e.b = vec![0.0; MAX_PARTS];
        e.cx = vec![0.0; MAX_PARTS];
        e.cy = vec![0.0; MAX_PARTS];
        e.cz = vec![0.0; MAX_PARTS];
        e.vx = vec![0.0; MAX_PARTS];
        e.vy = vec![0.0; MAX_PARTS];
        e.vz = vec![0.0; MAX_PARTS];

        e.dc_sect = vec![0; MAX_PARTS];
        e.dc_trk = vec![0; MAX_PARTS];
        e.dc_stat = vec![0; MAX_PARTS];
        e.dc_vx = vec![0.0; MAX_PARTS];
        e.dc_vy = vec![0.0; MAX_PARTS];
        e.dc_vz = vec![0.0; MAX_PARTS];
        e.dc_vr = vec![0.0; MAX_PARTS];
        e.dc_xsc = vec![0.0; MAX_PARTS];
        e.dc_ysc = vec![0.0; MAX_PARTS];
        e.dc_zsc = vec![0.0; MAX_PARTS];
        e.dc_cxsc = vec![0.0; MAX_PARTS];
        e.dc_cysc = vec![0.0; MAX_PARTS];
        e.dc_czsc = vec![0.0; MAX_PARTS];
        e.dc_c2 = vec![0.0; MAX_PARTS];

        e.ec_stat = vec![0; MAX_PARTS];
        e.ec_sect = vec![0; MAX_PARTS];
        e.ec_whol = vec![0; MAX_PARTS];
        e.ec_inst = vec![0; MAX_PARTS];
        e.ec_oust = vec![0; MAX_PARTS];
        e.etot = vec![0.0; MAX_PARTS];
        e.ec_ei = vec![0.0; MAX_PARTS];
        e.ec_eo = vec![0.0; MAX_PARTS];
        e.ec_t = vec![0.0; MAX_PARTS];
        e.ec_r = vec![0.0; MAX_PARTS];
        e.ech_x = vec![0.0; MAX_PARTS];
        e.ech_y = vec![0.0; MAX_PARTS];
        e.ech_z = vec![0.0; MAX_PARTS];
        e.ec_m2 = vec![0.0; MAX_PARTS];
        e.ec_m3 = vec![0.0; MAX_PARTS];
        e.ec_m4 = vec![0.0; MAX_PARTS];
        e.ec_c2 = vec![0.0; MAX_PARTS];

        e.sc_sect = vec![0; MAX_PARTS];
        e.sc_hit = vec![0; MAX_PARTS];
        e.sc_pd = vec![0; MAX_PARTS];
        e.sc_stat = vec![0; MAX_PARTS];
        e.edep = vec![0.0; MAX_PARTS];
        e.sc_t = vec![0.0; MAX_PARTS];
        e.sc_r = vec![0.0; MAX_PARTS];
        e.sc_c2 = vec![0.0; MAX_PARTS];

        e.cc_sect = vec![0; MAX_PARTS];
        e.cc_hit = vec![0; MAX_PARTS];
        e.cc_segm = vec![0; MAX_PARTS];
        e.nphe = vec![0; MAX_PARTS];
        e.cc_t = vec![0.0; MAX_PARTS];
        e.cc_r = vec![0.0; MAX_PARTS];
        e.cc_c2 = vec![0.0; MAX_PARTS];

        e.pidpart = vec![0; MAX_PARTS];
        e.xpart = vec![0.0; MAX_PARTS];
        e.ypart = vec![0.0; MAX_PARTS];
        e.zpart = vec![0.0; MAX_PARTS];
        e.epart = vec![0.0; MAX_PARTS];
        e.pxpart = vec![0.0; MAX_PARTS];
        e.pypart = vec![0.0; MAX_PARTS];
        e.pzpart = vec![0.0; MAX_PARTS];
        e.qpart = vec![0.0; MAX_PARTS];
        e.flagspart = vec![0; MAX_PARTS];

        e
    }

    // Accessor methods matching the C++ Branches class interface
    pub fn id(&self, i: usize) -> i32 {
        self.id[i]
    }
    pub fn stat(&self, i: usize) -> i32 {
        self.stat[i]
    }
    pub fn dc(&self, i: usize) -> i32 {
        self.dc[i]
    }
    pub fn cc(&self, i: usize) -> i32 {
        self.cc[i]
    }
    pub fn sc(&self, i: usize) -> i32 {
        self.sc[i]
    }
    pub fn ec(&self, i: usize) -> i32 {
        self.ec[i]
    }
    pub fn lec(&self, i: usize) -> i32 {
        self.lec[i]
    }
    pub fn ccst(&self, i: usize) -> i32 {
        self.ccst[i]
    }
    pub fn p(&self, i: usize) -> f32 {
        self.p[i]
    }
    pub fn q(&self, i: usize) -> i32 {
        self.q[i]
    }
    pub fn b(&self, i: usize) -> f32 {
        self.b[i]
    }
    pub fn cx(&self, i: usize) -> f32 {
        self.cx[i]
    }
    pub fn cy(&self, i: usize) -> f32 {
        self.cy[i]
    }
    pub fn cz(&self, i: usize) -> f32 {
        self.cz[i]
    }
    pub fn vx(&self, i: usize) -> f32 {
        self.vx[i]
    }
    pub fn vy(&self, i: usize) -> f32 {
        self.vy[i]
    }
    pub fn vz(&self, i: usize) -> f32 {
        self.vz[i]
    }

    // DC bank - accessed via dc[i]-1 indirection
    pub fn dc_sect(&self, i: usize) -> i32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_sect.len() {
            self.dc_sect[idx as usize]
        } else {
            0
        }
    }
    pub fn dc_trk(&self, i: usize) -> i32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_trk.len() {
            self.dc_trk[idx as usize]
        } else {
            0
        }
    }
    pub fn dc_stat(&self, i: usize) -> i32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_stat.len() {
            self.dc_stat[idx as usize]
        } else {
            0
        }
    }
    pub fn dc_vx(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_vx.len() {
            let val = self.dc_vx[idx as usize];
            if val != 0.0 { val } else { f32::NAN }
        } else {
            f32::NAN
        }
    }
    pub fn dc_vy(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_vy.len() {
            let val = self.dc_vy[idx as usize];
            if val != 0.0 { val } else { f32::NAN }
        } else {
            f32::NAN
        }
    }
    pub fn dc_vz(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_vz.len() {
            let val = self.dc_vz[idx as usize];
            if val != 0.0 { val } else { f32::NAN }
        } else {
            f32::NAN
        }
    }
    pub fn dc_vr(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_vr.len() {
            let val = self.dc_vr[idx as usize];
            if val != 0.0 { val } else { f32::NAN }
        } else {
            f32::NAN
        }
    }
    pub fn dc_xsc(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_xsc.len() {
            self.dc_xsc[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn dc_ysc(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_ysc.len() {
            self.dc_ysc[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn dc_zsc(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_zsc.len() {
            self.dc_zsc[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn dc_cxsc(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_cxsc.len() {
            self.dc_cxsc[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn dc_cysc(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_cysc.len() {
            self.dc_cysc[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn dc_czsc(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_czsc.len() {
            self.dc_czsc[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn dc_c2(&self, i: usize) -> f32 {
        let idx = self.dc[i] - 1;
        if idx >= 0 && (idx as usize) < self.dc_c2.len() {
            self.dc_c2[idx as usize]
        } else {
            f32::NAN
        }
    }

    // EC bank - accessed via ec[i]-1 indirection
    pub fn ec_stat(&self, i: usize) -> i32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_stat.len() {
            self.ec_stat[idx as usize]
        } else {
            0
        }
    }
    pub fn ec_sect(&self, i: usize) -> i32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_sect.len() {
            self.ec_sect[idx as usize]
        } else {
            0
        }
    }
    pub fn ec_whol(&self, i: usize) -> i32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_whol.len() {
            self.ec_whol[idx as usize]
        } else {
            0
        }
    }
    pub fn ec_inst(&self, i: usize) -> i32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_inst.len() {
            self.ec_inst[idx as usize]
        } else {
            0
        }
    }
    pub fn ec_oust(&self, i: usize) -> i32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_oust.len() {
            self.ec_oust[idx as usize]
        } else {
            0
        }
    }
    pub fn etot(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.etot.len() {
            self.etot[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_ei(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_ei.len() {
            self.ec_ei[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_eo(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_eo.len() {
            self.ec_eo[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_t(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_t.len() {
            self.ec_t[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_r(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_r.len() {
            self.ec_r[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ech_x(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ech_x.len() {
            self.ech_x[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ech_y(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ech_y.len() {
            self.ech_y[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ech_z(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ech_z.len() {
            self.ech_z[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_m2(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_m2.len() {
            self.ec_m2[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_m3(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_m3.len() {
            self.ec_m3[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_m4(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_m4.len() {
            self.ec_m4[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn ec_c2(&self, i: usize) -> f32 {
        let idx = self.ec[i] - 1;
        if idx >= 0 && (idx as usize) < self.ec_c2.len() {
            self.ec_c2[idx as usize]
        } else {
            f32::NAN
        }
    }

    // SC bank - accessed via sc[i]-1 indirection
    pub fn sc_sect(&self, i: usize) -> i32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_sect.len() {
            self.sc_sect[idx as usize]
        } else {
            0
        }
    }
    pub fn sc_hit(&self, i: usize) -> i32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_hit.len() {
            self.sc_hit[idx as usize]
        } else {
            0
        }
    }
    pub fn sc_pd(&self, i: usize) -> i32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_pd.len() {
            self.sc_pd[idx as usize]
        } else {
            0
        }
    }
    pub fn sc_stat(&self, i: usize) -> i32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_stat.len() {
            self.sc_stat[idx as usize]
        } else {
            0
        }
    }
    pub fn edep(&self, i: usize) -> f32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.edep.len() {
            self.edep[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn sc_t(&self, i: usize) -> f32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_t.len() {
            self.sc_t[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn sc_r(&self, i: usize) -> f32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_r.len() {
            self.sc_r[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn sc_c2(&self, i: usize) -> f32 {
        let idx = self.sc[i] - 1;
        if idx >= 0 && (idx as usize) < self.sc_c2.len() {
            self.sc_c2[idx as usize]
        } else {
            f32::NAN
        }
    }

    // CC bank - accessed via cc[i]-1 indirection
    pub fn cc_sect(&self, i: usize) -> i32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.cc_sect.len() {
            self.cc_sect[idx as usize]
        } else {
            0
        }
    }
    pub fn cc_hit(&self, i: usize) -> i32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.cc_hit.len() {
            self.cc_hit[idx as usize]
        } else {
            0
        }
    }
    pub fn cc_segm(&self, i: usize) -> i32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.cc_segm.len() {
            self.cc_segm[idx as usize]
        } else {
            0
        }
    }
    pub fn nphe(&self, i: usize) -> i32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.nphe.len() {
            self.nphe[idx as usize]
        } else {
            0
        }
    }
    pub fn cc_t(&self, i: usize) -> f32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.cc_t.len() {
            self.cc_t[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn cc_r(&self, i: usize) -> f32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.cc_r.len() {
            self.cc_r[idx as usize]
        } else {
            f32::NAN
        }
    }
    pub fn cc_c2(&self, i: usize) -> f32 {
        let idx = self.cc[i] - 1;
        if idx >= 0 && (idx as usize) < self.cc_c2.len() {
            self.cc_c2[idx as usize]
        } else {
            f32::NAN
        }
    }

    pub fn cc_x(&self, i: usize) -> f32 {
        self.cc_r[i]
    }
    pub fn cc_y(&self, i: usize) -> f32 {
        self.cc_r[i]
    }

    // Computed momentum components (p * direction cosine)
    pub fn px(&self, i: usize) -> f32 {
        self.p[i] * self.cx[i]
    }
    pub fn py(&self, i: usize) -> f32 {
        self.p[i] * self.cy[i]
    }
    pub fn pz(&self, i: usize) -> f32 {
        self.p[i] * self.cz[i]
    }

    // MC thrown-level accessors
    pub fn pidpart(&self, i: usize) -> i32 {
        self.pidpart[i]
    }
    pub fn pxpart(&self, i: usize) -> f32 {
        self.pxpart[i]
    }
    pub fn pypart(&self, i: usize) -> f32 {
        self.pypart[i]
    }
    pub fn pzpart(&self, i: usize) -> f32 {
        self.pzpart[i]
    }
    pub fn epart(&self, i: usize) -> f32 {
        self.epart[i]
    }
}
