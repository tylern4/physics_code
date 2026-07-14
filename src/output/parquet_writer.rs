use anyhow::Result;
use arrow::array::{BooleanArray, Float32Array, Int32Array, RecordBatch};
use arrow::datatypes::{DataType, Field, Schema};
use parquet::arrow::ArrowWriter;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::sync::Arc;

/// Per-particle row data
#[derive(Debug, Clone, Default)]
pub struct ParticleRow {
    pub part_index: i32,
    pub part_id: i32,
    pub part_q: i32,
    pub part_p: f32,
    pub part_beta: f32,
    pub part_theta: f32,
    pub part_phi: f32,
    pub part_sector: i32,
    pub part_sc_sector: i32,
    pub part_sc_pd: i32,
    pub part_delta_t_p: f32,
    pub part_delta_t_pip: f32,
    pub part_delta_t_e: f32,
    pub part_delta_t_k: f32,
    pub part_nphe: i32,
    pub part_cc_segm: i32,
    pub part_etot: f32,
    pub part_ec_ei: f32,
    pub part_ec_eo: f32,
    pub part_dc_xsc: f32,
    pub part_dc_ysc: f32,
    pub part_dc_zsc: f32,
    pub part_edep: f32,
    pub part_vx: f32,
    pub part_vy: f32,
    pub part_vz: f32,
    pub part_fid_chern: bool,
    pub part_elec_fid: bool,
    pub part_hadron_fid: bool,
    pub part_dt_p_pass: bool,
    pub part_dt_pip_pass: bool,
    pub part_is_pip: bool,
    pub part_is_prot: bool,
    pub part_is_pim: bool,
    pub part_is_electron: bool,
}

/// Analysis result for a single event - one row per particle (long format)
#[derive(Debug, Clone)]
pub struct EventRow {
    // Event-level (repeated for each particle)
    pub w: f32,
    pub q2: f32,
    pub xb: f32,
    pub theta_star: f32,
    pub phi_star: f32,
    pub mm: f32,
    pub mm2: f32,
    pub e_sector: i32,
    pub event_type: i32,
    pub beam_energy: f32,
    pub e_prime: f32,
    pub w_thrown: f32,
    pub q2_thrown: f32,
    pub mm_thrown: f32,
    pub mm2_thrown: f32,
    pub e_dc_vx: f32,
    pub e_dc_vy: f32,
    pub e_dc_vz: f32,
    pub num_pip: i32,
    pub num_prot: i32,
    pub num_pim: i32,
    pub num_pos: i32,
    pub num_neg: i32,
    pub num_neutral: i32,
    pub num_photons: i32,
    // Particle-level
    pub particles: Vec<ParticleRow>,
}

/// Accumulates all event data for writing to parquet
#[derive(Debug, Clone)]
pub struct AnalysisResult {
    // Event-level columns
    pub w: Vec<f32>,
    pub q2: Vec<f32>,
    pub xb: Vec<f32>,
    pub theta_star: Vec<f32>,
    pub phi_star: Vec<f32>,
    pub mm: Vec<f32>,
    pub mm2: Vec<f32>,
    pub e_sector: Vec<i32>,
    pub event_type: Vec<i32>,
    pub beam_energy: Vec<f32>,
    pub e_prime: Vec<f32>,
    pub w_thrown: Vec<f32>,
    pub q2_thrown: Vec<f32>,
    pub mm_thrown: Vec<f32>,
    pub mm2_thrown: Vec<f32>,
    pub e_dc_vx: Vec<f32>,
    pub e_dc_vy: Vec<f32>,
    pub e_dc_vz: Vec<f32>,
    pub num_pip: Vec<i32>,
    pub num_prot: Vec<i32>,
    pub num_pim: Vec<i32>,
    pub num_pos: Vec<i32>,
    pub num_neg: Vec<i32>,
    pub num_neutral: Vec<i32>,
    pub num_photons: Vec<i32>,
    pub pi0_mass: Vec<f32>,
    pub pi0_mass2: Vec<f32>,
    // Particle-level columns
    pub part_index: Vec<i32>,
    pub part_id: Vec<i32>,
    pub part_q: Vec<i32>,
    pub part_p: Vec<f32>,
    pub part_beta: Vec<f32>,
    pub part_theta: Vec<f32>,
    pub part_phi: Vec<f32>,
    pub part_sector: Vec<i32>,
    pub part_sc_sector: Vec<i32>,
    pub part_sc_pd: Vec<i32>,
    pub part_delta_t_p: Vec<f32>,
    pub part_delta_t_pip: Vec<f32>,
    pub part_delta_t_e: Vec<f32>,
    pub part_delta_t_k: Vec<f32>,
    pub part_nphe: Vec<i32>,
    pub part_cc_segm: Vec<i32>,
    pub part_etot: Vec<f32>,
    pub part_ec_ei: Vec<f32>,
    pub part_ec_eo: Vec<f32>,
    pub part_dc_xsc: Vec<f32>,
    pub part_dc_ysc: Vec<f32>,
    pub part_dc_zsc: Vec<f32>,
    pub part_edep: Vec<f32>,
    pub part_vx: Vec<f32>,
    pub part_vy: Vec<f32>,
    pub part_vz: Vec<f32>,
    pub part_fid_chern: Vec<bool>,
    pub part_elec_fid: Vec<bool>,
    pub part_hadron_fid: Vec<bool>,
    pub part_dt_p_pass: Vec<bool>,
    pub part_dt_pip_pass: Vec<bool>,
    pub part_is_pip: Vec<bool>,
    pub part_is_prot: Vec<bool>,
    pub part_is_pim: Vec<bool>,
    pub part_is_electron: Vec<bool>,
    pub n_events: usize,
    pub n_particles: usize,
}

impl AnalysisResult {
    pub fn new() -> Self {
        Self {
            w: Vec::new(),
            q2: Vec::new(),
            xb: Vec::new(),
            theta_star: Vec::new(),
            phi_star: Vec::new(),
            mm: Vec::new(),
            mm2: Vec::new(),
            e_sector: Vec::new(),
            event_type: Vec::new(),
            beam_energy: Vec::new(),
            e_prime: Vec::new(),
            w_thrown: Vec::new(),
            q2_thrown: Vec::new(),
            mm_thrown: Vec::new(),
            mm2_thrown: Vec::new(),
            e_dc_vx: Vec::new(),
            e_dc_vy: Vec::new(),
            e_dc_vz: Vec::new(),
            num_pip: Vec::new(),
            num_prot: Vec::new(),
            num_pim: Vec::new(),
            num_pos: Vec::new(),
            num_neg: Vec::new(),
            num_neutral: Vec::new(),
            num_photons: Vec::new(),
            pi0_mass: Vec::new(),
            pi0_mass2: Vec::new(),
            part_index: Vec::new(),
            part_id: Vec::new(),
            part_q: Vec::new(),
            part_p: Vec::new(),
            part_beta: Vec::new(),
            part_theta: Vec::new(),
            part_phi: Vec::new(),
            part_sector: Vec::new(),
            part_sc_sector: Vec::new(),
            part_sc_pd: Vec::new(),
            part_delta_t_p: Vec::new(),
            part_delta_t_pip: Vec::new(),
            part_delta_t_e: Vec::new(),
            part_delta_t_k: Vec::new(),
            part_nphe: Vec::new(),
            part_cc_segm: Vec::new(),
            part_etot: Vec::new(),
            part_ec_ei: Vec::new(),
            part_ec_eo: Vec::new(),
            part_dc_xsc: Vec::new(),
            part_dc_ysc: Vec::new(),
            part_dc_zsc: Vec::new(),
            part_edep: Vec::new(),
            part_vx: Vec::new(),
            part_vy: Vec::new(),
            part_vz: Vec::new(),
            part_fid_chern: Vec::new(),
            part_elec_fid: Vec::new(),
            part_hadron_fid: Vec::new(),
            part_dt_p_pass: Vec::new(),
            part_dt_pip_pass: Vec::new(),
            part_is_pip: Vec::new(),
            part_is_prot: Vec::new(),
            part_is_pim: Vec::new(),
            part_is_electron: Vec::new(),
            n_events: 0,
            n_particles: 0,
        }
    }

    /// Push a single event with its particles
    pub fn push_event(
        &mut self,
        w: f32,
        q2: f32,
        xb: f32,
        theta_star: f32,
        phi_star: f32,
        mm: f32,
        mm2: f32,
        e_sector: i32,
        event_type: i32,
        beam_energy: f32,
        e_prime: f32,
        w_thrown: f32,
        q2_thrown: f32,
        mm_thrown: f32,
        mm2_thrown: f32,
        e_dc_vx: f32,
        e_dc_vy: f32,
        e_dc_vz: f32,
        num_pip: i32,
        num_prot: i32,
        num_pim: i32,
        num_pos: i32,
        num_neg: i32,
        num_neutral: i32,
        num_photons: i32,
        pi0_mass: f32,
        pi0_mass2: f32,
        particles: &[ParticleRow],
    ) {
        for part in particles {
            self.w.push(w);
            self.q2.push(q2);
            self.xb.push(xb);
            self.theta_star.push(theta_star);
            self.phi_star.push(phi_star);
            self.mm.push(mm);
            self.mm2.push(mm2);
            self.e_sector.push(e_sector);
            self.event_type.push(event_type);
            self.beam_energy.push(beam_energy);
            self.e_prime.push(e_prime);
            self.w_thrown.push(w_thrown);
            self.q2_thrown.push(q2_thrown);
            self.mm_thrown.push(mm_thrown);
            self.mm2_thrown.push(mm2_thrown);
            self.e_dc_vx.push(e_dc_vx);
            self.e_dc_vy.push(e_dc_vy);
            self.e_dc_vz.push(e_dc_vz);
            self.num_pip.push(num_pip);
            self.num_prot.push(num_prot);
            self.num_pim.push(num_pim);
            self.num_pos.push(num_pos);
            self.num_neg.push(num_neg);
            self.num_neutral.push(num_neutral);
            self.num_photons.push(num_photons);
            self.pi0_mass.push(pi0_mass);
            self.pi0_mass2.push(pi0_mass2);

            self.part_index.push(part.part_index);
            self.part_id.push(part.part_id);
            self.part_q.push(part.part_q);
            self.part_p.push(part.part_p);
            self.part_beta.push(part.part_beta);
            self.part_theta.push(part.part_theta);
            self.part_phi.push(part.part_phi);
            self.part_sector.push(part.part_sector);
            self.part_sc_sector.push(part.part_sc_sector);
            self.part_sc_pd.push(part.part_sc_pd);
            self.part_delta_t_p.push(part.part_delta_t_p);
            self.part_delta_t_pip.push(part.part_delta_t_pip);
            self.part_delta_t_e.push(part.part_delta_t_e);
            self.part_delta_t_k.push(part.part_delta_t_k);
            self.part_nphe.push(part.part_nphe);
            self.part_cc_segm.push(part.part_cc_segm);
            self.part_etot.push(part.part_etot);
            self.part_ec_ei.push(part.part_ec_ei);
            self.part_ec_eo.push(part.part_ec_eo);
            self.part_dc_xsc.push(part.part_dc_xsc);
            self.part_dc_ysc.push(part.part_dc_ysc);
            self.part_dc_zsc.push(part.part_dc_zsc);
            self.part_edep.push(part.part_edep);
            self.part_vx.push(part.part_vx);
            self.part_vy.push(part.part_vy);
            self.part_vz.push(part.part_vz);
            self.part_fid_chern.push(part.part_fid_chern);
            self.part_elec_fid.push(part.part_elec_fid);
            self.part_hadron_fid.push(part.part_hadron_fid);
            self.part_dt_p_pass.push(part.part_dt_p_pass);
            self.part_dt_pip_pass.push(part.part_dt_pip_pass);
            self.part_is_pip.push(part.part_is_pip);
            self.part_is_prot.push(part.part_is_prot);
            self.part_is_pim.push(part.part_is_pim);
            self.part_is_electron.push(part.part_is_electron);

            self.n_particles += 1;
        }
        self.n_events += 1;
    }

    /// Merge another AnalysisResult into this one
    pub fn merge(&mut self, other: AnalysisResult) {
        self.w.extend(other.w);
        self.q2.extend(other.q2);
        self.xb.extend(other.xb);
        self.theta_star.extend(other.theta_star);
        self.phi_star.extend(other.phi_star);
        self.mm.extend(other.mm);
        self.mm2.extend(other.mm2);
        self.e_sector.extend(other.e_sector);
        self.event_type.extend(other.event_type);
        self.beam_energy.extend(other.beam_energy);
        self.e_prime.extend(other.e_prime);
        self.w_thrown.extend(other.w_thrown);
        self.q2_thrown.extend(other.q2_thrown);
        self.mm_thrown.extend(other.mm_thrown);
        self.mm2_thrown.extend(other.mm2_thrown);
        self.e_dc_vx.extend(other.e_dc_vx);
        self.e_dc_vy.extend(other.e_dc_vy);
        self.e_dc_vz.extend(other.e_dc_vz);
        self.num_pip.extend(other.num_pip);
        self.num_prot.extend(other.num_prot);
        self.num_pim.extend(other.num_pim);
        self.num_pos.extend(other.num_pos);
        self.num_neg.extend(other.num_neg);
        self.num_neutral.extend(other.num_neutral);
        self.num_photons.extend(other.num_photons);
        self.pi0_mass.extend(other.pi0_mass);
        self.pi0_mass2.extend(other.pi0_mass2);
        self.part_index.extend(other.part_index);
        self.part_id.extend(other.part_id);
        self.part_q.extend(other.part_q);
        self.part_p.extend(other.part_p);
        self.part_beta.extend(other.part_beta);
        self.part_theta.extend(other.part_theta);
        self.part_phi.extend(other.part_phi);
        self.part_sector.extend(other.part_sector);
        self.part_sc_sector.extend(other.part_sc_sector);
        self.part_sc_pd.extend(other.part_sc_pd);
        self.part_delta_t_p.extend(other.part_delta_t_p);
        self.part_delta_t_pip.extend(other.part_delta_t_pip);
        self.part_delta_t_e.extend(other.part_delta_t_e);
        self.part_delta_t_k.extend(other.part_delta_t_k);
        self.part_nphe.extend(other.part_nphe);
        self.part_cc_segm.extend(other.part_cc_segm);
        self.part_etot.extend(other.part_etot);
        self.part_ec_ei.extend(other.part_ec_ei);
        self.part_ec_eo.extend(other.part_ec_eo);
        self.part_dc_xsc.extend(other.part_dc_xsc);
        self.part_dc_ysc.extend(other.part_dc_ysc);
        self.part_dc_zsc.extend(other.part_dc_zsc);
        self.part_edep.extend(other.part_edep);
        self.part_vx.extend(other.part_vx);
        self.part_vy.extend(other.part_vy);
        self.part_vz.extend(other.part_vz);
        self.part_fid_chern.extend(other.part_fid_chern);
        self.part_elec_fid.extend(other.part_elec_fid);
        self.part_hadron_fid.extend(other.part_hadron_fid);
        self.part_dt_p_pass.extend(other.part_dt_p_pass);
        self.part_dt_pip_pass.extend(other.part_dt_pip_pass);
        self.part_is_pip.extend(other.part_is_pip);
        self.part_is_prot.extend(other.part_is_prot);
        self.part_is_pim.extend(other.part_is_pim);
        self.part_is_electron.extend(other.part_is_electron);
        self.n_events += other.n_events;
        self.n_particles += other.n_particles;
    }
}

fn make_schema() -> Arc<Schema> {
    Arc::new(Schema::new(vec![
        // Event-level
        Field::new("w", DataType::Float32, false),
        Field::new("q2", DataType::Float32, false),
        Field::new("xb", DataType::Float32, false),
        Field::new("theta_star", DataType::Float32, false),
        Field::new("phi_star", DataType::Float32, false),
        Field::new("mm", DataType::Float32, false),
        Field::new("mm2", DataType::Float32, false),
        Field::new("e_sector", DataType::Int32, false),
        Field::new("event_type", DataType::Int32, false),
        Field::new("beam_energy", DataType::Float32, false),
        Field::new("e_prime", DataType::Float32, false),
        Field::new("w_thrown", DataType::Float32, false),
        Field::new("q2_thrown", DataType::Float32, false),
        Field::new("mm_thrown", DataType::Float32, false),
        Field::new("mm2_thrown", DataType::Float32, false),
        Field::new("e_dc_vx", DataType::Float32, false),
        Field::new("e_dc_vy", DataType::Float32, false),
        Field::new("e_dc_vz", DataType::Float32, false),
        Field::new("num_pip", DataType::Int32, false),
        Field::new("num_prot", DataType::Int32, false),
        Field::new("num_pim", DataType::Int32, false),
        Field::new("num_pos", DataType::Int32, false),
        Field::new("num_neg", DataType::Int32, false),
        Field::new("num_neutral", DataType::Int32, false),
        Field::new("num_photons", DataType::Int32, false),
        Field::new("pi0_mass", DataType::Float32, false),
        Field::new("pi0_mass2", DataType::Float32, false),
        // Particle-level
        Field::new("part_index", DataType::Int32, false),
        Field::new("part_id", DataType::Int32, false),
        Field::new("part_q", DataType::Int32, false),
        Field::new("part_p", DataType::Float32, false),
        Field::new("part_beta", DataType::Float32, false),
        Field::new("part_theta", DataType::Float32, false),
        Field::new("part_phi", DataType::Float32, false),
        Field::new("part_sector", DataType::Int32, false),
        Field::new("part_sc_sector", DataType::Int32, false),
        Field::new("part_sc_pd", DataType::Int32, false),
        Field::new("part_delta_t_p", DataType::Float32, false),
        Field::new("part_delta_t_pip", DataType::Float32, false),
        Field::new("part_delta_t_e", DataType::Float32, false),
        Field::new("part_delta_t_k", DataType::Float32, false),
        Field::new("part_nphe", DataType::Int32, false),
        Field::new("part_cc_segm", DataType::Int32, false),
        Field::new("part_etot", DataType::Float32, false),
        Field::new("part_ec_ei", DataType::Float32, false),
        Field::new("part_ec_eo", DataType::Float32, false),
        Field::new("part_dc_xsc", DataType::Float32, false),
        Field::new("part_dc_ysc", DataType::Float32, false),
        Field::new("part_dc_zsc", DataType::Float32, false),
        Field::new("part_edep", DataType::Float32, false),
        Field::new("part_vx", DataType::Float32, false),
        Field::new("part_vy", DataType::Float32, false),
        Field::new("part_vz", DataType::Float32, false),
        Field::new("part_fid_chern", DataType::Boolean, false),
        Field::new("part_elec_fid", DataType::Boolean, false),
        Field::new("part_hadron_fid", DataType::Boolean, false),
        Field::new("part_dt_p_pass", DataType::Boolean, false),
        Field::new("part_dt_pip_pass", DataType::Boolean, false),
        Field::new("part_is_pip", DataType::Boolean, false),
        Field::new("part_is_prot", DataType::Boolean, false),
        Field::new("part_is_pim", DataType::Boolean, false),
        Field::new("part_is_electron", DataType::Boolean, false),
    ]))
}

/// Write analysis results to a CSV file
pub fn write_csv(result: &AnalysisResult, filename: &str) -> Result<()> {
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "w,q2,xb,theta_star,phi_star,mm,mm2,e_sector,event_type,beam_energy,e_prime,\
         w_thrown,q2_thrown,mm_thrown,mm2_thrown,e_dc_vx,e_dc_vy,e_dc_vz,\
         num_pip,num_prot,num_pim,num_pos,num_neg,num_neutral,num_photons,pi0_mass,pi0_mass2,\
         part_index,part_id,part_q,part_p,part_beta,part_theta,part_phi,\
         part_sector,part_sc_sector,part_sc_pd,\
         part_delta_t_p,part_delta_t_pip,part_delta_t_e,part_delta_t_k,\
         part_nphe,part_cc_segm,part_etot,part_ec_ei,part_ec_eo,\
         part_dc_xsc,part_dc_ysc,part_dc_zsc,part_edep,part_vx,part_vy,part_vz,\
         part_fid_chern,part_elec_fid,part_hadron_fid,part_dt_p_pass,part_dt_pip_pass,\
         part_is_pip,part_is_prot,part_is_pim,part_is_electron"
    )?;

    for i in 0..result.n_particles {
        let fields = [
            result.w[i].to_string(),
            result.q2[i].to_string(),
            result.xb[i].to_string(),
            result.theta_star[i].to_string(),
            result.phi_star[i].to_string(),
            result.mm[i].to_string(),
            result.mm2[i].to_string(),
            result.e_sector[i].to_string(),
            result.event_type[i].to_string(),
            result.beam_energy[i].to_string(),
            result.e_prime[i].to_string(),
            result.w_thrown[i].to_string(),
            result.q2_thrown[i].to_string(),
            result.mm_thrown[i].to_string(),
            result.mm2_thrown[i].to_string(),
            result.e_dc_vx[i].to_string(),
            result.e_dc_vy[i].to_string(),
            result.e_dc_vz[i].to_string(),
            result.num_pip[i].to_string(),
            result.num_prot[i].to_string(),
            result.num_pim[i].to_string(),
            result.num_pos[i].to_string(),
            result.num_neg[i].to_string(),
            result.num_neutral[i].to_string(),
            result.num_photons[i].to_string(),
            result.pi0_mass[i].to_string(),
            result.pi0_mass2[i].to_string(),
            result.part_index[i].to_string(),
            result.part_id[i].to_string(),
            result.part_q[i].to_string(),
            result.part_p[i].to_string(),
            result.part_beta[i].to_string(),
            result.part_theta[i].to_string(),
            result.part_phi[i].to_string(),
            result.part_sector[i].to_string(),
            result.part_sc_sector[i].to_string(),
            result.part_sc_pd[i].to_string(),
            result.part_delta_t_p[i].to_string(),
            result.part_delta_t_pip[i].to_string(),
            result.part_delta_t_e[i].to_string(),
            result.part_delta_t_k[i].to_string(),
            result.part_nphe[i].to_string(),
            result.part_cc_segm[i].to_string(),
            result.part_etot[i].to_string(),
            result.part_ec_ei[i].to_string(),
            result.part_ec_eo[i].to_string(),
            result.part_dc_xsc[i].to_string(),
            result.part_dc_ysc[i].to_string(),
            result.part_dc_zsc[i].to_string(),
            result.part_edep[i].to_string(),
            result.part_vx[i].to_string(),
            result.part_vy[i].to_string(),
            result.part_vz[i].to_string(),
            result.part_fid_chern[i].to_string(),
            result.part_elec_fid[i].to_string(),
            result.part_hadron_fid[i].to_string(),
            result.part_dt_p_pass[i].to_string(),
            result.part_dt_pip_pass[i].to_string(),
            result.part_is_pip[i].to_string(),
            result.part_is_prot[i].to_string(),
            result.part_is_pim[i].to_string(),
            result.part_is_electron[i].to_string(),
        ];
        writeln!(writer, "{}", fields.join(","))?;
    }

    Ok(())
}

/// Write analysis results to a Parquet file
pub fn write_parquet(result: &AnalysisResult, filename: &str) -> Result<()> {
    let schema = make_schema();

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(Float32Array::from(result.w.clone())),
            Arc::new(Float32Array::from(result.q2.clone())),
            Arc::new(Float32Array::from(result.xb.clone())),
            Arc::new(Float32Array::from(result.theta_star.clone())),
            Arc::new(Float32Array::from(result.phi_star.clone())),
            Arc::new(Float32Array::from(result.mm.clone())),
            Arc::new(Float32Array::from(result.mm2.clone())),
            Arc::new(Int32Array::from(result.e_sector.clone())),
            Arc::new(Int32Array::from(result.event_type.clone())),
            Arc::new(Float32Array::from(result.beam_energy.clone())),
            Arc::new(Float32Array::from(result.e_prime.clone())),
            Arc::new(Float32Array::from(result.w_thrown.clone())),
            Arc::new(Float32Array::from(result.q2_thrown.clone())),
            Arc::new(Float32Array::from(result.mm_thrown.clone())),
            Arc::new(Float32Array::from(result.mm2_thrown.clone())),
            Arc::new(Float32Array::from(result.e_dc_vx.clone())),
            Arc::new(Float32Array::from(result.e_dc_vy.clone())),
            Arc::new(Float32Array::from(result.e_dc_vz.clone())),
            Arc::new(Int32Array::from(result.num_pip.clone())),
            Arc::new(Int32Array::from(result.num_prot.clone())),
            Arc::new(Int32Array::from(result.num_pim.clone())),
            Arc::new(Int32Array::from(result.num_pos.clone())),
            Arc::new(Int32Array::from(result.num_neg.clone())),
            Arc::new(Int32Array::from(result.num_neutral.clone())),
            Arc::new(Int32Array::from(result.num_photons.clone())),
            Arc::new(Float32Array::from(result.pi0_mass.clone())),
            Arc::new(Float32Array::from(result.pi0_mass2.clone())),
            Arc::new(Int32Array::from(result.part_index.clone())),
            Arc::new(Int32Array::from(result.part_id.clone())),
            Arc::new(Int32Array::from(result.part_q.clone())),
            Arc::new(Float32Array::from(result.part_p.clone())),
            Arc::new(Float32Array::from(result.part_beta.clone())),
            Arc::new(Float32Array::from(result.part_theta.clone())),
            Arc::new(Float32Array::from(result.part_phi.clone())),
            Arc::new(Int32Array::from(result.part_sector.clone())),
            Arc::new(Int32Array::from(result.part_sc_sector.clone())),
            Arc::new(Int32Array::from(result.part_sc_pd.clone())),
            Arc::new(Float32Array::from(result.part_delta_t_p.clone())),
            Arc::new(Float32Array::from(result.part_delta_t_pip.clone())),
            Arc::new(Float32Array::from(result.part_delta_t_e.clone())),
            Arc::new(Float32Array::from(result.part_delta_t_k.clone())),
            Arc::new(Int32Array::from(result.part_nphe.clone())),
            Arc::new(Int32Array::from(result.part_cc_segm.clone())),
            Arc::new(Float32Array::from(result.part_etot.clone())),
            Arc::new(Float32Array::from(result.part_ec_ei.clone())),
            Arc::new(Float32Array::from(result.part_ec_eo.clone())),
            Arc::new(Float32Array::from(result.part_dc_xsc.clone())),
            Arc::new(Float32Array::from(result.part_dc_ysc.clone())),
            Arc::new(Float32Array::from(result.part_dc_zsc.clone())),
            Arc::new(Float32Array::from(result.part_edep.clone())),
            Arc::new(Float32Array::from(result.part_vx.clone())),
            Arc::new(Float32Array::from(result.part_vy.clone())),
            Arc::new(Float32Array::from(result.part_vz.clone())),
            Arc::new(BooleanArray::from(result.part_fid_chern.clone())),
            Arc::new(BooleanArray::from(result.part_elec_fid.clone())),
            Arc::new(BooleanArray::from(result.part_hadron_fid.clone())),
            Arc::new(BooleanArray::from(result.part_dt_p_pass.clone())),
            Arc::new(BooleanArray::from(result.part_dt_pip_pass.clone())),
            Arc::new(BooleanArray::from(result.part_is_pip.clone())),
            Arc::new(BooleanArray::from(result.part_is_prot.clone())),
            Arc::new(BooleanArray::from(result.part_is_pim.clone())),
            Arc::new(BooleanArray::from(result.part_is_electron.clone())),
        ],
    )?;

    let file = File::create(filename)?;
    let mut writer = ArrowWriter::try_new(file, schema, None)?;
    writer.write(&batch)?;
    writer.close()?;

    Ok(())
}
