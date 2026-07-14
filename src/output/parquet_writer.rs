use anyhow::Result;
use std::fs::File;
use std::io::{BufWriter, Write};

/// Analysis result for a single event
#[derive(Debug, Clone)]
pub struct AnalysisResult {
    pub w: Vec<f32>,
    pub q2: Vec<f32>,
    pub xb: Vec<f32>,
    pub theta_star: Vec<f32>,
    pub phi_star: Vec<f32>,
    pub mm: Vec<f32>,
    pub mm2: Vec<f32>,
    pub sector: Vec<i32>,
    pub event_type: Vec<i32>,
    pub beam_energy: Vec<f32>,
    pub e_prime: Vec<f32>,
    // MC fields (NaN if not MC)
    pub w_thrown: Vec<f32>,
    pub q2_thrown: Vec<f32>,
    pub mm_thrown: Vec<f32>,
    pub mm2_thrown: Vec<f32>,
    pub n_events: usize,
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
            sector: Vec::new(),
            event_type: Vec::new(),
            beam_energy: Vec::new(),
            e_prime: Vec::new(),
            w_thrown: Vec::new(),
            q2_thrown: Vec::new(),
            mm_thrown: Vec::new(),
            mm2_thrown: Vec::new(),
            n_events: 0,
        }
    }

    pub fn push(
        &mut self,
        w: f32,
        q2: f32,
        xb: f32,
        theta_star: f32,
        phi_star: f32,
        mm: f32,
        mm2: f32,
        sector: i32,
        event_type: i32,
        beam_energy: f32,
        e_prime: f32,
    ) {
        self.w.push(w);
        self.q2.push(q2);
        self.xb.push(xb);
        self.theta_star.push(theta_star);
        self.phi_star.push(phi_star);
        self.mm.push(mm);
        self.mm2.push(mm2);
        self.sector.push(sector);
        self.event_type.push(event_type);
        self.beam_energy.push(beam_energy);
        self.e_prime.push(e_prime);
        self.w_thrown.push(f32::NAN);
        self.q2_thrown.push(f32::NAN);
        self.mm_thrown.push(f32::NAN);
        self.mm2_thrown.push(f32::NAN);
        self.n_events += 1;
    }

    pub fn push_mc(
        &mut self,
        w: f32,
        q2: f32,
        xb: f32,
        theta_star: f32,
        phi_star: f32,
        mm: f32,
        mm2: f32,
        sector: i32,
        event_type: i32,
        beam_energy: f32,
        e_prime: f32,
        w_thrown: f32,
        q2_thrown: f32,
        mm_thrown: f32,
        mm2_thrown: f32,
    ) {
        self.w.push(w);
        self.q2.push(q2);
        self.xb.push(xb);
        self.theta_star.push(theta_star);
        self.phi_star.push(phi_star);
        self.mm.push(mm);
        self.mm2.push(mm2);
        self.sector.push(sector);
        self.event_type.push(event_type);
        self.beam_energy.push(beam_energy);
        self.e_prime.push(e_prime);
        self.w_thrown.push(w_thrown);
        self.q2_thrown.push(q2_thrown);
        self.mm_thrown.push(mm_thrown);
        self.mm2_thrown.push(mm2_thrown);
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
        self.sector.extend(other.sector);
        self.event_type.extend(other.event_type);
        self.beam_energy.extend(other.beam_energy);
        self.e_prime.extend(other.e_prime);
        self.w_thrown.extend(other.w_thrown);
        self.q2_thrown.extend(other.q2_thrown);
        self.mm_thrown.extend(other.mm_thrown);
        self.mm2_thrown.extend(other.mm2_thrown);
        self.n_events += other.n_events;
    }
}

/// Write analysis results to a CSV file
pub fn write_csv(result: &AnalysisResult, filename: &str) -> Result<()> {
    let file = File::create(filename)?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "w,q2,xb,theta_star,phi_star,mm,mm2,sector,event_type,beam_energy,e_prime,w_thrown,q2_thrown,mm_thrown,mm2_thrown"
    )?;

    for i in 0..result.n_events {
        writeln!(
            writer,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            result.w[i],
            result.q2[i],
            result.xb[i],
            result.theta_star[i],
            result.phi_star[i],
            result.mm[i],
            result.mm2[i],
            result.sector[i],
            result.event_type[i],
            result.beam_energy[i],
            result.e_prime[i],
            result.w_thrown[i],
            result.q2_thrown[i],
            result.mm_thrown[i],
            result.mm2_thrown[i],
        )?;
    }

    Ok(())
}
