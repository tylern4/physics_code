mod analysis;
mod output;
mod physics;
mod reader;

use pyo3::prelude::*;
use indicatif::{ProgressBar, ProgressStyle};
use rayon::prelude::*;
use tracing::{info, warn, debug, trace};
use tracing_subscriber::EnvFilter;
use std::sync::OnceLock;

use analysis::cuts::{Cuts, E1dCuts, E1fCuts, E16Cuts};
use analysis::delta_t;
use analysis::reaction::{Reaction, MCReaction};
use output::parquet_writer::{write_csv, write_parquet, AnalysisResult, ParticleRow};
use physics::constants;
use physics::four_momentum::FourMomentum;
use physics::kinematics;
use reader::root_reader;

static THREAD_POOL: OnceLock<rayon::ThreadPool> = OnceLock::new();

fn get_thread_pool() -> &'static rayon::ThreadPool {
    THREAD_POOL.get_or_init(|| {
        let num_threads = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(4);
        info!("Creating rayon thread pool with {} worker threads", num_threads);
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build()
            .unwrap()
    })
}

/// Initialize tracing/logging with env-filter support
fn init_logging() {
    let _ = tracing_subscriber::fmt()
        .with_env_filter(EnvFilter::from_default_env())
        .with_target(false)
        .with_thread_ids(true)
        .try_init();
}

/// Compute per-particle row data
fn make_particle_row(
    event: &reader::event::Event,
    part_num: usize,
    cuts: &Cuts,
    dt: &delta_t::DeltaT,
    is_electron: bool,
    is_pip: bool,
    is_prot: bool,
    is_pim: bool,
    fid_chern: bool,
    elec_fid: bool,
    hadron_fid: bool,
) -> ParticleRow {
    let p = event.p(part_num);
    let beta = event.b(part_num);
    let theta = kinematics::theta_calc(event.cz(part_num));
    let phi = kinematics::phi_calc(event.cx(part_num), event.cy(part_num));
    let sector = event.dc_sect(part_num);
    let sc_sector = event.sc_sect(part_num);

    ParticleRow {
        part_index: part_num as i32,
        part_id: event.id(part_num),
        part_q: event.q(part_num),
        part_p: p,
        part_beta: beta,
        part_theta: theta,
        part_phi: phi,
        part_sector: sector,
        part_sc_sector: sc_sector,
        part_sc_pd: event.sc_pd(part_num),
        part_delta_t_p: dt.get_dt_p(part_num),
        part_delta_t_pip: dt.get_dt_pi(part_num),
        part_delta_t_e: dt.get_dt_e(part_num),
        part_delta_t_k: dt.get_dt_k(part_num),
        part_nphe: event.nphe(part_num),
        part_cc_segm: event.cc_segm(part_num),
        part_etot: event.etot(part_num),
        part_ec_ei: event.ec_ei(part_num),
        part_ec_eo: event.ec_eo(part_num),
        part_dc_xsc: event.dc_xsc(part_num),
        part_dc_ysc: event.dc_ysc(part_num),
        part_dc_zsc: event.dc_zsc(part_num),
        part_edep: event.edep(part_num),
        part_vx: event.vx(part_num),
        part_vy: event.vy(part_num),
        part_vz: event.vz(part_num),
        part_fid_chern: fid_chern,
        part_elec_fid: elec_fid,
        part_hadron_fid: hadron_fid,
        part_dt_p_pass: cuts.dt_p_cut(part_num),
        part_dt_pip_pass: cuts.dt_pip_cut(part_num),
        part_is_pip: is_pip,
        part_is_prot: is_prot,
        part_is_pim: is_pim,
        part_is_electron: is_electron,
    }
}

/// Process a single event and return an AnalysisResult entry
fn process_event(
    event: &reader::event::Event,
    experiment: &str,
    beam_energy: f32,
    mc: bool,
) -> Option<(Reaction, Option<MCReaction>, Vec<ParticleRow>)> {
    if event.gpart < 1 {
        return None;
    }

    // Apply electron cuts based on experiment
    let electron_pass = match experiment {
        "e1d" => {
            let cuts = E1dCuts::new(event);
            cuts.is_electron()
        }
        "e1f" => {
            let cuts = E1fCuts::new(event);
            cuts.is_electron()
        }
        "e16" => {
            let cuts = E16Cuts::new(event);
            cuts.is_electron()
        }
        _ => {
            let cuts = Cuts::new(event);
            cuts.is_electron()
        }
    };

    if !electron_pass {
        return None;
    }

    // Compute electron-level quantities for cuts
    let cuts_base = Cuts::new(event);
    let fid_chern_elec = cuts_base.fid_chern_cut();
    let elec_fid = cuts_base.elec_fid_cut();

    // Create reaction and apply hadron cuts
    let mut reaction = if mc {
        MCReaction::new(event, beam_energy).base
    } else {
        Reaction::new(event, beam_energy)
    };

    // Build per-particle rows
    let mut particles = Vec::new();

    // Particle 0 is always the electron
    particles.push(make_particle_row(
        event, 0, &cuts_base, &cuts_base.dt,
        true, false, false, false,
        fid_chern_elec, elec_fid, false,
    ));

    for part_num in 1..event.gpart as usize {
        if part_num >= event.p.len() {
            break;
        }

        let is_pip = cuts_base.pip(part_num);
        let is_prot = cuts_base.prot(part_num);
        let is_pim = cuts_base.pim(part_num);
        let hadron_fid = cuts_base.hadron_fid_arjun(part_num);

        if is_pip {
            reaction.set_pip(event.px(part_num), event.py(part_num), event.pz(part_num));
        } else if is_prot {
            reaction.set_proton(event.px(part_num), event.py(part_num), event.pz(part_num));
        } else if is_pim {
            reaction.set_pim(event.px(part_num), event.py(part_num), event.pz(part_num));
        } else {
            reaction.set_other(
                event.id(part_num),
                event.px(part_num),
                event.py(part_num),
                event.pz(part_num),
            );
        }

        particles.push(make_particle_row(
            event, part_num, &cuts_base, &cuts_base.dt,
            false, is_pip, is_prot, is_pim,
            false, false, hadron_fid,
        ));
    }

    // Only keep events that pass the channel selection
    if !reaction.channel() {
        return None;
    }

    let mc_react = if mc {
        Some(MCReaction::new(event, beam_energy))
    } else {
        None
    };

    Some((reaction, mc_react, particles))
}

/// Process a chunk of events and return an AnalysisResult
fn process_chunk(
    events: &[reader::event::Event],
    experiment: &str,
    beam_energy: f32,
    mc: bool,
) -> AnalysisResult {
    trace!("processing chunk of {} events on {:?}", events.len(), std::thread::current().id());
    let mut result = AnalysisResult::new();
    for event in events {
        if let Some((mut reaction, mc_react, particles)) = process_event(event, experiment, beam_energy, mc) {
            let e_sector = event.dc_sect(0);
            let event_type = reaction.event_type();
            let theta_star = reaction.theta_star();
            let phi_star = reaction.phi_star();

            let (w_thrown, q2_thrown, mm_thrown, mm2_thrown) = if mc {
                let mc_r = mc_react.unwrap();
                (mc_r.w_thrown(), mc_r.q2_thrown(), mc_r.mm_thrown(), mc_r.mm2_thrown())
            } else {
                (f32::NAN, f32::NAN, f32::NAN, f32::NAN)
            };

            result.push_event(
                reaction.w(),
                reaction.q2(),
                reaction.xb(),
                theta_star,
                phi_star,
                reaction.mm(),
                reaction.mm2(),
                e_sector,
                event_type,
                beam_energy,
                reaction.e_prime(),
                w_thrown,
                q2_thrown,
                mm_thrown,
                mm2_thrown,
                event.dc_vx(0),
                event.dc_vy(0),
                event.dc_vz(0),
                reaction.num_pip,
                reaction.num_prot,
                reaction.num_pim,
                reaction.num_pos,
                reaction.num_neg,
                reaction.num_neutral,
                reaction.num_photons,
                reaction.pi0_mass(),
                reaction.pi0_mass2(),
                &particles,
            );
        }
    }
    result
}

/// Process a single ROOT file with cuts and return analysis results
#[pyfunction]
fn process_file(
    py: Python,
    filename: String,
    experiment: String,
    beam_energy: f32,
    mc: bool,
) -> PyResult<PyObject> {
    init_logging();
    info!("Processing file: {}", filename);

    let events = root_reader::read_root_file(&filename)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    let total_events = events.len();
    info!("Read {} events from {}", total_events, filename);

    let pb = ProgressBar::new(total_events as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    let num_threads = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4);
    info!("Using {} threads", num_threads);

    let chunk_size = (total_events + num_threads - 1) / num_threads;
    let pool = get_thread_pool();
    let experiment = experiment.clone();

    let results: Vec<AnalysisResult> = pool.install(|| {
        events.par_chunks(chunk_size).enumerate().map(|(i, chunk)| {
            debug!("chunk {} ({} events) on {:?}", i, chunk.len(), std::thread::current().id());
            process_chunk(chunk, &experiment, beam_energy, mc)
        }).collect()
    });

    let mut result = AnalysisResult::new();
    for r in results {
        result.merge(r);
    }

    pb.finish_with_message("done");
    info!("Processed {} events, {} passed cuts ({} particles)",
          total_events, result.n_events, result.n_particles);

    // Convert to Python dict
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("n_events", result.n_events)?;
    dict.set_item("n_particles", result.n_particles)?;
    dict.set_item("w", &result.w)?;
    dict.set_item("q2", &result.q2)?;
    dict.set_item("xb", &result.xb)?;
    dict.set_item("theta_star", &result.theta_star)?;
    dict.set_item("phi_star", &result.phi_star)?;
    dict.set_item("mm", &result.mm)?;
    dict.set_item("mm2", &result.mm2)?;
    dict.set_item("e_sector", &result.e_sector)?;
    dict.set_item("event_type", &result.event_type)?;
    dict.set_item("beam_energy", &result.beam_energy)?;
    dict.set_item("e_prime", &result.e_prime)?;
    dict.set_item("w_thrown", &result.w_thrown)?;
    dict.set_item("q2_thrown", &result.q2_thrown)?;
    dict.set_item("mm_thrown", &result.mm_thrown)?;
    dict.set_item("mm2_thrown", &result.mm2_thrown)?;
    dict.set_item("pi0_mass", &result.pi0_mass)?;
    dict.set_item("pi0_mass2", &result.pi0_mass2)?;

    Ok(dict.into())
}

/// Process multiple ROOT files and write results to CSV
#[pyfunction]
#[pyo3(signature = (filenames, experiment, beam_energy, output, mc, num_threads=0, batch_size=16, output_format="parquet"))]
fn process_files_to_parquet(
    py: Python,
    filenames: Vec<String>,
    experiment: String,
    beam_energy: f32,
    output: String,
    mc: bool,
    num_threads: usize,
    batch_size: usize,
    output_format: &str,
) -> PyResult<usize> {
    init_logging();
    info!("Processing {} files with {} cuts (E_beam = {} GeV){}",
          filenames.len(), experiment, beam_energy, if mc { " [MC]" } else { "" });

    let total_files = filenames.len();
    let pb = ProgressBar::new(total_files as u64);
    pb.set_style(
        ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} files ({eta})")
            .unwrap()
            .progress_chars("#>-"),
    );

    // Release the GIL for parallel processing
    let result = py.allow_threads(|| {
        let effective_threads = if num_threads > 0 { num_threads } else {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4)
        };
        info!("Creating rayon thread pool with {} worker threads, batch_size={}", effective_threads, batch_size);
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(effective_threads)
            .build()
            .unwrap();

        let experiment_clone = experiment.clone();

        let mut final_result = AnalysisResult::new();

        for batch in filenames.chunks(batch_size) {
            let batch_results: Vec<AnalysisResult> = pool.install(|| {
                batch.par_iter().filter_map(|filename| {
                    let events = match root_reader::read_root_file(filename) {
                        Ok(events) => events,
                        Err(e) => {
                            warn!("Failed to read {}: {}", filename, e);
                            pb.inc(1);
                            return None;
                        }
                    };
                    let n_events = events.len();
                    debug!("Read {} events from {} on {:?}", n_events, filename, std::thread::current().id());
                    pb.inc(1);

                    let result = process_chunk(&events, &experiment_clone, beam_energy, mc);
                    info!("{}: {} events, {} passed cuts ({} particles) on {:?}",
                          filename, n_events, result.n_events, result.n_particles, std::thread::current().id());
                    Some(result)
                }).collect()
            });

            for r in batch_results {
                final_result.merge(r);
            }
        }

        final_result
    });

    pb.finish_with_message("done");
    info!("Total events passing cuts: {}, total particles: {}", result.n_events, result.n_particles);

    match output_format {
        "csv" => write_csv(&result, &output)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?,
        _ => write_parquet(&result, &output)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?,
    }

    Ok(result.n_events)
}

/// Read a ROOT file and return the number of events
#[pyfunction]
fn read_root_file_py(filename: &str) -> PyResult<usize> {
    init_logging();
    info!("Reading file: {}", filename);
    let events = root_reader::read_root_file(filename)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    info!("Read {} events", events.len());
    Ok(events.len())
}

/// Read a ROOT file and return event data as a dictionary
#[pyfunction]
fn read_root_file_summary(
    py: Python,
    filename: &str,
) -> PyResult<Vec<PyObject>> {
    init_logging();
    info!("Reading summary from: {}", filename);
    let events = root_reader::read_root_file(filename)
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

    let mut results = Vec::new();
    for event in events.iter().take(10) {
        let dict = pyo3::types::PyDict::new_bound(py);
        dict.set_item("gpart", event.gpart)?;
        dict.set_item("npart", event.npart)?;
        dict.set_item("evntid", event.evntid)?;
        results.push(dict.into());
    }
    Ok(results)
}

// Physics functions exposed to Python
#[pyfunction]
fn q2_calc(
    e_beam_px: f32,
    e_beam_py: f32,
    e_beam_pz: f32,
    e_scattered_px: f32,
    e_scattered_py: f32,
    e_scattered_pz: f32,
) -> f32 {
    let beam = FourMomentum::new(e_beam_px, e_beam_py, e_beam_pz, constants::MASS_E);
    let scattered = FourMomentum::new(
        e_scattered_px,
        e_scattered_py,
        e_scattered_pz,
        constants::MASS_E,
    );
    kinematics::q2_calc(&beam, &scattered)
}

#[pyfunction]
fn w_calc(
    e_beam_px: f32,
    e_beam_py: f32,
    e_beam_pz: f32,
    e_scattered_px: f32,
    e_scattered_py: f32,
    e_scattered_pz: f32,
) -> f32 {
    let beam = FourMomentum::new(e_beam_px, e_beam_py, e_beam_pz, constants::MASS_E);
    let scattered = FourMomentum::new(
        e_scattered_px,
        e_scattered_py,
        e_scattered_pz,
        constants::MASS_E,
    );
    kinematics::w_calc(&beam, &scattered)
}

#[pyfunction]
fn get_sector(phi: f32) -> i32 {
    kinematics::get_sector(phi)
}

#[pyfunction]
fn get_mass(pid: i32) -> f32 {
    constants::get_mass(pid)
}

#[pyfunction]
fn theta_calc(cz: f32) -> f32 {
    kinematics::theta_calc(cz)
}

#[pyfunction]
fn phi_calc(cx: f32, cy: f32) -> f32 {
    kinematics::phi_calc(cx, cy)
}

/// Python module definition
#[pymodule]
fn _lib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(read_root_file_py, m)?)?;
    m.add_function(wrap_pyfunction!(read_root_file_summary, m)?)?;
    m.add_function(wrap_pyfunction!(process_file, m)?)?;
    m.add_function(wrap_pyfunction!(process_files_to_parquet, m)?)?;
    m.add_function(wrap_pyfunction!(q2_calc, m)?)?;
    m.add_function(wrap_pyfunction!(w_calc, m)?)?;
    m.add_function(wrap_pyfunction!(get_sector, m)?)?;
    m.add_function(wrap_pyfunction!(get_mass, m)?)?;
    m.add_function(wrap_pyfunction!(theta_calc, m)?)?;
    m.add_function(wrap_pyfunction!(phi_calc, m)?)?;

    Ok(())
}
