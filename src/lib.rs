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
use analysis::reaction::{Reaction, MCReaction};
use output::parquet_writer::{write_csv, write_parquet, AnalysisResult};
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

/// Process a single event and return an AnalysisResult entry
fn process_event(
    event: &reader::event::Event,
    experiment: &str,
    beam_energy: f32,
    mc: bool,
) -> Option<(Reaction, Option<MCReaction>)> {
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

    // Create reaction and apply hadron cuts
    let mut reaction = if mc {
        MCReaction::new(event, beam_energy).base
    } else {
        Reaction::new(event, beam_energy)
    };

    // Create cuts once per event, not per particle
    let cuts = Cuts::new(event);

    for part_num in 1..event.gpart as usize {
        if part_num >= event.p.len() {
            break;
        }

        if cuts.pip(part_num) {
            reaction.set_pip(event.px(part_num), event.py(part_num), event.pz(part_num));
        } else if cuts.prot(part_num) {
            reaction.set_proton(event.px(part_num), event.py(part_num), event.pz(part_num));
        } else if cuts.pim(part_num) {
            reaction.set_pim(event.px(part_num), event.py(part_num), event.pz(part_num));
        } else {
            reaction.set_other(
                event.id(part_num),
                event.px(part_num),
                event.py(part_num),
                event.pz(part_num),
            );
        }
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

    Some((reaction, mc_react))
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
        if let Some((mut reaction, mc_react)) = process_event(event, experiment, beam_energy, mc) {
            let sector = event.dc_sect(0);
            let event_type = reaction.event_type();
            let theta_star = reaction.theta_star();
            let phi_star = reaction.phi_star();

            if mc {
                let mc_r = mc_react.unwrap();
                result.push_mc(
                    reaction.w(),
                    reaction.q2(),
                    reaction.xb(),
                    theta_star,
                    phi_star,
                    reaction.mm(),
                    reaction.mm2(),
                    sector,
                    event_type,
                    beam_energy,
                    reaction.e_prime(),
                    mc_r.w_thrown(),
                    mc_r.q2_thrown(),
                    mc_r.mm_thrown(),
                    mc_r.mm2_thrown(),
                );
            } else {
                result.push(
                    reaction.w(),
                    reaction.q2(),
                    reaction.xb(),
                    theta_star,
                    phi_star,
                    reaction.mm(),
                    reaction.mm2(),
                    sector,
                    event_type,
                    beam_energy,
                    reaction.e_prime(),
                );
            }
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
    info!("Processed {} events, {} passed cuts", total_events, result.n_events);

    // Convert to Python dict
    let dict = pyo3::types::PyDict::new_bound(py);
    dict.set_item("n_events", result.n_events)?;
    dict.set_item("w", &result.w)?;
    dict.set_item("q2", &result.q2)?;
    dict.set_item("xb", &result.xb)?;
    dict.set_item("theta_star", &result.theta_star)?;
    dict.set_item("phi_star", &result.phi_star)?;
    dict.set_item("mm", &result.mm)?;
    dict.set_item("mm2", &result.mm2)?;
    dict.set_item("sector", &result.sector)?;
    dict.set_item("event_type", &result.event_type)?;
    dict.set_item("beam_energy", &result.beam_energy)?;
    dict.set_item("e_prime", &result.e_prime)?;
    dict.set_item("w_thrown", &result.w_thrown)?;
    dict.set_item("q2_thrown", &result.q2_thrown)?;
    dict.set_item("mm_thrown", &result.mm_thrown)?;
    dict.set_item("mm2_thrown", &result.mm2_thrown)?;

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
                    info!("{}: {} events, {} passed cuts on {:?}", filename, n_events, result.n_events, std::thread::current().id());
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
    info!("Total events passing cuts: {}", result.n_events);

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
