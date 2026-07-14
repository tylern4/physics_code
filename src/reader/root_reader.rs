use anyhow::{Context, Result};
use oxyroot::{RootFile, Slice};

use super::event::Event;

/// Helper to read a vector branch as Slice<T> -> Vec<T>
fn read_vec_branch_i32(tree: &oxyroot::ReaderTree, name: &str) -> Vec<Vec<i32>> {
    tree.branch(name)
        .and_then(|b| b.as_iter::<Slice<i32>>().ok())
        .map(|iter| {
            iter.map(|s| s.into_vec()).collect()
        })
        .unwrap_or_default()
}

fn read_vec_branch_f32(tree: &oxyroot::ReaderTree, name: &str) -> Vec<Vec<f32>> {
    tree.branch(name)
        .and_then(|b| b.as_iter::<Slice<f32>>().ok())
        .map(|iter| {
            iter.map(|s| s.into_vec()).collect()
        })
        .unwrap_or_default()
}

/// Read a ROOT file and return a vector of events from the h10 TTree.
pub fn read_root_file(filename: &str) -> Result<Vec<Event>> {
    let mut file = RootFile::open(filename)
        .with_context(|| format!("Failed to open ROOT file: {}", filename))?;

    let tree_name = "h10";
    let tree = file
        .get_tree(tree_name)
        .with_context(|| format!("Failed to find tree '{}' in {}", tree_name, filename))?;

    // Read scalar branches (one value per entry)
    let gpart_vec: Vec<i32> = tree.branch("gpart")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();

    let n_entries = gpart_vec.len();
    if n_entries == 0 {
        return Ok(Vec::new());
    }

    let npart_vec: Vec<i32> = tree.branch("npart")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let evstat_vec: Vec<i32> = tree.branch("evstat")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let intt_vec: Vec<i32> = tree.branch("intt")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let evntid_vec: Vec<i32> = tree.branch("evntid")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let evtype_vec: Vec<i32> = tree.branch("evtype")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let evntclas_vec: Vec<i32> = tree.branch("evntclas")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let evthel_vec: Vec<i32> = tree.branch("evthel")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let evntclas2_vec: Vec<i32> = tree.branch("evntclas2")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let q_l_vec: Vec<f32> = tree.branch("q_l")
        .and_then(|b| b.as_iter::<f32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let t_l_vec: Vec<f32> = tree.branch("t_l")
        .and_then(|b| b.as_iter::<f32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let tr_time_vec: Vec<f32> = tree.branch("tr_time")
        .and_then(|b| b.as_iter::<f32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let rf_time1_vec: Vec<f32> = tree.branch("rf_time1")
        .and_then(|b| b.as_iter::<f32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let rf_time2_vec: Vec<f32> = tree.branch("rf_time2")
        .and_then(|b| b.as_iter::<f32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();

    // Read per-particle vector branches
    let id_vecs = read_vec_branch_i32(&tree, "id");
    let stat_vecs = read_vec_branch_i32(&tree, "stat");
    let dc_vecs = read_vec_branch_i32(&tree, "dc");
    let cc_vecs = read_vec_branch_i32(&tree, "cc");
    let sc_vecs = read_vec_branch_i32(&tree, "sc");
    let ec_vecs = read_vec_branch_i32(&tree, "ec");
    let lec_vecs = read_vec_branch_i32(&tree, "lec");
    let ccst_vecs = read_vec_branch_i32(&tree, "ccst");
    let p_vecs = read_vec_branch_f32(&tree, "p");
    let q_vecs = read_vec_branch_i32(&tree, "q");
    let b_vecs = read_vec_branch_f32(&tree, "b");
    let cx_vecs = read_vec_branch_f32(&tree, "cx");
    let cy_vecs = read_vec_branch_f32(&tree, "cy");
    let cz_vecs = read_vec_branch_f32(&tree, "cz");
    let vx_vecs = read_vec_branch_f32(&tree, "vx");
    let vy_vecs = read_vec_branch_f32(&tree, "vy");
    let vz_vecs = read_vec_branch_f32(&tree, "vz");

    // DC bank scalar
    let dc_part_vec: Vec<i32> = tree.branch("dc_part")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    // DC bank vectors
    let dc_sect_vecs = read_vec_branch_i32(&tree, "dc_sect");
    let dc_trk_vecs = read_vec_branch_i32(&tree, "dc_trk");
    let dc_stat_vecs = read_vec_branch_i32(&tree, "dc_stat");
    let dc_vx_vecs = read_vec_branch_f32(&tree, "dc_vx");
    let dc_vy_vecs = read_vec_branch_f32(&tree, "dc_vy");
    let dc_vz_vecs = read_vec_branch_f32(&tree, "dc_vz");
    let dc_vr_vecs = read_vec_branch_f32(&tree, "dc_vr");
    let dc_xsc_vecs = read_vec_branch_f32(&tree, "dc_xsc");
    let dc_ysc_vecs = read_vec_branch_f32(&tree, "dc_ysc");
    let dc_zsc_vecs = read_vec_branch_f32(&tree, "dc_zsc");
    let dc_cxsc_vecs = read_vec_branch_f32(&tree, "dc_cxsc");
    let dc_cysc_vecs = read_vec_branch_f32(&tree, "dc_cysc");
    let dc_czsc_vecs = read_vec_branch_f32(&tree, "dc_czsc");
    let dc_c2_vecs = read_vec_branch_f32(&tree, "dc_c2");

    // EC bank scalar
    let ec_part_vec: Vec<i32> = tree.branch("ec_part")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    // EC bank vectors
    let ec_stat_vecs = read_vec_branch_i32(&tree, "ec_stat");
    let ec_sect_vecs = read_vec_branch_i32(&tree, "ec_sect");
    let ec_whol_vecs = read_vec_branch_i32(&tree, "ec_whol");
    let ec_inst_vecs = read_vec_branch_i32(&tree, "ec_inst");
    let ec_oust_vecs = read_vec_branch_i32(&tree, "ec_oust");
    let etot_vecs = read_vec_branch_f32(&tree, "etot");
    let ec_ei_vecs = read_vec_branch_f32(&tree, "ec_ei");
    let ec_eo_vecs = read_vec_branch_f32(&tree, "ec_eo");
    let ec_t_vecs = read_vec_branch_f32(&tree, "ec_t");
    let ec_r_vecs = read_vec_branch_f32(&tree, "ec_r");
    let ech_x_vecs = read_vec_branch_f32(&tree, "ech_x");
    let ech_y_vecs = read_vec_branch_f32(&tree, "ech_y");
    let ech_z_vecs = read_vec_branch_f32(&tree, "ech_z");
    let ec_m2_vecs = read_vec_branch_f32(&tree, "ec_m2");
    let ec_m3_vecs = read_vec_branch_f32(&tree, "ec_m3");
    let ec_m4_vecs = read_vec_branch_f32(&tree, "ec_m4");
    let ec_c2_vecs = read_vec_branch_f32(&tree, "ec_c2");

    // SC bank scalar
    let sc_part_vec: Vec<i32> = tree.branch("sc_part")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    // SC bank vectors
    let sc_sect_vecs = read_vec_branch_i32(&tree, "sc_sect");
    let sc_hit_vecs = read_vec_branch_i32(&tree, "sc_hit");
    let sc_pd_vecs = read_vec_branch_i32(&tree, "sc_pd");
    let sc_stat_vecs = read_vec_branch_i32(&tree, "sc_stat");
    let edep_vecs = read_vec_branch_f32(&tree, "edep");
    let sc_t_vecs = read_vec_branch_f32(&tree, "sc_t");
    let sc_r_vecs = read_vec_branch_f32(&tree, "sc_r");
    let sc_c2_vecs = read_vec_branch_f32(&tree, "sc_c2");

    // CC bank scalar
    let cc_part_vec: Vec<i32> = tree.branch("cc_part")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    // CC bank vectors
    let cc_sect_vecs = read_vec_branch_i32(&tree, "cc_sect");
    let cc_hit_vecs = read_vec_branch_i32(&tree, "cc_hit");
    let cc_segm_vecs = read_vec_branch_i32(&tree, "cc_segm");
    let nphe_vecs = read_vec_branch_i32(&tree, "nphe");
    let cc_t_vecs = read_vec_branch_f32(&tree, "cc_t");
    let cc_r_vecs = read_vec_branch_f32(&tree, "cc_r");
    let cc_c2_vecs = read_vec_branch_f32(&tree, "cc_c2");

    // MC branches (may not exist for data)
    let nprt_vec: Vec<i32> = tree.branch("nprt")
        .and_then(|b| b.as_iter::<i32>().ok())
        .map(|iter| iter.collect())
        .unwrap_or_default();
    let pidpart_vecs = read_vec_branch_i32(&tree, "pidpart");
    let pxpart_vecs = read_vec_branch_f32(&tree, "pxpart");
    let pypart_vecs = read_vec_branch_f32(&tree, "pypart");
    let pzpart_vecs = read_vec_branch_f32(&tree, "pzpart");
    let epart_vecs = read_vec_branch_f32(&tree, "epart");
    let xpart_vecs = read_vec_branch_f32(&tree, "xpart");
    let ypart_vecs = read_vec_branch_f32(&tree, "ypart");
    let zpart_vecs = read_vec_branch_f32(&tree, "zpart");
    let qpart_vecs = read_vec_branch_f32(&tree, "qpart");
    let flagspart_vecs = read_vec_branch_i32(&tree, "flagspart");

    // Build events
    let mut events = Vec::with_capacity(n_entries);

    for i in 0..n_entries {
        let gpart = gpart_vec[i];
        let mut event = Event::new();

        event.npart = npart_vec.get(i).copied().unwrap_or(0);
        event.evstat = evstat_vec.get(i).copied().unwrap_or(0);
        event.intt = intt_vec.get(i).copied().unwrap_or(0);
        event.evntid = evntid_vec.get(i).copied().unwrap_or(0);
        event.evtype = evtype_vec.get(i).copied().unwrap_or(0);
        event.evntclas = evntclas_vec.get(i).copied().unwrap_or(0);
        event.evthel = evthel_vec.get(i).copied().unwrap_or(0);
        event.evntclas2 = evntclas2_vec.get(i).copied().unwrap_or(0);
        event.q_l = q_l_vec.get(i).copied().unwrap_or(0.0);
        event.t_l = t_l_vec.get(i).copied().unwrap_or(0.0);
        event.tr_time = tr_time_vec.get(i).copied().unwrap_or(0.0);
        event.rf_time1 = rf_time1_vec.get(i).copied().unwrap_or(0.0);
        event.rf_time2 = rf_time2_vec.get(i).copied().unwrap_or(0.0);
        event.gpart = gpart;
        event.dc_part = dc_part_vec.get(i).copied().unwrap_or(0);
        event.ec_part = ec_part_vec.get(i).copied().unwrap_or(0);
        event.sc_part = sc_part_vec.get(i).copied().unwrap_or(0);
        event.cc_part = cc_part_vec.get(i).copied().unwrap_or(0);
        event.nprt = nprt_vec.get(i).copied().unwrap_or(0);

        copy_vec(&id_vecs, i, &mut event.id);
        copy_vec(&stat_vecs, i, &mut event.stat);
        copy_vec(&dc_vecs, i, &mut event.dc);
        copy_vec(&cc_vecs, i, &mut event.cc);
        copy_vec(&sc_vecs, i, &mut event.sc);
        copy_vec(&ec_vecs, i, &mut event.ec);
        copy_vec(&lec_vecs, i, &mut event.lec);
        copy_vec(&ccst_vecs, i, &mut event.ccst);
        copy_vec_f(&p_vecs, i, &mut event.p);
        copy_vec(&q_vecs, i, &mut event.q);
        copy_vec_f(&b_vecs, i, &mut event.b);
        copy_vec_f(&cx_vecs, i, &mut event.cx);
        copy_vec_f(&cy_vecs, i, &mut event.cy);
        copy_vec_f(&cz_vecs, i, &mut event.cz);
        copy_vec_f(&vx_vecs, i, &mut event.vx);
        copy_vec_f(&vy_vecs, i, &mut event.vy);
        copy_vec_f(&vz_vecs, i, &mut event.vz);

        copy_vec(&dc_sect_vecs, i, &mut event.dc_sect);
        copy_vec(&dc_trk_vecs, i, &mut event.dc_trk);
        copy_vec(&dc_stat_vecs, i, &mut event.dc_stat);
        copy_vec_f(&dc_vx_vecs, i, &mut event.dc_vx);
        copy_vec_f(&dc_vy_vecs, i, &mut event.dc_vy);
        copy_vec_f(&dc_vz_vecs, i, &mut event.dc_vz);
        copy_vec_f(&dc_vr_vecs, i, &mut event.dc_vr);
        copy_vec_f(&dc_xsc_vecs, i, &mut event.dc_xsc);
        copy_vec_f(&dc_ysc_vecs, i, &mut event.dc_ysc);
        copy_vec_f(&dc_zsc_vecs, i, &mut event.dc_zsc);
        copy_vec_f(&dc_cxsc_vecs, i, &mut event.dc_cxsc);
        copy_vec_f(&dc_cysc_vecs, i, &mut event.dc_cysc);
        copy_vec_f(&dc_czsc_vecs, i, &mut event.dc_czsc);
        copy_vec_f(&dc_c2_vecs, i, &mut event.dc_c2);

        copy_vec(&ec_stat_vecs, i, &mut event.ec_stat);
        copy_vec(&ec_sect_vecs, i, &mut event.ec_sect);
        copy_vec(&ec_whol_vecs, i, &mut event.ec_whol);
        copy_vec(&ec_inst_vecs, i, &mut event.ec_inst);
        copy_vec(&ec_oust_vecs, i, &mut event.ec_oust);
        copy_vec_f(&etot_vecs, i, &mut event.etot);
        copy_vec_f(&ec_ei_vecs, i, &mut event.ec_ei);
        copy_vec_f(&ec_eo_vecs, i, &mut event.ec_eo);
        copy_vec_f(&ec_t_vecs, i, &mut event.ec_t);
        copy_vec_f(&ec_r_vecs, i, &mut event.ec_r);
        copy_vec_f(&ech_x_vecs, i, &mut event.ech_x);
        copy_vec_f(&ech_y_vecs, i, &mut event.ech_y);
        copy_vec_f(&ech_z_vecs, i, &mut event.ech_z);
        copy_vec_f(&ec_m2_vecs, i, &mut event.ec_m2);
        copy_vec_f(&ec_m3_vecs, i, &mut event.ec_m3);
        copy_vec_f(&ec_m4_vecs, i, &mut event.ec_m4);
        copy_vec_f(&ec_c2_vecs, i, &mut event.ec_c2);

        copy_vec(&sc_sect_vecs, i, &mut event.sc_sect);
        copy_vec(&sc_hit_vecs, i, &mut event.sc_hit);
        copy_vec(&sc_pd_vecs, i, &mut event.sc_pd);
        copy_vec(&sc_stat_vecs, i, &mut event.sc_stat);
        copy_vec_f(&edep_vecs, i, &mut event.edep);
        copy_vec_f(&sc_t_vecs, i, &mut event.sc_t);
        copy_vec_f(&sc_r_vecs, i, &mut event.sc_r);
        copy_vec_f(&sc_c2_vecs, i, &mut event.sc_c2);

        copy_vec(&cc_sect_vecs, i, &mut event.cc_sect);
        copy_vec(&cc_hit_vecs, i, &mut event.cc_hit);
        copy_vec(&cc_segm_vecs, i, &mut event.cc_segm);
        copy_vec(&nphe_vecs, i, &mut event.nphe);
        copy_vec_f(&cc_t_vecs, i, &mut event.cc_t);
        copy_vec_f(&cc_r_vecs, i, &mut event.cc_r);
        copy_vec_f(&cc_c2_vecs, i, &mut event.cc_c2);

        copy_vec(&pidpart_vecs, i, &mut event.pidpart);
        copy_vec_f(&xpart_vecs, i, &mut event.xpart);
        copy_vec_f(&ypart_vecs, i, &mut event.ypart);
        copy_vec_f(&zpart_vecs, i, &mut event.zpart);
        copy_vec_f(&epart_vecs, i, &mut event.epart);
        copy_vec_f(&pxpart_vecs, i, &mut event.pxpart);
        copy_vec_f(&pypart_vecs, i, &mut event.pypart);
        copy_vec_f(&pzpart_vecs, i, &mut event.pzpart);
        copy_vec_f(&qpart_vecs, i, &mut event.qpart);
        copy_vec(&flagspart_vecs, i, &mut event.flagspart);

        events.push(event);
    }

    Ok(events)
}

/// Copy from a Vec<Vec<i32>> at index into a fixed-size array
fn copy_vec(src: &Vec<Vec<i32>>, idx: usize, dst: &mut [i32]) {
    if let Some(row) = src.get(idx) {
        let n = row.len().min(dst.len());
        dst[..n].copy_from_slice(&row[..n]);
    }
}

fn copy_vec_f(src: &Vec<Vec<f32>>, idx: usize, dst: &mut [f32]) {
    if let Some(row) = src.get(idx) {
        let n = row.len().min(dst.len());
        dst[..n].copy_from_slice(&row[..n]);
    }
}
