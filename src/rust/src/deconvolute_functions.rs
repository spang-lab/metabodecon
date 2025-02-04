use extendr_api::prelude::*;
use metabodecon::{deconvolution, spectrum};
use std::collections::HashMap;

fn convert_to_spectrum(spectrum: List) -> spectrum::Spectrum {
    let spectrum = spectrum.into_hashmap();
    let chemical_shifts = match spectrum.get("cs") {
        Some(chemical_shifts) => chemical_shifts.as_real_vector().unwrap(),
        None => throw_r_error("cs attribute not found"),
    };
    let intensities = match spectrum.get("si") {
        Some(intensities) => intensities.as_real_vector().unwrap(),
        None => throw_r_error("intensities attribute not found"),
    };
    let signal_boundaries = match spectrum.get("sfr") {
        Some(signal_boundaries) => {
            if signal_boundaries.len() != 2 {
                throw_r_error("signal_boundaries must be a vector of length 2");
            }
            let tmp: Vec<f64> = signal_boundaries.as_real_vector().unwrap();

            (tmp[0], tmp[1])
        }
        None => throw_r_error("signal_boundaries attribute not found"),
    };
    let spectrum = match spectrum::Spectrum::new(chemical_shifts, intensities, signal_boundaries) {
        Ok(rust_spectrum) => rust_spectrum,
        Err(e) => throw_r_error(format!("{}", e)),
    };

    spectrum
}

fn convert_from_deconvolution(deconvolution: deconvolution::Deconvolution) -> Result<List> {
    let len = deconvolution.lorentzians().len();
    let mut sf = Vec::<f64>::with_capacity(len);
    let mut hw = Vec::<f64>::with_capacity(len);
    let mut maxp = Vec::<f64>::with_capacity(len);
    deconvolution.lorentzians().iter().for_each(|lorentzian| {
        sf.push(lorentzian.sf());
        hw.push(lorentzian.hw());
        maxp.push(lorentzian.maxp());
    });

    let mut result = HashMap::<&str, Robj>::new();
    result.insert("A", sf.into());
    result.insert("lambda", hw.into());
    result.insert("x0", maxp.into());

    List::from_hashmap(result)
}

fn set_up_deconvoluter(ignore_regions: Vec<f64>) -> deconvolution::Deconvoluter {
    if ignore_regions.len() % 2 != 0 {
        throw_r_error("ignore_regions must be a list of pairs");
    }
    let mut deconvoluter = deconvolution::Deconvoluter::default();
    let ignore_regions = ignore_regions
        .chunks(2)
        .map(|region| (region[0], region[1]))
        .collect::<Vec<(f64, f64)>>();
    ignore_regions.iter().for_each(|region| {
        match deconvoluter.add_ignore_region((region.0, region.1)) {
            Ok(_) => (),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    });

    deconvoluter
}

/// Deconvolutes a spectrum.
/// @export
#[extendr]
fn backend_deconvolute_spectrum(spectrum: List, ignore_regions: Vec<f64>) -> Result<List> {
    let spectrum = convert_to_spectrum(spectrum);
    let deconvoluter = set_up_deconvoluter(ignore_regions);
    let deconvolution = match deconvoluter.deconvolute_spectrum(&spectrum) {
        Ok(deconvolution) => deconvolution,
        Err(e) => throw_r_error(format!("{}", e)),
    };

    convert_from_deconvolution(deconvolution)
}

/// Deconvolutes a spectrum in parallel.
/// @export
#[extendr]
fn backend_par_deconvolute_spectrum(spectrum: List, ignore_regions: Vec<f64>) -> Result<List> {
    let spectrum = convert_to_spectrum(spectrum);
    let deconvoluter = set_up_deconvoluter(ignore_regions);
    let deconvolution = match deconvoluter.par_deconvolute_spectrum(&spectrum) {
        Ok(deconvolution) => deconvolution,
        Err(e) => throw_r_error(format!("{}", e)),
    };

    convert_from_deconvolution(deconvolution)
}

extendr_module! {
    mod deconvolute_functions;
    fn backend_deconvolute_spectrum;
    fn backend_par_deconvolute_spectrum;
}
