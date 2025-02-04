use crate::deconvolution::Deconvolution;
use crate::spectrum::Spectrum;
use extendr_api::prelude::*;
use metabodecon::deconvolution;

#[derive(Clone, Debug, Default)]
pub(crate) struct Deconvoluter {
    inner: deconvolution::Deconvoluter,
}

#[extendr]
impl Deconvoluter {
    pub(crate) fn new() -> Self {
        Self::default()
    }

    pub(crate) fn set_moving_average_smoother(&mut self, iterations: usize, window_size: usize) {
        match self
            .inner
            .set_smoothing_settings(deconvolution::SmoothingSettings::MovingAverage {
                iterations,
                window_size,
            }) {
            Ok(_) => (),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn set_noise_score_selector(&mut self, threshold: f64) {
        match self.inner.set_selection_settings(
            deconvolution::SelectionSettings::NoiseScoreFilter {
                scoring_method: deconvolution::ScoringMethod::MinimumSum,
                threshold,
            },
        ) {
            Ok(_) => (),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn set_analytical_fitter(&mut self, iterations: usize) {
        match self
            .inner
            .set_fitting_settings(deconvolution::FittingSettings::Analytical { iterations })
        {
            Ok(_) => (),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn add_ignore_region(&mut self, start: f64, end: f64) {
        match self.inner.add_ignore_region((start, end)) {
            Ok(_) => (),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn clear_ignore_regions(&mut self) {
        self.inner.clear_ignore_regions();
    }

    pub(crate) fn deconvolute_spectrum(&self, spectrum: &Spectrum) -> Deconvolution {
        match self.inner.deconvolute_spectrum(spectrum.as_ref()) {
            Ok(deconvolution) => deconvolution.into(),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn par_deconvolute_spectrum(&self, spectrum: &Spectrum) -> Deconvolution {
        match self.inner.par_deconvolute_spectrum(spectrum.as_ref()) {
            Ok(deconvolution) => deconvolution.into(),
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }
}

extendr_module! {
    mod deconvoluter;
    impl Deconvoluter;
}
