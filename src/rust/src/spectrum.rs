use extendr_api::prelude::*;
use metabodecon::spectrum;

#[derive(Clone, Debug)]
pub(crate) struct Spectrum {
    inner: spectrum::Spectrum,
}

impl AsRef<spectrum::Spectrum> for Spectrum {
    fn as_ref(&self) -> &spectrum::Spectrum {
        &self.inner
    }
}

#[extendr]
impl Spectrum {
    pub(crate) fn new(
        chemical_shifts: Vec<f64>,
        intensities: Vec<f64>,
        signal_boundaries: Vec<f64>,
    ) -> Self {
        if signal_boundaries.len() != 2 {
            throw_r_error("signal_boundaries must be a vector of length 2");
        }
        let signal_boundaries = (signal_boundaries[0], signal_boundaries[1]);

        match spectrum::Spectrum::new(chemical_shifts, intensities, signal_boundaries) {
            Ok(inner) => Self { inner },
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn from_bruker(
        path: &str,
        experiment: u32,
        processing: u32,
        signal_boundaries: Vec<f64>,
    ) -> Self {
        if signal_boundaries.len() != 2 {
            throw_r_error("signal_boundaries must be a vector of length 2");
        }
        let signal_boundaries = (signal_boundaries[0], signal_boundaries[1]);

        match spectrum::Bruker::read_spectrum(path, experiment, processing, signal_boundaries) {
            Ok(inner) => Self { inner },
            Err(e) => throw_r_error(format!("{}", e)),
        }
    }

    pub(crate) fn from_bruker_set(
        path: &str,
        experiment: u32,
        processing: u32,
        signal_boundaries: Vec<f64>,
    ) -> List {
        if signal_boundaries.len() != 2 {
            throw_r_error("signal_boundaries must be a vector of length 2");
        }
        let signal_boundaries = (signal_boundaries[0], signal_boundaries[1]);
        let spectra =
            match spectrum::Bruker::read_spectra(path, experiment, processing, signal_boundaries) {
                Ok(spectra) => spectra,
                Err(e) => throw_r_error(format!("{}", e)),
            };

        List::from_values(
            spectra
                .into_iter()
                .map(|spectrum| Robj::from(Spectrum { inner: spectrum })),
        )
    }

    pub(crate) fn chemical_shifts(&self) -> Vec<f64> {
        self.inner.chemical_shifts().to_vec()
    }

    pub(crate) fn intensities(&self) -> Vec<f64> {
        self.inner.intensities().to_vec()
    }

    pub(crate) fn signal_boundaries(&self) -> Vec<f64> {
        vec![
            self.inner.signal_boundaries().0,
            self.inner.signal_boundaries().1,
        ]
    }
}

extendr_module! {
    mod spectrum;
    impl Spectrum;
}
