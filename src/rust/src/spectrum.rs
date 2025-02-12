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

impl TryFrom<&Robj> for Spectrum {
    type Error = Error;

    fn try_from(value: &Robj) -> Result<Self> {
        if let Some(class) = value.class() {
            let class = class.collect::<String>();
            match class.as_str() {
                "Spectrum" => (),
                _ => return Err(Error::from(format!("Expected Spectrum, got {:?}", class))),
            }
        } else {
            return Err(Error::from(format!("Expected Spectrum, got {:?}", value)));
        }
        let ptr: ExternalPtr<Spectrum> = value.try_into()?;

        Ok(ptr.as_ref().clone())
    }
}

impl Spectrum {
    pub(crate) fn recover_list(spectra: &List) -> Result<Vec<Spectrum>> {
        spectra
            .to_vec()
            .iter()
            .map(|r_obj| r_obj.try_into())
            .collect::<Result<Vec<Spectrum>>>()
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
