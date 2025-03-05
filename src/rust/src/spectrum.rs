use extendr_api::prelude::*;
use metabodecon::spectrum;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub(crate) struct Spectrum {
    inner: spectrum::Spectrum,
}

impl AsRef<spectrum::Spectrum> for Spectrum {
    fn as_ref(&self) -> &spectrum::Spectrum {
        &self.inner
    }
}

impl From<spectrum::Spectrum> for Spectrum {
    fn from(value: spectrum::Spectrum) -> Self {
        Self { inner: value }
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
            Ok(spectrum) => spectrum.into(),
            Err(error) => throw_r_error(format!("{}", error)),
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

    pub(crate) fn nucleus(&self) -> String {
        self.inner.nucleus().to_string()
    }

    pub(crate) fn frequency(&self) -> f64 {
        self.inner.frequency()
    }

    pub(crate) fn reference_compound(&self) -> Result<List> {
        let reference = self.inner.reference_compound();
        let chemical_shift = reference.chemical_shift();
        let index = reference.index();
        let name = reference.name();
        let referencing_method = reference
            .referencing_method()
            .map(|method| method.to_string());
        let mut result = HashMap::<&str, Robj>::new();
        result.insert("chemical_shift", chemical_shift.into());
        result.insert("index", index.into());
        result.insert("name", Nullable::from(name).into());
        result.insert(
            "referencing_method",
            Nullable::from(referencing_method).into(),
        );

        List::from_hashmap(result)
    }

    pub(crate) fn set_signal_boundaries(&mut self, signal_boundaries: Vec<f64>) {
        if signal_boundaries.len() != 2 {
            throw_r_error("signal_boundaries must be a vector of length 2");
        }
        let signal_boundaries = (signal_boundaries[0], signal_boundaries[1]);

        match self.inner.set_signal_boundaries(signal_boundaries) {
            Ok(_) => (),
            Err(error) => throw_r_error(format!("{}", error)),
        }
    }

    pub(crate) fn set_nucleus(&mut self, nucleus: &str) {
        self.inner
            .set_nucleus(std::str::FromStr::from_str(nucleus).unwrap());
    }

    pub(crate) fn set_frequency(&mut self, frequency: f64) {
        self.inner.set_frequency(frequency)
    }

    pub(crate) fn set_reference_compound(&mut self, reference: List) {
        let reference = reference.into_hashmap();
        let chemical_shift = reference
            .get("chemical_shift")
            .unwrap_or_else(|| throw_r_error("missing chemical_shift"))
            .as_real()
            .unwrap_or_else(|| throw_r_error("chemical_shift must be a numeric"));
        let index = reference
            .get("index")
            .unwrap_or_else(|| throw_r_error("missing index"))
            .as_integer()
            .unwrap_or_else(|| throw_r_error("index must be an integer"));
        let name = reference
            .get("name")
            .map(|name| name.as_str().unwrap().to_string());
        let referencing_method = reference
            .get("referencing_method")
            .map(|method| std::str::FromStr::from_str(method.as_str().unwrap()).ok())
            .flatten();
        let reference = spectrum::meta::ReferenceCompound::new(
            chemical_shift,
            index as usize,
            name,
            referencing_method,
        );
        self.inner.set_reference_compound(reference);
    }

    pub(crate) fn read_bruker(
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
            Ok(spectrum) => spectrum.into(),
            Err(error) => throw_r_error(format!("{}", error)),
        }
    }

    pub(crate) fn read_bruker_set(
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
                Ok(spectra) => spectra
                    .into_iter()
                    .map(|spectrum| spectrum.into())
                    .collect::<Vec<Spectrum>>(),
                Err(error) => throw_r_error(format!("{}", error)),
            };

        List::from_values(spectra)
    }

    pub(crate) fn read_jcampdx(path: &str, signal_boundaries: Vec<f64>) -> Self {
        if signal_boundaries.len() != 2 {
            throw_r_error("signal_boundaries must be a vector of length 2");
        }
        let signal_boundaries = (signal_boundaries[0], signal_boundaries[1]);

        match spectrum::JcampDx::read_spectrum(path, signal_boundaries) {
            Ok(spectrum) => spectrum.into(),
            Err(error) => throw_r_error(format!("{}", error)),
        }
    }

    pub(crate) fn write_json(&self, path: &str) {
        let serialized = match serde_json::to_string_pretty(self.as_ref()) {
            Ok(serialized) => serialized,
            Err(error) => throw_r_error(format!("{}", error)),
        };
        std::fs::write(path, serialized).unwrap();
    }

    pub(crate) fn read_json(path: &str) -> Self {
        let serialized = std::fs::read_to_string(path).unwrap();

        match serde_json::from_str::<spectrum::Spectrum>(&serialized) {
            Ok(deserialized) => deserialized.into(),
            Err(error) => throw_r_error(format!("{}", error)),
        }
    }

    pub(crate) fn write_bin(&self, path: &str) {
        let serialized = match rmp_serde::to_vec(self.as_ref()) {
            Ok(serialized) => serialized,
            Err(error) => throw_r_error(format!("{}", error)),
        };
        std::fs::write(path, serialized).unwrap();
    }

    pub(crate) fn read_bin(path: &str) -> Self {
        let serialized = std::fs::read(path).unwrap();

        match rmp_serde::from_slice::<spectrum::Spectrum>(&serialized) {
            Ok(deserialized) => deserialized.into(),
            Err(error) => throw_r_error(format!("{}", error)),
        }
    }
}

extendr_module! {
    mod spectrum;
    impl Spectrum;
}
