use extendr_api::prelude::*;
use metabodecon::deconvolution;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub(crate) struct Deconvolution {
    inner: deconvolution::Deconvolution,
}

impl AsRef<deconvolution::Deconvolution> for Deconvolution {
    fn as_ref(&self) -> &deconvolution::Deconvolution {
        &self.inner
    }
}

impl From<deconvolution::Deconvolution> for Deconvolution {
    fn from(value: deconvolution::Deconvolution) -> Self {
        Self { inner: value }
    }
}

#[extendr]
impl Deconvolution {
    pub(crate) fn lorentzians(&self) -> Result<List> {
        let len = self.inner.lorentzians().len();
        let mut sf = Vec::<f64>::with_capacity(len);
        let mut hw = Vec::<f64>::with_capacity(len);
        let mut maxp = Vec::<f64>::with_capacity(len);
        self.inner.lorentzians().iter().for_each(|lorentzian| {
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

    pub(crate) fn mse(&self) -> f64 {
        self.inner.mse()
    }

    pub(crate) fn superposition(&self, chemical_shift: f64) -> f64 {
        deconvolution::Lorentzian::superposition(chemical_shift, self.inner.lorentzians())
    }

    pub(crate) fn superposition_vec(&self, chemical_shifts: Vec<f64>) -> Vec<f64> {
        deconvolution::Lorentzian::superposition_vec(&*chemical_shifts, self.inner.lorentzians())
    }

    pub(crate) fn par_superposition_vec(&self, chemical_shifts: Vec<f64>) -> Vec<f64> {
        deconvolution::Lorentzian::par_superposition_vec(
            &*chemical_shifts,
            self.inner.lorentzians(),
        )
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

        match serde_json::from_str::<deconvolution::Deconvolution>(&serialized) {
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

        match rmp_serde::from_slice::<deconvolution::Deconvolution>(&serialized) {
            Ok(deserialized) => deserialized.into(),
            Err(error) => throw_r_error(format!("{}", error)),
        }
    }
}

extendr_module! {
    mod deconvolution;
    impl Deconvolution;
}
