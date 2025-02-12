use extendr_api::prelude::*;
use metabodecon::deconvolution;

#[derive(Copy, Clone, Debug)]
pub(crate) struct Lorentzian {
    inner: deconvolution::Lorentzian,
}

impl AsRef<deconvolution::Lorentzian> for Lorentzian {
    fn as_ref(&self) -> &deconvolution::Lorentzian {
        &self.inner
    }
}

impl From<deconvolution::Lorentzian> for Lorentzian {
    fn from(inner: deconvolution::Lorentzian) -> Self {
        Self { inner }
    }
}

impl Lorentzian {
    pub(crate) fn inner_from_parameters(
        sf: Vec<f64>,
        hw: Vec<f64>,
        maxp: Vec<f64>,
    ) -> Vec<deconvolution::Lorentzian> {
        if sf.len() != hw.len() || sf.len() != maxp.len() {
            throw_r_error("Length of sf, hw, and maxp must be equal.");
        }

        sf.iter()
            .zip(hw.iter())
            .zip(maxp.iter())
            .map(|((sf, hw), maxp)| deconvolution::Lorentzian::new(sf * hw, hw.powi(2), *maxp))
            .collect()
    }
}

#[extendr]
impl Lorentzian {
    pub(crate) fn new(sf: f64, hw: f64, maxp: f64) -> Self {
        Self {
            inner: deconvolution::Lorentzian::new(sf * hw, hw.powi(2), maxp),
        }
    }

    pub(crate) fn sf(&self) -> f64 {
        self.inner.sf()
    }

    pub(crate) fn hw(&self) -> f64 {
        self.inner.hw()
    }

    pub(crate) fn maxp(&self) -> f64 {
        self.inner.maxp()
    }

    pub(crate) fn set_sf(&mut self, sf: f64) {
        self.inner.set_sf(sf);
    }

    pub(crate) fn set_hw(&mut self, hw: f64) {
        self.inner.set_hw(hw);
    }

    pub(crate) fn set_maxp(&mut self, maxp: f64) {
        self.inner.set_maxp(maxp);
    }

    pub(crate) fn evaluate(&self, x: f64) -> f64 {
        self.inner.evaluate(x)
    }

    pub(crate) fn evaluate_vec(&self, x: Vec<f64>) -> Vec<f64> {
        self.inner.evaluate_vec(&x)
    }

    pub(crate) fn superposition(x: f64, sf: Vec<f64>, hw: Vec<f64>, maxp: Vec<f64>) -> f64 {
        let lorentzians = Self::inner_from_parameters(sf, hw, maxp);

        deconvolution::Lorentzian::superposition(x, &lorentzians)
    }

    pub(crate) fn superposition_vec(
        x: Vec<f64>,
        sf: Vec<f64>,
        hw: Vec<f64>,
        maxp: Vec<f64>,
    ) -> Vec<f64> {
        let lorentzians = Self::inner_from_parameters(sf, hw, maxp);

        deconvolution::Lorentzian::superposition_vec(&x, &lorentzians)
    }

    pub(crate) fn par_superposition_vec(
        x: Vec<f64>,
        sf: Vec<f64>,
        hw: Vec<f64>,
        maxp: Vec<f64>,
    ) -> Vec<f64> {
        let lorentzians = Self::inner_from_parameters(sf, hw, maxp);

        deconvolution::Lorentzian::par_superposition_vec(&x, &lorentzians)
    }
}

extendr_module! {
    mod lorentzian;
    impl Lorentzian;
}
