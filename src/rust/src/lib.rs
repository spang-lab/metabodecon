use extendr_api::prelude::*;

mod deconvoluter;
mod deconvolution;
mod lorentzian;
mod spectrum;

extendr_module! {
    mod metabodecon;
    use deconvoluter;
    use deconvolution;
    use lorentzian;
    use spectrum;
}
