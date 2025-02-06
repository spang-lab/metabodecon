use extendr_api::prelude::*;

mod deconvoluter;
mod deconvolution;
mod spectrum;

extendr_module! {
    mod metabodecon;
    use deconvoluter;
    use deconvolution;
    use spectrum;
}
