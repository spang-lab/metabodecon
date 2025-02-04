use extendr_api::prelude::*;

mod deconvolute_functions;
mod deconvoluter;
mod deconvolution;
mod spectrum;

extendr_module! {
    mod metabodecon;
    use deconvolute_functions;
    use deconvoluter;
    use deconvolution;
    use spectrum;
}
