# The SAP Dataset

The SAP Dataset consists of a single 'Simple-As-Possible' (SAP)
spectrum. The purpose of the SAP spectrum is to provide a
straightforward example that can be used to test and understand the
deconvolution algorithm in detail.

## Usage

``` r
sap
```

## Format

An object of class `spectra` of length 1.

## Details

The first (and only) spectrum within the SAP dataset contains 128
datapoints ranging from -6.4 to 6.4 ppm with four peaks. A rough sketch
of the spectrum is shown below:

    -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~
    |      SFR      |               w               |     SFR      |
    |               |  x           www       p      |              |
    |~-~-~-~-~-~-~-~|~-|-|-~-~-~-~-~|~-~-~-~-|-~-~-~-~-~-~-~-~-~-~-~
    |               |  | |          |        |      |
    6.4             |  | 2.24       0.047    -2.22  -3.2
                    |  2.61
                    3.2
