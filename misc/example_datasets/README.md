* `blood`: Sixteen 1D CPMG NMR spectra of blood plasma in Bruker format
* `urine`: Two NOESY NMR spectra of urine.
* `test`: Same as Blood dataset. Renamed to `blood` in 2023/12/04 to improve
  docs ("blood" is more consistent with "urine" and avoids confusion with
  test-datasets and internal test files for the package). The old folder is kept
  around for backwards compatibility (because some people might use the old
  vignette that downloads and refers to test).
* `aki`: 106 Bruker urine spectra from MetaboLights study `MTBLS24` (72
  controls, 34 AKI cases, all sampled at 24 h). Full study download:
  <https://www.ebi.ac.uk/metabolights/MTBLS24>. The copy included here is a
  convenience subset containing `s_MTBLS24.txt` plus only the Bruker files
  needed for reading spectra (`1r`, `procs`, `acqus`); other spectra-related
  outputs and early analysis/classifier result files from Zacharias et al.
  (2012) are excluded to keep size and download time low.