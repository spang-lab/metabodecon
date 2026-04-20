
# Review and fix R CMD check findings

1. The follwowing check fail on Windows-latest:

   > ══ Failed tests ════════════════════════════════════════════════════════════════
   > ── Failure ('test-datadir.R:36:5'): datadir works if datadir_persistent=filled ──
   > Expected `x` to equal `y`.
   > Differences:
   > actual vs expected
   > - "C:/Users/runneradmin/AppData/Local/Temp/RtmpqkDWz3/working_dir/RtmpYjzTAM/metabodecon/data"
   > + "C:/Users/runneradmin/AppData/Local/Temp/RtmpqkDWz3/working_dir/RtmpYjzTAM/metabodecon/mocks/datadir/persistent/filled"

   Analyze how mocking of the datadir works inside evalwidth and then see if
   you can understand why the test fails. If you cannot find the cause easily
   within one minute, ignore this error. I will analyze it on a Windows machine
   later on.

2. Vignette building fails for R versions < 4.2:

   > -- R CMD build -----------------------------------------------------------------
   > * checking for file 'D:\a\metabodecon\metabodecon/DESCRIPTION' ... OK
   > * preparing 'metabodecon':
   > * checking DESCRIPTION meta-information ... OK
   > * installing the package to process help pages
   > * saving partial Rd database
   > * creating vignettes ... ERROR
   > Error: --- re-building 'Get_Started.Rmd' using rmarkdown
   > --- finished re-building 'Get_Started.Rmd'
   > --- re-building 'MDM.Rmd' using rmarkdown
   > Quitting from MDM.Rmd:124-131 [fit-default]
   > ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   > <error/rlang_error>
   > Error:
   > ! Using the Rust backend requires mdrb >= 0.0.1.
   > To install or upgrade mdrb run: install_mdrb()
   > To check system requirements run: check_mdrb_deps()
   > For more information see: https://github.com/spang-lab/mdrb
   > ---
   > Backtrace:
   >     x
   >  1. \-metabodecon::fit_mdm(...)
   >  2.   \-metabodecon::deconvolute(...) at metabodecon/R/mdm.R:275:5
   >  3.     \-metabodecon::check_mdrb(stop_on_fail = TRUE) at metabodecon/R/decon.R:113:5
   > ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   > Error: Error: processing vignette 'MDM.Rmd' failed with diagnostics:
   > Using the Rust backend requires mdrb >= 0.0.1.
   > To install or upgrade mdrb run: install_mdrb()
   > To check system requirements run: check_mdrb_deps()
   > For more information see: https://github.com/spang-lab/mdrb
   > --- failed re-building 'MDM.Rmd'
   > SUMMARY: processing the following file failed:
   >   'MDM.Rmd'
   > Error: Error: Vignette re-building failed.
   > Execution halted
   > Error: Error in proc$get_built_file() : Build process failed
   > Calls: <Anonymous> ... build_package -> with_envvar -> force -> <Anonymous>
   > Execution halted
   > Error: Process completed with exit code 1.
   > Run ## --------------------------------------------------------------------
   > Show testthat output
   > Run actions/upload-artifact@v4
   > Warning: No files were found with the provided path: D:/a/metabodecon/metabodecon/check. No artifacts will be uploaded.

   To fix this we should replace `use_rust=TRUE` with `use_rust=0.5`, which will
   use the newest R implementation, which is almost as fast as the Rust
   implementation.

3. The following tests are currently skipped during CI runs:

   > ══ Skipped tests (13) ══════════════════════════════════════════════════════════
   > • Manual Test for r4 only (1): 'test-align-sim.R:3:1'
   > • Manual execution only (1): 'test-benchmark_mdlm.R:3:1'
   > • On CI (4): 'test-draw_spectrum.R:2:1', 'test-plot_sfr.R:2:1',
   >   'test-plot_spectrum.R:2:1', 'test-plot_ws.R:2:1'
   > • benchmark_extension_table not yet implemented (2): 'test-mdm.R:56:3',
   >   'test-mdm.R:79:3'
   > • find_best_params not yet implemented (1): 'test-find_best_params.R:2:5'
   > • mdm(X, y, model=...) API not yet implemented (4): 'test-mdm.R:5:3',
   >   'test-mdm.R:25:3', 'test-mdm.R:42:3', 'test-mdm.R:131:3'

   This is fine:

   - On CI (4): 'test-draw_spectrum.R:2:1', 'test-plot_sfr.R:2:1',
        'test-plot_spectrum.R:2:1', 'test-plot_ws.R:2:1'

   These should probably be deleted (unless they really need to be checked in the
   future again, in which case we must not skip them during automated tests, which
   means we need to ensure they are fast)

   - benchmark_extension_table not yet implemented (2): 'test-mdm.R:56:3',
     'test-mdm.R:79:3'
   - find_best_params not yet implemented (1): 'test-find_best_params.R:2:5'
   - Manual Test for r4 only (1): 'test-align-sim.R:3:1'
   - Manual execution only (1): 'test-benchmark_mdlm.R:3:1'

   These should probably be activated:

   - mdm(X, y, model=...) API not yet implemented (4): 'test-mdm.R:5:3',
     'test-mdm.R:25:3', 'test-mdm.R:42:3', 'test-mdm.R:131:3'


# Remove backwards compatibility code

- Remove support for `ask`
- Remove support for `wshw`
