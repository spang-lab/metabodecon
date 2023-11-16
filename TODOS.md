# Todos

## Feature: Batch Mode

We should have a batch mode, that does all the above steps truly automatically and creates a pdf containing all quality control images. The pdf can be inspected later on and based on the findings the function call can be adjusted. Batch mode can also run in parallel to speed up calculations. Instead of waiting 1h we need to wait 3 or 6 minutes then.

## Refactor: Text Output

The output should be improved. The License should not be printed after every function execution, unless there is a strong reason to do so. Timestamps should be added to the output, so the user automatically has a rough idea how long the function will take to finish.

## Refactor: Plotting defaults

Function `plot_triplets` should not store to file by default. If we offer an option to store to a file for convenience, it shouldn't be png and we should print the file path.

## Refactor: Plotting speed

Function `plot_lorentz_curves_save_as_png` is suuuuper slow. We should try to make this quicker.
