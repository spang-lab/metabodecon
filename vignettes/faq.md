# FAQ

## Water Artifact

__Question:__ When I choose `n` for `Water artefact fully inside red vertical lines? (y/n)`, MetaboDecon1D asks me to `Choose another half width range (in ppm) for the water artefact`, instead of left and right borders. Why is that and what should I enter here?

__Answer:__ the water signal is always centered by Bruker, so you only have to specify the half width range to make the water area around the midpoint wider or narrower.

## Parameter Optimization

__Question:__ Why do you do exactly 10 iterations of parameter optimization?

__Answer:__ This value was determined empirically. We found that 10 iterations are enough to get a good fit and more iterations would not improve the fit significantly, but would take longer. Obviously, a more objective stopping criterion would be better and will likely be implemented in a future version.
