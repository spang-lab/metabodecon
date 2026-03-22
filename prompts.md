Analyze the function the `lorentz_sup()`. It takes a vector of x positions
(usually of length ~130k) plus vectors of parameters of lorentzian functions.

Then it iterates over each position, evaluates each lorentz function at that
position, and then sums up the results.

This is quite slow.

Please think about ways to speed this up. I have three strategies in mind:

1.  See if we can simplify the sum of lorentzians analytically to a single
    function and calculate that. This might just not be possible.

2.  Build the sum of lorentzians as R expression at runtime, using `substitute`
    and similar functions for creating expression/functions at runtime. E.g.
    for a list of three lorentzians, we could build an expression like this:

    abs(A[1] * (lambda[1] / (lambda[1]^2 + (x - x0[1])^2))) +
    abs(A[2] * (lambda[2] / (lambda[2]^2 + (x - x0[2])^2))) +
    abs(A[3] * (lambda[3] / (lambda[3]^2 + (x - x0[3])^2)))

    This expression can then be applied in one go to the whole `x` vector, making
    use of R's vectorized operators. This should be much faster than using sapply
    over each position in x. However, I'm not sure how big R-expression are
    allowed to be. Maybe there is a limit?

3.  Write a C function that iterates over lorentzians and datapoints and
    calculates the matrix directly. The function should be as simple as possible
    and have as few dependencies as possible. The number of operations should be
    kept as small as possible, e.g. by calculating lambda[j]^2 only once
    and not within each inner loop. Also, the C code should operate on the R
    objects directly, so I don't need to copy the inputs and/or outputs
    of the function in RAM.

Evaluate all three strategies, their advantages and weaknesses and then give
me a summary of what you think is feasible and what's the best strategy.