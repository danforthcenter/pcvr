# Class `conjugate` for output from the `pcvr::conjugate` function.

Comparisons made by the `conjugate` function return objects of this
class containing parameters of the prior and posterior distributions,
hypothesis tests, ROPE tests, Bayes Factors, and plots of the posterior.

## Details

See `methods(class = "conjugate")` for an overview of available methods.

## Slots

- `summary`:

  Summary data frame of results

- `posterior`:

  Posterior distribution as a list of named lists

- `prior`:

  Prior distribution as a list of named lists

- `plot`:

  Optionally a plot of the distributions and their differences

- `data`:

  The data from s1 and s2 arguments to
  [conjugate](https://danforthcenter.github.io/pcvr/reference/conjugate.md).

- `call`:

  Matched call to
  [conjugate](https://danforthcenter.github.io/pcvr/reference/conjugate.md).

## See also

[`conjugate`](https://danforthcenter.github.io/pcvr/reference/conjugate.md)
