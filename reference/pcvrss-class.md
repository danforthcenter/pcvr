# Class `pcvrss` for models specified in `pcvr`.

Models specified by
[growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
or [mvSS](https://danforthcenter.github.io/pcvr/reference/mvSS.md) are
represented by a `pcvrss` object, which contains the model type,
formulas, starting values or priors, the data for the model to use, and
the model backend to use.

## Details

See `methods(class = "pcvrss")` for an overview of available methods.

## Slots

- `formula`:

  The formula that will be used to fit the model.

- `prior`:

  Priors if the model is a Bayesian model (ie using the brms backend).

- `initfun`:

  Initialization function if the model is a Bayesian model.

- `df`:

  The data that will be used to fit the model.

- `family`:

  The model family, currently only used in the brms backend.

- `pcvrForm`:

  The formula that was specified in
  [growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  and used in other pcvr functions.

- `type`:

  The model backend.

- `model`:

  The name of the main growth formula.

- `call`:

  The call to
  [growthSS](https://danforthcenter.github.io/pcvr/reference/growthSS.md)
  or [mvSS](https://danforthcenter.github.io/pcvr/reference/mvSS.md).

- `start`:

  Starting values for frequentist models.

- `taus`:

  Quantiles for nlrq/rq models.

## See also

[`growthSS`](https://danforthcenter.github.io/pcvr/reference/growthSS.md),
[`mvSS`](https://danforthcenter.github.io/pcvr/reference/mvSS.md)
