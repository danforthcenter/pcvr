# Calculate pseudo water use efficiency from phenotype and watering data

Rate based water use efficiency (WUE) is the change in biomass per unit
of water metabolized. Using image based phenotypes and watering data we
can calculate pseudo-WUE (pwue) over time. Here area_pixels is used as a
proxy for biomass and transpiration is approximated using watering data.
The equation is then \\ \frac{P\_{t} -
P\_{t-1}}{W\_{t\_{end-1}}-W\_{t\_{start}} }\\, where P is the phenotype
and W is the weight before watering.

Absolute value based WUE is the amount of water used to sustain a plants
biomass over a given period. The equation is then
\\\frac{P\_{t}}{W\_{t\_{end-1}}-W\_{t\_{start}} }\\

## Usage

``` r
pwue(
  df,
  w = NULL,
  pheno = "area_pixels",
  time = "timestamp",
  id = "barcode",
  offset = 0,
  pre_watering = "weight_before",
  post_watering = "weight_after",
  method = "rate"
)
```

## Arguments

- df:

  Dataframe containing wide single-value phenotype data. This should
  already be aggregated to one row per plant per day (angles/rotations
  combined).

- w:

  Watering data as returned from pcv.water.

- pheno:

  Phenotype column name, defaults to "area_pixels"

- time:

  Variable(s) that identify a plant on a given day. If variable names
  are different between df and w then this can be a vector of two names.

- id:

  Variable(s) that identify a plant over time. Defaults to `"barcode"`.
  If variable names are different between df and w then this can be a
  vector of two names.

- offset:

  Optionally you can specify how long before imaging a watering should
  not be taken into account. This defaults to 0, meaning that if a plant
  were watered directly before being imaged then that water would be
  counted towards WUE between the current image and the prior one. This
  argument is taken to be in seconds.

- pre_watering:

  Column containing weight before watering in `w`, defaults to
  "weight_before".

- post_watering:

  Column containing weight after watering in `w`, defaults to
  "weight_after".

- method:

  Which method to use, options are "rate", "abs", and "ndt". The "rate"
  method considers WUE as the change in a phenotype divided by the
  amount of water added. The "abs" method considers WUE as the amount of
  water used by a plant given its absolute size. The "ndt" method
  calculates normalized daily transpiration, which is the reciprocal of
  the "abs" method. The "rate" method is for questions more related to
  efficiency in using water to grow while "abs"/"ndt" are more suited to
  questions about how efficient a plant is at maintaining size given
  some amount of water or how much water it uses at a given size.

## Value

A data frame containing the watering data and to phenotype data with new
columns for change in the phenotype (`pheno_diff`), amount of water used
(`total_water`) over the interval between phenotype measurements (water
added post to pre phenotype measurement), `start` and `end` times for
the interval as well as their difference (`timeLengthSeconds`), and
pseudo water use efficiency (`pWUE`).

## Examples

``` r
set.seed(123)
weight_before <- sort(round(rnorm(20, 100, 10), 0))
weight_after <- sort(round(rnorm(20, 120, 10), 0))
df <- data.frame(
  time = seq_along(weight_before),
  area_pixels = round(130 / (1 + exp( (12 - seq_along(weight_before))/3) ), 0),
  weight_before,
  weight_after,
  barcode = 1,
  other = "x"
)
ex <- pwue(df, time = "time", method = "rate", id = c("barcode", "other"))
w <- df[, c("time", "weight_before", "weight_after", "barcode")]
ex2 <- pwue(df, w, time = c("time", "time"), method = "abs")
ex3 <- pwue(df, w, time = c("time", "time"), method = "ndt")
```
