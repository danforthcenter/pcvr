# Plot a `conjugate` object.

Plot a `conjugate` object.

## Usage

``` r
# S3 method for class 'conjugate'
plot(x, ...)
```

## Arguments

- x:

  An object of class `conjugate`.

- ...:

  further arguments, ignored.

## Examples

``` r
x <- conjugate(
  s1 = rnorm(10, 50, 10), s2 = rnorm(10, 60, 12), method = "t",
  priors = list(list(mu = 40, sd = 10), list(mu = 45, sd = 8)),
  rope_range = c(-5, 8), rope_ci = 0.89,
  cred.int.level = 0.89, hypothesis = "equal"
)
plot(x)

```
