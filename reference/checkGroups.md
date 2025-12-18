# Helper function to check groups in data.

Helper function to check groups in data.

## Usage

``` r
checkGroups(df, group)
```

## Arguments

- df:

  Data frame to use.

- group:

  Set of variables to use in grouping observations. These taken together
  should identify a unique plant (or unique plant at a unique angle)
  across time.

## Value

If there are duplicates in the grouping then this will return a message
with code to start checking the duplicates in your data.

## Examples

``` r
df <- growthSim("linear",
  n = 10, t = 10,
  params = list("A" = c(2, 1.5))
)
checkGroups(df, c("time", "id", "group"))
#> Grouping is unique
df$time[12] <- 3
checkGroups(df, c("time", "id", "group"))
#> There are 1 observations that are not uniquely identified.
#> The max number of duplicates is 2.
#> Run `df[duplicated(interaction(df$time, df$id, df$group)),]` to see the duplicated rows,
#>  or df[interaction(df$time, df$id, df$group)=='3.id_2.a',] to see the first duplicated instance.
```
