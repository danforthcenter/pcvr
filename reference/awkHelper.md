# subset helper function for use reading in large data, called in pcv.sub.read

subset helper function for use reading in large data, called in
pcv.sub.read

## Usage

``` r
awkHelper(inputFile, filters, awk = NULL)
```

## Arguments

- inputFile:

  Path to csv file of plantCV output, should be provided internally in
  read.pcv

- filters:

  filtering conditions, see read.pcv for details. Format as list("trait
  in area, perimeter", "other contains stringToMatch")

- awk:

  Optional awk command to use instead.

## Value

Returns a character string representing a unix style awk statement which
is typically passed to `pipe` or used as a connection in
[`data.table::fread`](https://rdatatable.gitlab.io/data.table/reference/fread.html).

## Details

awkHelper attempts to make awk commands from human readable input.
Currently when filters are supplied the input file has quotes removed by
\`sed\` then is piped into awk, so an equivalent command line statement
may be:
`sed 's/\"//g' pcvrTest2.csv | awk -F ',' '{ if (NR==1 || $18=="area") { print } }'`

## Examples

``` r
tryCatch(
  { # in case offline
    link1 <- "https://gist.githubusercontent.com/seankross/"
    link2 <- "a412dfbd88b3db70b74b/raw/5f23f993cd87c283ce766e7ac6b329ee7cc2e1d1/mtcars.csv"
    file <- paste0(link1, link2)
    awkHelper(file, list("gear in 4, 3"), awk = NULL)
    awkHelper(file, "gear contains 3", awk = NULL)
    # note that to be filtered the file has to exist on your local system, this example only shows
    # the output of awkHelper, which would then be executed by read.pcv on a unix system
    awkHelper(file, list("gear in 4, 3"), awk = "existing_command")
  },
  error = function(e) {
    message(e)
  }
)
#> [1] "existing_command"
```
