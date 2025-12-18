# Read in plantCV csv output in wide or long format

Read in plantCV csv output in wide or long format

## Usage

``` r
read.pcv(
  filepath,
  mode = NULL,
  traitCol = "trait",
  labelCol = "label",
  valueCol = "value",
  reader = NULL,
  filters = NULL,
  awk = NULL,
  ...
)
```

## Arguments

- filepath:

  Path to csv file of plantCV output.

- mode:

  NULL (the default) or one of "wide" or "long", partial string matching
  is supported. This controls whether data is **returned** in long or
  wide format. If left NULL then the output format will be the same as
  the input format.

- traitCol:

  Column with phenotype names, defaults to "trait". This should
  generally not need to be changed from the default. This, labelCol, and
  valueCol are used to determine if data are in long format in their raw
  state (the csv file itself).

- labelCol:

  Column with phenotype labels (units), defaults to "label". This should
  generally not need to be changed from the default. This is used with
  traitCol when `mode="wide"` to identify unique traits since some may
  be ambiguous (ellipseCenter.x vs ellipseCenter.y, bins of histograms,
  etc)

- valueCol:

  Column with phenotype values, defaults to "value". This should
  generally not need to be changed from the default.

- reader:

  The function to use to read in data, defaults to NULL in which case
  [`data.table::fread`](https://rdatatable.gitlab.io/data.table/reference/fread.html)
  is used if filters are in place and `read.csv` is used otherwise. Note
  that if you use `read.csv` with filters in place then you will need to
  specify `header=FALSE` so that the piped output from awk is read
  correctly. If fread is too slow for your needs then `vroom::vroom()`
  may be useful.

- filters:

  If a very large pcv output file is read then it may be desireable to
  subset it before reading it into R, either for ease of use or because
  of RAM limitations. The filter argument works with "COLUMN in VALUES"
  syntax. This can either be a character vector or a list of character
  vectors. In these vectors there needs to be a column name, one of " in
  ", " is ", or " = " to match the string exactly, or "contains" to
  match with awk style regex, then a set of comma delimited values to
  filter that column for (see examples). Note that this and awk both use
  awk through [`pipe()`](https://rdrr.io/r/base/connections.html). This
  functionality will not work on a windows system.

- awk:

  As an alternative to filters a direct call to awk can be supplied
  here, in which case that call will be used through
  [`pipe()`](https://rdrr.io/r/base/connections.html).

- ...:

  Other arguments passed to the reader function. In the case of 'fread'
  there are several defaults provided already which can be overwritten
  with these extra arguments.

## Value

Returns a data.frame in wide or long format.

## Details

In plantCV version 4 the single value traits are returned in wide format
from `json2csv` and the multi value traits are returned in long format.
Briefly plantCV data was returned as one long table which sparked the
emphasis in this function on reading data quickly and parsing it outside
of R. With the current plantCV output these options are largely
unnecessary. When data is read in using read.pcv the traitCol, valueCol,
and labelCol arguments are checked to determine if the data is in long
format. This is done to keep compatibility with interim versions of
plantcv output where all outputs were in a single long format file.

With the current implementation and plantcv output you can read wide or
long format files into wide or long format in R. Keep in mind that the
'mode' argument controls the format that will be returned in R, not the
format that the data saved as in your csv file.

## Examples

``` r
# \donttest{
tryCatch(
  {
    mv <- paste0(
      "https://media.githubusercontent.com/media/joshqsumner/",
      "pcvrTestData/main/pcv4-multi-value-traits.csv"
    )
    sv <- paste0(
      "https://raw.githubusercontent.com/joshqsumner/",
      "pcvrTestData/main/pcv4-single-value-traits.csv"
    )

    w2w <- read.pcv(sv, mode = "wide", reader = "fread")
    dim(w2w)

    w2l <- read.pcv(sv, mode = "long", reader = "fread")
    dim(w2l)

    l2w <- read.pcv(mv, mode = "wide", reader = "fread")
    dim(l2w)

    l2l <- read.pcv(mv, mode = "long", reader = "fread")
    dim(l2l)
  },
  error = function(e) {
    message(e)
  }
)
#> [1] 513720      6
# }
```
