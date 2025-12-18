# possible that ss is not a pcvrss object for compatibility with other brms models if "df" is part of it then work with that otherwise use general data. NOTE ggplot2:::Stat\$compute_layer can use the default from ggproto \`make plot within a given panel of the ggplot (a facet)\` this is mostly the same as the default ggproto compute_panel function, but it takes more named args and passes them to compute_group and avoids warning about individual/time columns. \`make data out of model per a given aes-group, should only be 1 per panel\` this is the heavily customized component which makes data for ribbons from the model and ss objects.

possible that ss is not a pcvrss object for compatibility with other
brms models if "df" is part of it then work with that otherwise use
general data. NOTE ggplot2:::Stat\$compute_layer can use the default
from ggproto \`make plot within a given panel of the ggplot (a facet)\`
this is mostly the same as the default ggproto compute_panel function,
but it takes more named args and passes them to compute_group and avoids
warning about individual/time columns. \`make data out of model per a
given aes-group, should only be 1 per panel\` this is the heavily
customized component which makes data for ribbons from the model and ss
objects.

## Usage

``` r
statBrmsMod
```

## Format

An object of class `StatBrm` (inherits from `Stat`, `ggproto`, `gg`) of
length 6.
