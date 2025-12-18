# Multi Value Trait Aggregation function

EMD can get very heavy with large datasets. For an example lemnatech
dataset filtering for images from every 5th day there are 6332^2 =
40,094,224 pairwise EMD values. In long format that's a 40 million row
dataframe, which is unwieldy. This function is to help reduce the size
of datasets before comparing histograms and moving on with matrix
methods or network analysis.

## Usage

``` r
mv_ag(
  df,
  group,
  mvCols = "frequencies",
  n_per_group = 1,
  outRows = NULL,
  keep = NULL,
  parallel = getOption("mc.cores", 1),
  traitCol = "trait",
  labelCol = "label",
  valueCol = "value",
  id = "image"
)
```

## Arguments

- df:

  A dataframe with multi value traits. This can be in wide or long
  format, data is assumed to be long if traitCol, valueCol, and labelCol
  are present.

- group:

  Vector of column names for variables which uniquely identify groups in
  the data to summarize data over. Typically this would be the design
  variables and a time variable.

- mvCols:

  Either a vector of column names/positions representing multi value
  traits or a character string that identifies the multi value trait
  columns as a regex pattern. Defaults to "frequencies".

- n_per_group:

  Number of rows to return for each group.

- outRows:

  Optionally this is a different way to specify how many rows to return.
  This will often not be exact so that groups have the same number of
  observations each.

- keep:

  A vector of single value traits to also average over groups, if there
  are a mix of single and multi value traits in your data.

- parallel:

  Optionally the groups can be run in parallel with this number of
  cores, defaults to 1 if the "mc.cores" option is not set globally.

- traitCol:

  Column with phenotype names, defaults to "trait".

- labelCol:

  Column with phenotype labels (units), defaults to "label".

- valueCol:

  Column with phenotype values, defaults to "value".

- id:

  Column that uniquely identifies images if the data is in long format.
  This is ignored when data is in wide format.

## Value

Returns a dataframe summarized by the specified groups over the
multi-value traits.

## Examples

``` r
s1 <- mvSim(
  dists = list(runif = list(min = 15, max = 150)),
  n_samples = 10,
  counts = 1000,
  min_bin = 1,
  max_bin = 180,
  wide = TRUE
)
mv_ag(s1, group = "group", mvCols = "sim_", n_per_group = 2)
#>             group sim_1 sim_2 sim_3 sim_4 sim_5 sim_6 sim_7 sim_8 sim_9 sim_10
#> runif_1.1 runif_1     0     0     0     0     0     0     0     0     0      0
#> runif_1.2 runif_1     0     0     0     0     0     0     0     0     0      0
#>           sim_11 sim_12 sim_13 sim_14 sim_15 sim_16 sim_17 sim_18 sim_19 sim_20
#> runif_1.1      0      0      0      0 0.0052 0.0096 0.0070  0.006 0.0052 0.0062
#> runif_1.2      0      0      0      0 0.0080 0.0068 0.0068  0.010 0.0062 0.0074
#>           sim_21 sim_22 sim_23 sim_24 sim_25 sim_26 sim_27 sim_28 sim_29 sim_30
#> runif_1.1 0.0052 0.0072 0.0068 0.0070 0.0074 0.0072 0.0066 0.0070 0.0084 0.0086
#> runif_1.2 0.0072 0.0066 0.0084 0.0094 0.0054 0.0068 0.0082 0.0082 0.0082 0.0088
#>           sim_31 sim_32 sim_33 sim_34 sim_35 sim_36 sim_37 sim_38 sim_39 sim_40
#> runif_1.1 0.0082 0.0072 0.0088 0.0078 0.0084 0.0058 0.0068 0.0084  0.007  0.006
#> runif_1.2 0.0090 0.0084 0.0064 0.0070 0.0076 0.0050 0.0082 0.0066  0.006  0.009
#>           sim_41 sim_42 sim_43 sim_44 sim_45 sim_46 sim_47 sim_48 sim_49 sim_50
#> runif_1.1  0.007 0.0094 0.0064 0.0066 0.0066 0.0064 0.0070 0.0086 0.0056 0.0086
#> runif_1.2  0.008 0.0082 0.0072 0.0068 0.0076 0.0064 0.0092 0.0074 0.0072 0.0062
#>           sim_51 sim_52 sim_53 sim_54 sim_55 sim_56 sim_57 sim_58 sim_59 sim_60
#> runif_1.1  0.007 0.0074 0.0084 0.0078 0.0084 0.0058 0.0072 0.0084 0.0074 0.0058
#> runif_1.2  0.006 0.0088 0.0086 0.0084 0.0068 0.0072 0.0064 0.0082 0.0080 0.0096
#>           sim_61 sim_62 sim_63 sim_64 sim_65 sim_66 sim_67 sim_68 sim_69 sim_70
#> runif_1.1 0.0060 0.0072 0.0058 0.0082 0.0086 0.0072 0.0076 0.0084 0.0080 0.0080
#> runif_1.2 0.0076 0.0084 0.0068 0.0076 0.0082 0.0072 0.0078 0.0066 0.0074 0.0072
#>           sim_71 sim_72 sim_73 sim_74 sim_75 sim_76 sim_77 sim_78 sim_79 sim_80
#> runif_1.1 0.0082 0.0062 0.0078 0.0068 0.0058 0.0096 0.0070 0.0066 0.0084 0.0062
#> runif_1.2 0.0064 0.0086 0.0086 0.0076 0.0080 0.0094 0.0094 0.0076 0.0072 0.0052
#>           sim_81 sim_82 sim_83 sim_84 sim_85 sim_86 sim_87 sim_88 sim_89 sim_90
#> runif_1.1 0.0078 0.0082 0.0072 0.0078 0.0072 0.0084 0.0074 0.0078 0.0054 0.0072
#> runif_1.2 0.0074 0.0080 0.0064 0.0074 0.0050 0.0064 0.0066 0.0060 0.0058 0.0072
#>           sim_91 sim_92 sim_93 sim_94 sim_95 sim_96 sim_97 sim_98 sim_99
#> runif_1.1 0.0066 0.0090 0.0062 0.0072 0.0062 0.0086  0.008 0.0076 0.0056
#> runif_1.2 0.0068 0.0078 0.0066 0.0066 0.0088 0.0102  0.008 0.0064 0.0066
#>           sim_100 sim_101 sim_102 sim_103 sim_104 sim_105 sim_106 sim_107
#> runif_1.1  0.0072  0.0060  0.0090  0.0052  0.0080  0.0082  0.0072  0.0074
#> runif_1.2  0.0070  0.0058  0.0066  0.0068  0.0076  0.0086  0.0046  0.0076
#>           sim_108 sim_109 sim_110 sim_111 sim_112 sim_113 sim_114 sim_115
#> runif_1.1  0.0080  0.0084  0.0062   0.006  0.0054  0.0098  0.0084  0.0074
#> runif_1.2  0.0082  0.0082  0.0092   0.006  0.0054  0.0070  0.0084  0.0066
#>           sim_116 sim_117 sim_118 sim_119 sim_120 sim_121 sim_122 sim_123
#> runif_1.1  0.0076  0.0086  0.0084  0.0082  0.0074  0.0072  0.0054   0.008
#> runif_1.2  0.0056  0.0052  0.0056  0.0072  0.0086  0.0084  0.0080   0.009
#>           sim_124 sim_125 sim_126 sim_127 sim_128 sim_129 sim_130 sim_131
#> runif_1.1  0.0074  0.0072  0.0068  0.0088  0.0070  0.0078  0.0082  0.0076
#> runif_1.2  0.0064  0.0090  0.0070  0.0082  0.0074  0.0066  0.0058  0.0070
#>           sim_132 sim_133 sim_134 sim_135 sim_136 sim_137 sim_138 sim_139
#> runif_1.1  0.0072  0.0066  0.0084  0.0084  0.0078  0.0062  0.0068  0.0080
#> runif_1.2  0.0088  0.0064  0.0100  0.0060  0.0092  0.0060  0.0094  0.0074
#>           sim_140 sim_141 sim_142 sim_143 sim_144 sim_145 sim_146 sim_147
#> runif_1.1  0.0086   0.008  0.0110  0.0094  0.0068  0.0086  0.0058  0.0086
#> runif_1.2  0.0072   0.006  0.0092  0.0072  0.0078  0.0074  0.0072  0.0080
#>           sim_148 sim_149 sim_150 sim_151 sim_152 sim_153 sim_154 sim_155
#> runif_1.1  0.0086  0.0080       0       0       0       0       0       0
#> runif_1.2  0.0052  0.0084       0       0       0       0       0       0
#>           sim_156 sim_157 sim_158 sim_159 sim_160 sim_161 sim_162 sim_163
#> runif_1.1       0       0       0       0       0       0       0       0
#> runif_1.2       0       0       0       0       0       0       0       0
#>           sim_164 sim_165 sim_166 sim_167 sim_168 sim_169 sim_170 sim_171
#> runif_1.1       0       0       0       0       0       0       0       0
#> runif_1.2       0       0       0       0       0       0       0       0
#>           sim_172 sim_173 sim_174 sim_175 sim_176 sim_177 sim_178 sim_179
#> runif_1.1       0       0       0       0       0       0       0       0
#> runif_1.2       0       0       0       0       0       0       0       0
#>           sim_180
#> runif_1.1       0
#> runif_1.2       0
```
