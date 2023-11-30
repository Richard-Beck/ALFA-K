Data in this folder was derived from longitudinal sequencing data
obtained from *Salehi, Sohrab, et al. “Clonal fitness inferred from
time-series modelling of single-cell cancer genomes.” Nature 595.7868
(2021): 585-590.* In particular, Extended Data Fig. 1 in their paper is
a useful reference for relating the data in this folder to the
underlying experimental procedures.

The data objects provided here contain karyotype frequency data
formatted as follows:

``` r
x <- readRDS("example_data/hTERTb.Rds")
head(x$x)
```

    ##                                              20  30  35  40  45  50  55
    ## 2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2 748 308 322 184 153 120 134
    ## 2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.3.2.2.3.2.2   6 235 387 194 127  82  15
    ## 2.2.2.1.2.2.2.2.2.2.2.2.2.2.2.2.2.2.1.2.2.2   0   0   4   7 203 297 238
    ## 2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.1  27  34  38  17  35  42 136
    ## 2.1.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2  74  63  76  14   3   3   2
    ## 2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.2.3.2.2  19  25  33  16  22  21  12

Note that the column names refer to the passage number as indicated in
Salehi et. al. The data we received contained copy number calls for a
multitude of genomic segments per chromosome. We reduced these data to
whole-chromosome resolution by assigning each chromosome the modal copy
number of all of it’s segments.

For more information about the formatting of these data, see example_1
