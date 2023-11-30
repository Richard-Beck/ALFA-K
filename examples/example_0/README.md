Visualizing karyotype frequencies in the input data and overlay ALFA-K
fit results.

``` r
library(ggplot2)
source("utils/visualisation_functions.R")
```

melt_for_plotting() plots karyotype frequencies over time, with option
to overlay alfa-k fitted frequency estimates. See example_1 for more
info on generating suitable input data for this function.

**Arguments**

**data_obj:** karyotype frequency data formatted as output from
proc_sim() function

**nclones:** number of karyotypes to plot (sorted by descending
frequency in input data)

**fit_obj:** (optional) fit resulting from alfak() function

**Value**

melt_for_plotting returns a list that is convenient for input into
ggplot

**$data:** data frame containing karyotype frequencies from input data

**$fit:** either NULL or a dataframe containing karyotype frequency
estimates from ALFA-K

example:

``` r
data_obj <- readRDS("example_data/SA535.Rds")
fit_obj <- readRDS("example_data/SA535_landscape.Rds")

x1 <- melt_for_plotting(data_obj,nclones=5)
x2 <- melt_for_plotting(data_obj,nclones=5,fit_obj)


p <- ggplot(x1$data,aes(x=time,y=frequency,color=karyotype))+
  geom_point(size=2)+
  geom_line()+
  theme_classic(base_size=12)
p
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
p <- ggplot(x2$data,aes(x=time,y=frequency,color=karyotype))+
  geom_point(size=2)+
  geom_line(data=x2$fit)+
  theme_classic(base_size=12)
p
```

![](README_files/figure-markdown_github/unnamed-chunk-2-2.png)
