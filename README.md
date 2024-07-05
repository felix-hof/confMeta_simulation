# Combined confidence curves for meta-analysis 

This repository contains:

  - **simulation/**: Simulation study comparing *p*-value combination methods for meta-analysis.
  - **shiny/**: A shiny app that visualizes the simulation results. Currently under construction.
  - **run_pipeline.sh**: Shell script that runs the simulation and pushes the resulting files to Github.

## sessionInfo output of the last simulation run

``` r
#> R version 4.4.0 (2024-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Debian GNU/Linux 12 (bookworm)
#> 
#> Matrix products: default
#> BLAS:   /usr/local/lib/R/lib/libRblas.so 
#> LAPACK: /usr/local/lib/R/lib/libRlapack.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_CH.UTF-8        LC_COLLATE=en_GB.UTF-8    
#>  [5] LC_MONETARY=de_CH.UTF-8    LC_MESSAGES=en_GB.UTF-8   
#>  [7] LC_PAPER=de_CH.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_CH.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Zurich
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] parallel  stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#> [1] RhpcBLASctl_0.23-42 doRNG_1.8.6         rngtools_1.5.2     
#> [4] doParallel_1.0.17   iterators_1.0.14    foreach_1.5.2      
#> [7] confMeta_0.3.1     
#> 
#> loaded via a namespace (and not attached):
#>  [1] vctrs_0.6.5      cli_3.6.3        rlang_1.1.4      generics_0.1.3  
#>  [5] glue_1.7.0       colorspace_2.1-0 scales_1.3.0     fansi_1.0.6     
#>  [9] grid_4.4.0       munsell_0.5.1    tibble_3.2.1     lifecycle_1.0.4 
#> [13] compiler_4.4.0   dplyr_1.1.4      codetools_0.2-20 pkgconfig_2.0.3 
#> [17] digest_0.6.35    R6_2.5.1         tidyselect_1.2.1 utf8_1.2.4      
#> [21] pillar_1.9.0     magrittr_2.0.3   tools_4.4.0      gtable_0.3.5    
#> [25] ggplot2_3.5.1    remotes_2.5.0 
```
