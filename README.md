# A comparison of combined p-value functions for meta-analysis

This repository contains:

  - **simulation/**: Simulation study comparing *p*-value combination methods for meta-analysis.
  - **run_pipeline.sh**: Shell script that runs the simulation and pushes the resulting files to Github.

## sessionInfo output of the last simulation run

``` r
#> R version 4.4.1 (2024-06-14)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Debian GNU/Linux trixie/sid
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
#> [7] confMeta_0.4.2     
#> 
#> loaded via a namespace (and not attached):
#>  [1] generics_0.1.3           xml2_1.3.6               stringi_1.8.4           
#>  [4] lattice_0.22-6           digest_0.6.35            lme4_1.1-36             
#>  [7] hms_1.1.3                magrittr_2.0.3           grid_4.4.1              
#> [10] meta_8.0-2               CompQuadForm_1.4.3       Matrix_1.7-0            
#> [13] purrr_1.0.2              scales_1.3.0             codetools_0.2-20        
#> [16] numDeriv_2016.8-1.1      reformulas_0.4.0         Rdpack_2.6.2            
#> [19] cli_3.6.3                rlang_1.1.5              rbibutils_2.3           
#> [22] ReplicationSuccess_1.3.3 munsell_0.5.1            splines_4.4.1           
#> [25] tools_4.4.1              tzdb_0.4.0               nloptr_2.1.1            
#> [28] minqa_1.2.8              metafor_4.6-0            dplyr_1.1.4             
#> [31] colorspace_2.1-1         ggplot2_3.5.1            mathjaxr_1.6-0          
#> [34] boot_1.3-30              vctrs_0.6.5              R6_2.5.1                
#> [37] lifecycle_1.0.4          stringr_1.5.1            MASS_7.3-60.2           
#> [40] pkgconfig_2.0.3          pillar_1.10.1            gtable_0.3.6            
#> [43] glue_1.8.0               Rcpp_1.0.14              tibble_3.2.1            
#> [46] tidyselect_1.2.1         farver_2.1.2             nlme_3.1-164            
#> [49] patchwork_1.3.0          readr_2.1.5              compiler_4.4.1          
#> [52] metadat_1.2-0
```
