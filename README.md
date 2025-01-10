# mi-using-mars

This repository contains the R code used to carry out the simulation study reported in the following paper:

**Sepin, Jerome. "Multiple Imputation Using Multivariate Adaptive Regression Splines."**

## Files and folder structure:

-   **simulation**

    -   `mice.impute.mars.R`: Contains the MARS function to be used within the `mice` package.

    -   `functions_for_sim.R`: Contains functions for data generation and missingness generation.

    -   `simulation.R`: Script to execute the simulation study.

    -   `visualization.R`: Script to generate figures for the manuscript.

    -   **results**: Folder containing simulation results and parameters for the study.

    -   **figures**: Folder containing the generated figures.

<!-- -->

-   **manuscript**: Contains a Quarto project for the manuscript.

## Reproducibility:

To reproduce the figures in the manuscript, run the following scripts (in this order):

-   `simulation.R`

-   `visualization.R`

## Software Requirements:

> sessionInfo() R version 4.4.2 (2024-10-31) Platform: x86_64-pc-linux-gnu Running under: Ubuntu 24.04.1 LTS

```         
Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0

locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8       
 [4] LC_COLLATE=de_DE.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Paris
tzcode source: system (glibc)

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.3.0    ggpubr_0.6.0    mixgb_1.0.2     zeallot_0.1.0   mvtnorm_1.3-2  
 [6] mice_3.17.0     lubridate_1.9.3 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
[11] purrr_1.0.2     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.1  
[16] tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] gtable_0.3.6       shape_1.4.6.1      rstatix_0.7.2      lattice_0.22-5    
 [5] tzdb_0.4.0         vctrs_0.6.5        tools_4.4.2        generics_0.1.3    
 [9] fansi_1.0.6        pan_1.9            pkgconfig_2.0.3    jomo_2.7-6        
[13] Matrix_1.7-1       data.table_1.16.2  RcppParallel_5.1.7 lifecycle_1.0.4   
[17] compiler_4.4.2     munsell_0.5.1      codetools_0.2-20   carData_3.0-5     
[21] glmnet_4.1-8       car_3.1-2          pillar_1.9.0       nloptr_2.1.1      
[25] MASS_7.3-61        iterators_1.0.14   abind_1.4-5        rpart_4.1.23      
[29] boot_1.3-30        foreach_1.5.2      mitml_0.4-5        nlme_3.1-165      
[33] tidyselect_1.2.1   stringi_1.8.4      splines_4.4.2      grid_4.4.2        
[37] colorspace_2.1-1   cli_3.6.3          magrittr_2.0.3     Rfast_2.1.0       
[41] survival_3.7-0     utf8_1.2.4         broom_1.0.7        withr_3.0.2       
[45] backports_1.5.0    RcppZiggurat_0.1.6 xgboost_1.7.8.1    timechange_0.3.0  
[49] nnet_7.3-19        lme4_1.1-35.1      ggsignif_0.6.4     hms_1.1.3         
[53] rlang_1.1.4        Rcpp_1.0.13-1      glue_1.8.0         rstudioapi_0.17.1 
[57] minqa_1.2.8        jsonlite_1.8.9     R6_2.5.1     
```
