[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/release/sung/POPS-Placenta-Transcriptome-2020.svg)](https://GitHub.com/sung/POPS-Placenta-Transcriptome-2020/releases/)
![GitHub all releases](https://img.shields.io/github/downloads/sung/POPS-Placenta-Transcriptome-2020/total)
![GitHub last commit ](https://img.shields.io/github/last-commit/sung/POPS-Placenta-Transcriptome-2020)

## Introduction
This is a collection of custom scripts used for the following paper: "The RNA landscape of the human placenta in health and disease, Gong *et al.* 2021. *Nat Comm*, 2021". A [Shiny app](https://www.obgyn.cam.ac.uk/placentome) to browse the placenta transcriptome is also available and the source code is [here](https://github.com/sung/ShinyPlacentome).

## Initial mapping of RNA-Seq data
The RNA-Seq reads were trimmed (using [cutadapt](https://github.com/marcelm/cutadapt) and [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)) and mapped to the primary chromosomal assemblies of the [GRCh38.p3 version](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.29/) of the human reference genome using [TopHat2](https://github.com/infphilo/tophat), a splice-aware mapper built on top of [Bowtie2](https://github.com/BenLangmead/bowtie2) short-read aligner. For more details, read [Gong et al., Epigenetics, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5989156/) and [Gong et al., JCI Insight, 2018](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6124516/).

## `SessionInfo()` shown below (or [here](sessionInfo.txt))
```r
R version 3.6.1 (2019-07-05)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Scientific Linux 7.8 (Nitrogen)

Matrix products: default
BLAS:   /usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/r-3.6.1-zrytncqvsnw5h4dl6t6njefj7otl4bg4/rlib/R/lib/libRblas.so
LAPACK: /usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/r-3.6.1-zrytncqvsnw5h4dl6t6njefj7otl4bg4/rlib/R/lib/libRlapack.so

locale:
 [1] LC_CTYPE=en_GB.UTF-8       LC_NUMERIC=C               LC_TIME=en_GB.UTF-8        LC_COLLATE=en_GB.UTF-8    
 [5] LC_MONETARY=en_GB.UTF-8    LC_MESSAGES=en_GB.UTF-8    LC_PAPER=en_GB.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_GB.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scales_1.1.0         ggthemes_4.2.0       GenomicRanges_1.38.0 GenomeInfoDb_1.22.0  IRanges_2.20.2      
 [6] S4Vectors_0.24.3     BiocGenerics_0.32.0  ggplot2_3.2.1        RColorBrewer_1.1-2   nvimcom_0.9-83      
[11] data.table_1.12.8    colorout_1.2-2      

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.3             pillar_1.4.3           compiler_3.6.1         XVector_0.26.0         bitops_1.0-6          
 [6] tools_3.6.1            zlibbioc_1.32.0        lifecycle_0.1.0        tibble_2.1.3           gtable_0.3.0          
[11] pkgconfig_2.0.3        rlang_0.4.4            GenomeInfoDbData_1.2.2 withr_2.1.2            dplyr_0.8.4           
[16] stringr_1.4.0          grid_3.6.1             tidyselect_1.0.0       glue_1.3.1             R6_2.4.1              
[21] purrr_0.3.3            magrittr_1.5           assertthat_0.2.1       colorspace_1.4-1       stringi_1.4.5         
[26] RCurl_1.98-1.1         lazyeval_0.2.2         munsell_0.5.0          crayon_1.3.4          
```

## Contacts
+ [Sung Gong](https://www.obgyn.cam.ac.uk/staff/research-staff/sung-gong/)
