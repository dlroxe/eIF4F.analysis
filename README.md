
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eIF4F.analysis

<!-- badges: start -->
<!-- badges: end -->

The goal of eIF4F.analysis is to …

## Installation

You can install the development version of eIF4F.analysis from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("a3609640/eIF4F.analysis")
```

## Before you begin

Perform the following steps to ensure the proper operation of this
package.

[System requirements](#system-requirements)

[Install RStudio/R](#install-rstudior)

[Download data-sets](#download-data-sets)

[File directories](#file-directories)

[Install dependent libraries](#install-dependent-libraries)

[Load package](#load-package)

[Session information](#session-information)

## System requirements

This project makes use of various resource-intensive R packages, which
carries relatively high demands for compute, RAM, and disk I/O resources
(but not for graphics resources). Nonetheless, the necessary hardware is
attainable in high-end consumer-grade systems.

### Description of development systems

The following systems have been used to execute the R scripts in this
project:

1.  (verified) System76 “Serval” mobile workstation

    -   Intel i7-8700k CPU
    -   64GB RAM (DDR4-3000, non-ECC)
    -   Samsung NVMe Pro SSD
    -   Pop!\_OS 20.04 LTS
    -   RStudio
    -   R 4.1.1

2.  (verified) PowerSpec G460 desktop computer

    -   Intel i7-8700k CPU
    -   64GB RAM (DDR4-3200, non-ECC)
    -   Intel M.2 SATA SSD
    -   Samsung NVMe Evo+ SSD
    -   Windows 10 Pro
    -   RStudio for Windows
    -   R 4.0.3

Additional details of these environments are provided in the “Session
Information” section below.

## Install RStudio/R

1.  Download & install R 4.1, if not already installed.
2.  Download & install RStudio, if not already installed.
    <https://www.rstudio.com/products/rstudio/download/>

## Download data-sets

Run `Download.sh` from `Script` folder of our GitHub repository.

Download.sh

``` bash
bash ~/github/eIF4F.analysis/Script/Download.sh
```

`Download.sh` is a bash script to download all needed datasets (TCGA,
GTEx, CPTAC, CCLE and etc) from URLs and unzip them. `Download.sh` will
create the `~/Downloads/EIF_data` directory to store all downloaded
datasets.

``` bash
#!/bin/sh


## download all datasets from the following weblinks

### create the directory to store all downloaded datasets
readonly DATA_FILE_DIRECTORY="${HOME}/Downloads/EIF_data"
mkdir -p "${DATA_FILE_DIRECTORY}"


### TCGA and GTEX DATA
#### TCGA CNV dataset (thresholded)
wget https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz -P "${DATA_FILE_DIRECTORY}"

#### TCGA CNV dataset
wget https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz -P "${DATA_FILE_DIRECTORY}"

#### TCGA CNV ratio dataset
wget https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz -P "${DATA_FILE_DIRECTORY}"

#### TCGA RNA-Seq dataset
wget https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz -P "${DATA_FILE_DIRECTORY}"

#### TCGA sample type annotation
wget https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz -P "${DATA_FILE_DIRECTORY}"

#### TCGA OS data ##
wget https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp -P "${DATA_FILE_DIRECTORY}"

#### TCGA and GTEX RNA-Seq dataset
wget https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz -P "${DATA_FILE_DIRECTORY}"

#### TCGA and GTEX sample type annotation
wget https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz -P "${DATA_FILE_DIRECTORY}"

### CPTAC DATA

#### CPTAC LUAD RNA-Seq data (Gillette et al., 2020)
wget https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/RNA.xlsx -P "${DATA_FILE_DIRECTORY}"

#### CPTAC LUAD Proteomics (Gillette et al., 2020)
wget https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/Protein.xlsx -P "${DATA_FILE_DIRECTORY}"

#### CPTAC LUAD Proteomics
wget https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_Proteome_CDAP_Protein_Report.r1/CPTAC3_Lung_Adeno_Carcinoma_Proteome.tmt10.tsv -P "${DATA_FILE_DIRECTORY}"

#### CPTAC LUAD Phosproteomics (Gillette et al., 2020)
wget https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/Phos.xlsx -P "${DATA_FILE_DIRECTORY}"

#### CPTAC LUAD Sample Annotation
wget https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx -P "${DATA_FILE_DIRECTORY}"

#### CPTAC Clinical Data
wget https://cptc-xfer.uis.georgetown.edu/publicData/Phase_III_Data/CPTAC_LUAD_S046/CPTAC_LUAD_metadata/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx -P "${DATA_FILE_DIRECTORY}"

### CCLE DATA

#### CCLE RNA-Seq data from DepMap Public 20Q4 20Q3
#wget -O CCLE_expression_full.csv https://ndownloader.figshare.com/files/#27902097 -P "${DATA_FILE_DIRECTORY}" #DepMap Public 21Q2

wget https://ndownloader.figshare.com/files/24613349 -O "${DATA_FILE_DIRECTORY}/CCLE_expression_full.csv" #DepMap Public 20Q3

#### CCLE annotation data
#wget -O sample_info.csv https://ndownloader.figshare.com/files/27902376 -P #"${DATA_FILE_DIRECTORY}" #DepMap Public 21Q2

wget https://ndownloader.figshare.com/files/24613394 -O "${DATA_FILE_DIRECTORY}/sample_info.csv" #DepMap Public 20Q3

#### CCLE proteomics data
wget https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz -P "${DATA_FILE_DIRECTORY}"


gunzip ${DATA_FILE_DIRECTORY}/*.gz
```

**CRITICAL**: If the root directory path `~/Downloads/EIF_data` does not
suit, they may be adjusted trivially in these lines near the top of the
`Download.sh` script.

## File directories

Confirm the directories for input and output files. The directories for
input and output files are defined in the `Load.R` as following.

    data_file_directory <- "~/Downloads/EIF_data" 
    output_directory <- "~/Documents/EIF_output"

**CRITICAL**: If the root directory paths `~/Download/EIF_data` and
`~/Documents/EIF_output` do not suit, they may be adjusted trivially in
these lines near the top of the `Download.sh` and `Load.R` scripts.

## Install dependent libraries

The work here depends upon many R libraries. The following command may
be a useful way to install them all.

``` r
# use Bioconductor version 3.14 for package installation
BiocManager::install(version = "3.14")

bio_pkgs <- c(
  "AnnotationDbi", "circlize",  "clusterProfiler", "ComplexHeatmap", "corrplot", 
  "data.table", "dplyr", "EnvStats", "eulerr", "factoextra", "FactoMineR", 
  "forcats", "forestplot", "ggfortify", "ggplot2", "ggpubr",  "graphics", 
  "grDevices", "grid", "limma", "missMDA", "org.Hs.eg.db", "purrr", 
  "RColorBrewer", "ReactomePA", "readr",  "readxl", "reshape2", "scales", 
  "stats", "stringr", "survival", "survivalAnalysis", "tibble", "tidyr", 
  "tidyselect")

# install required packages
BiocManager::install(bio_pkgs)

# load required packages
lapply(bio_pkgs, require, character.only = TRUE)
```

## Load package

Load `eIF4F.analysis` in the R console.

``` r
library(eIF4F.analysis)
```

## Session information

The version information of R, Linux and attached or loaded packages for
developing this package is the following.

``` r
> sessionInfo(package = "eIF4F.analysis")
R version 4.1.1 (2021-08-10)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Pop!_OS 20.04 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] AnnotationDbi_1.56.2 IRanges_2.28.0       S4Vectors_0.32.3     Biobase_2.54.0       BiocGenerics_0.40.0  eIF4F.analysis_0.1.0

loaded via a namespace (and not attached):
  [1] utf8_1.2.2             tidyselect_1.1.1       htmlwidgets_1.5.4      RSQLite_2.2.9          FactoMineR_2.4        
  [6] grid_4.1.1             BiocParallel_1.28.3    devtools_2.4.3         scatterpie_0.1.7       munsell_0.5.0         
 [11] codetools_0.2-18       DT_0.20                withr_2.4.3            colorspace_2.0-2       GOSemSim_2.20.0       
 [16] ggfortify_0.4.13       highr_0.9              knitr_1.37             leaps_3.1              rstudioapi_0.13       
 [21] ggsignif_0.6.3         DOSE_3.20.1            labeling_0.4.2         GenomeInfoDbData_1.2.7 KMsurv_0.1-5          
 [26] polyclip_1.10-0        bit64_4.0.5            farver_2.1.0           rprojroot_2.0.2        downloader_0.4        
 [31] survivalAnalysis_0.2.0 vctrs_0.3.8            treeio_1.18.1          generics_0.1.1         xfun_0.29             
 [36] R6_2.5.1               doParallel_1.0.16      GenomeInfoDb_1.30.0    clue_0.3-60            graphlayouts_0.7.2    
 [41] bitops_1.0-7           cachem_1.0.6           fgsea_1.20.0           gridGraphics_0.5-1     assertthat_0.2.1      
 [46] vroom_1.5.7            scales_1.1.1           ggraph_2.0.5           enrichplot_1.14.1      gtable_0.3.0          
 [51] processx_3.5.2         tidygraph_1.2.0        rlang_0.4.12           scatterplot3d_0.3-41   eulerr_6.1.1          
 [56] GlobalOptions_0.1.2    splines_4.1.1          rstatix_0.7.0          lazyeval_0.2.2         broom_0.7.10          
 [61] checkmate_2.0.0        reshape2_1.4.4         abind_1.4-5            backports_1.4.1        qvalue_2.26.0         
 [66] clusterProfiler_4.2.1  tools_4.1.1            usethis_2.1.5          tidytidbits_0.2.3      ggplotify_0.1.0       
 [71] ggplot2_3.3.5          ellipsis_0.3.2         RColorBrewer_1.1-2     sessioninfo_1.2.2      Rcpp_1.0.7            
 [76] plyr_1.8.6             zlibbioc_1.40.0        purrr_0.3.4            RCurl_1.98-1.5         ps_1.6.0              
 [81] prettyunits_1.1.1      ggpubr_0.4.0           GetoptLong_1.0.5       viridis_0.6.2          cowplot_1.1.1         
 [86] zoo_1.8-9              ggrepel_0.9.1          cluster_2.1.2          factoextra_1.0.7       fs_1.5.2              
 [91] magrittr_2.0.1         data.table_1.14.2      DO.db_2.9              forestplot_2.0.1       circlize_0.4.13       
 [96] reactome.db_1.77.0     survminer_0.4.9        mvtnorm_1.1-3          matrixStats_0.61.0     pkgload_1.2.4         
[101] evaluate_0.14          hms_1.1.1              patchwork_1.1.1        xtable_1.8-4           readxl_1.3.1          
[106] missMDA_1.18           gridExtra_2.3          shape_1.4.6            testthat_3.1.1         compiler_4.1.1        
[111] mice_3.14.0            tibble_3.1.6           crayon_1.4.2           shadowtext_0.0.9       htmltools_0.5.2       
[116] ggfun_0.0.4            tzdb_0.2.0             tidyr_1.1.4            aplot_0.1.1            ReactomePA_1.38.0     
[121] DBI_1.1.1              tweenr_1.0.2           corrplot_0.92          EnvStats_2.4.0         ComplexHeatmap_2.10.0 
[126] rappdirs_0.3.3         MASS_7.3-54            Matrix_1.4-0           car_3.0-12             readr_2.1.1           
[131] cli_3.1.0              parallel_4.1.1         igraph_1.2.10          forcats_0.5.1          pkgconfig_2.0.3       
[136] km.ci_0.5-2            flashClust_1.01-2      foreach_1.5.1          ggtree_3.2.1           XVector_0.34.0        
[141] yulab.utils_0.0.4      stringr_1.4.0          callr_3.7.0            digest_0.6.29          graph_1.72.0          
[146] Biostrings_2.62.0      rmarkdown_2.11         polylabelr_0.2.0       cellranger_1.1.0       fastmatch_1.1-3       
[151] survMisc_0.5.5         tidytree_0.3.6         graphite_1.40.0        rjson_0.2.20           lifecycle_1.0.1       
[156] nlme_3.1-153           jsonlite_1.7.2         carData_3.0-4          desc_1.4.0             viridisLite_0.4.0     
[161] limma_3.50.0           fansi_0.5.0            pillar_1.6.4           lattice_0.20-45        KEGGREST_1.34.0       
[166] fastmap_1.1.0          httr_1.4.2             pkgbuild_1.3.0         survival_3.2-13        GO.db_3.14.0          
[171] glue_1.6.0             remotes_2.4.2          png_0.1-7              iterators_1.0.13       bit_4.0.4             
[176] ggforce_0.3.3          stringi_1.7.6          blob_1.2.2             org.Hs.eg.db_3.14.0    memoise_2.0.1         
[181] dplyr_1.0.7            ape_5.5
```

## **Tutorials**

Open the `Analysis.R` from `Script` folder from our GitHub repository.
This script contains the command lines to execute all analyses.
