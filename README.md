
<!-- README.md is generated from README.Rmd. Please edit that file -->

# eIF4F.analysis

<!-- badges: start -->
<!-- badges: end -->

The goal of eIF4F.analysis is to understand function and regulation of
interactions among translation initiation complex proteins across tumor
types

## Before you begin

Perform the following steps to ensure the proper operation of this
package.

[System requirements](#system-requirements)

[Install RStudio/R](#install-rstudior)

[Install dependent libraries](#install-dependent-libraries)

[Install eIF4F analysis package](#install-eif4f-analysis-package)

[Download datasets](#download-datasets)

[File directories](#file-directories)

[Session information](#session-information)

[Tutorials](#tutorials)

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

## Install dependent libraries

The work here depends upon many R libraries. The following command may
be a useful way to install them all.

``` r
# use Bioconductor version 3.14 for package installation
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")

# install required packages
bio_pkgs <- c(
  "AnnotationDbi", "circlize",  "clusterProfiler", "ComplexHeatmap", "corrplot", 
  "data.table", "devtools", "dplyr", "EnvStats", "eulerr", "factoextra", 
  "FactoMineR", "forcats", "forestplot", "ggfortify", "ggplot2", "ggpubr", 
  "graphics", "grDevices", "grid", "limma", "missMDA", "org.Hs.eg.db", "purrr", 
  "RColorBrewer", "ReactomePA", "readr",  "readxl", "reshape2", "scales", 
  "stats", "stringr", "survival", "survivalAnalysis", "tibble", "tidyr", 
  "tidyselect")
BiocManager::install(bio_pkgs)

# load required packages
lapply(bio_pkgs, require, character.only = TRUE)
```

## Install eIF4F analysis package

You can install the development version of `eIF4F.analysis` from
[GitHub](https://github.com/) and load it in the R console.

``` r
# Install eIF4F.analysis package 
devtools::install_github("a3609640/eIF4F.analysis")

# Load eIF4F.analysis package 
library(eIF4F.analysis)
```

## Download datasets

Run `Download.sh` from `Script` folder of our GitHub repository.

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

#### TCGA OS data
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
wget https://ndownloader.figshare.com/files/24613349 -O "${DATA_FILE_DIRECTORY}/CCLE_expression_full.csv" #DepMap Public 20Q3

#### CCLE annotation data
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

``` r
data_file_directory <- "~/Downloads/EIF_data" 
output_directory <- "~/Documents/EIF_output"
```

**CRITICAL**: If the root directory paths `~/Download/EIF_data` and
`~/Documents/EIF_output` do not suit, they may be adjusted trivially in
these lines near the top of the `Download.sh` and `Load.R` scripts.

## Session information

The version information of R, Linux and attached or loaded packages for
developing this package is the following.

``` r
─ Session info ─────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.1.1 (2021-08-10)
 os       Pop!_OS 20.04 LTS
 system   x86_64, linux-gnu
 ui       RStudio
 language en_US:en
 collate  en_US.UTF-8
 ctype    en_US.UTF-8
 tz       America/New_York
 date     2021-12-29
 rstudio  2021.09.0-preview+341 Ghost Orchid (desktop)
 pandoc   2.14.0.3 @ /usr/lib/rstudio/bin/pandoc/ (via rmarkdown)

─ Packages ─────────────────────────────────────────────────────────────────────
 ! package          * version  date (UTC) lib source
   abind              1.4-5    2016-07-21 [1] CRAN (R 4.1.1)
   AnnotationDbi    * 1.56.2   2021-11-09 [1] Bioconductor
   ape                5.6      2021-12-21 [1] CRAN (R 4.1.1)
   aplot              0.1.1    2021-09-22 [1] CRAN (R 4.1.1)
   assertthat         0.2.1    2019-03-21 [1] CRAN (R 4.1.1)
   backports          1.4.1    2021-12-13 [1] CRAN (R 4.1.1)
   Biobase          * 2.54.0   2021-10-26 [1] Bioconductor
   BiocGenerics     * 0.40.0   2021-10-26 [1] Bioconductor
   BiocParallel       1.28.3   2021-12-09 [1] Bioconductor
   Biostrings         2.62.0   2021-10-26 [1] Bioconductor
   bit                4.0.4    2020-08-04 [1] CRAN (R 4.1.1)
   bit64              4.0.5    2020-08-30 [1] CRAN (R 4.1.1)
   bitops             1.0-7    2021-04-24 [1] CRAN (R 4.1.1)
   blob               1.2.2    2021-07-23 [1] CRAN (R 4.1.1)
   broom              0.7.10   2021-10-31 [1] CRAN (R 4.1.1)
   cachem             1.0.6    2021-08-19 [1] CRAN (R 4.1.1)
   callr              3.7.0    2021-04-20 [1] CRAN (R 4.1.1)
   car                3.0-12   2021-11-06 [1] CRAN (R 4.1.1)
   carData            3.0-4    2020-05-22 [1] CRAN (R 4.1.1)
   cellranger         1.1.0    2016-07-27 [1] CRAN (R 4.1.1)
   checkmate          2.0.0    2020-02-06 [1] CRAN (R 4.1.1)
   circlize           0.4.13   2021-06-09 [1] CRAN (R 4.1.1)
   cli                3.1.0    2021-10-27 [1] CRAN (R 4.1.1)
   clipr              0.7.1    2020-10-08 [1] CRAN (R 4.1.1)
   clue               0.3-60   2021-10-11 [1] CRAN (R 4.1.1)
   cluster            2.1.2    2021-04-17 [4] CRAN (R 4.0.5)
   clusterProfiler    4.2.1    2021-12-14 [1] Bioconductor
   codetools          0.2-18   2020-11-04 [4] CRAN (R 4.0.3)
   colorspace         2.0-2    2021-06-24 [1] CRAN (R 4.1.1)
   ComplexHeatmap     2.10.0   2021-10-26 [1] Bioconductor
   corrplot           0.92     2021-11-18 [1] CRAN (R 4.1.1)
   cowplot            1.1.1    2020-12-30 [1] CRAN (R 4.1.1)
   crayon             1.4.2    2021-10-29 [1] CRAN (R 4.1.1)
   data.table         1.14.2   2021-09-27 [1] CRAN (R 4.1.1)
   DBI                1.1.2    2021-12-20 [1] CRAN (R 4.1.1)
   desc               1.4.0    2021-09-28 [1] CRAN (R 4.1.1)
   details          * 0.2.1    2020-01-12 [1] CRAN (R 4.1.1)
   devtools           2.4.3    2021-11-30 [1] CRAN (R 4.1.1)
   digest             0.6.29   2021-12-01 [1] CRAN (R 4.1.1)
   DO.db              2.9      2021-09-14 [1] Bioconductor
   doParallel         1.0.16   2020-10-16 [1] CRAN (R 4.1.1)
   DOSE               3.20.1   2021-11-18 [1] Bioconductor
   downloader         0.4      2015-07-09 [1] CRAN (R 4.1.1)
   dplyr              1.0.7    2021-06-18 [1] CRAN (R 4.1.1)
   DT                 0.20     2021-11-15 [1] CRAN (R 4.1.1)
   ellipsis           0.3.2    2021-04-29 [1] CRAN (R 4.1.1)
   enrichplot         1.14.1   2021-10-31 [1] Bioconductor
   EnvStats           2.4.0    2020-10-21 [1] CRAN (R 4.1.1)
   eulerr             6.1.1    2021-09-06 [1] CRAN (R 4.1.1)
   evaluate           0.14     2019-05-28 [1] CRAN (R 4.1.1)
   factoextra         1.0.7    2020-04-01 [1] CRAN (R 4.1.1)
   FactoMineR         2.4      2020-12-11 [1] CRAN (R 4.1.1)
   fansi              0.5.0    2021-05-25 [1] CRAN (R 4.1.1)
   farver             2.1.0    2021-02-28 [1] CRAN (R 4.1.1)
   fastmap            1.1.0    2021-01-25 [1] CRAN (R 4.1.1)
   fastmatch          1.1-3    2021-07-23 [1] CRAN (R 4.1.1)
   fgsea              1.20.0   2021-10-26 [1] Bioconductor
   flashClust         1.01-2   2012-08-21 [1] CRAN (R 4.1.1)
   forcats            0.5.1    2021-01-27 [1] CRAN (R 4.1.1)
   foreach            1.5.1    2020-10-15 [1] CRAN (R 4.1.1)
   forestplot         2.0.1    2021-09-03 [1] CRAN (R 4.1.1)
   fs                 1.5.2    2021-12-08 [1] CRAN (R 4.1.1)
   generics           0.1.1    2021-10-25 [1] CRAN (R 4.1.1)
   GenomeInfoDb       1.30.0   2021-10-26 [1] Bioconductor
   GenomeInfoDbData   1.2.7    2021-12-17 [1] Bioconductor
   GetoptLong         1.0.5    2020-12-15 [1] CRAN (R 4.1.1)
   ggforce            0.3.3    2021-03-05 [1] CRAN (R 4.1.1)
   ggfortify          0.4.13   2021-10-25 [1] CRAN (R 4.1.1)
   ggfun              0.0.4    2021-09-17 [1] CRAN (R 4.1.1)
   ggplot2            3.3.5    2021-06-25 [1] CRAN (R 4.1.1)
   ggplotify          0.1.0    2021-09-02 [1] CRAN (R 4.1.1)
   ggpubr             0.4.0    2020-06-27 [1] CRAN (R 4.1.1)
   ggraph             2.0.5    2021-02-23 [1] CRAN (R 4.1.1)
   ggrepel            0.9.1    2021-01-15 [1] CRAN (R 4.1.1)
   ggsignif           0.6.3    2021-09-09 [1] CRAN (R 4.1.1)
   ggtree             3.2.1    2021-11-16 [1] Bioconductor
   GlobalOptions      0.1.2    2020-06-10 [1] CRAN (R 4.1.1)
   glue               1.6.0    2021-12-17 [1] CRAN (R 4.1.1)
   GO.db              3.14.0   2021-12-17 [1] Bioconductor
   GOSemSim           2.20.0   2021-10-26 [1] Bioconductor
   graph              1.72.0   2021-10-26 [1] Bioconductor
   graphite           1.40.0   2021-10-26 [1] Bioconductor
   graphlayouts       0.7.2    2021-11-21 [1] CRAN (R 4.1.1)
   gridExtra          2.3      2017-09-09 [1] CRAN (R 4.1.1)
   gridGraphics       0.5-1    2020-12-13 [1] CRAN (R 4.1.1)
   gtable             0.3.0    2019-03-25 [1] CRAN (R 4.1.1)
   hms                1.1.1    2021-09-26 [1] CRAN (R 4.1.1)
   htmltools          0.5.2    2021-08-25 [1] CRAN (R 4.1.1)
   htmlwidgets        1.5.4    2021-09-08 [1] CRAN (R 4.1.1)
   httr               1.4.2    2020-07-20 [1] CRAN (R 4.1.1)
   igraph             1.2.10   2021-12-15 [1] CRAN (R 4.1.1)
   IRanges          * 2.28.0   2021-10-26 [1] Bioconductor
   iterators          1.0.13   2020-10-15 [1] CRAN (R 4.1.1)
   jsonlite           1.7.2    2020-12-09 [1] CRAN (R 4.1.1)
   KEGGREST           1.34.0   2021-10-26 [1] Bioconductor
   km.ci              0.5-2    2009-08-30 [1] CRAN (R 4.1.1)
   KMsurv             0.1-5    2012-12-03 [1] CRAN (R 4.1.1)
   knitr              1.37     2021-12-16 [1] CRAN (R 4.1.1)
   labeling           0.4.2    2020-10-20 [1] CRAN (R 4.1.1)
   lattice            0.20-45  2021-09-22 [1] CRAN (R 4.1.1)
   lazyeval           0.2.2    2019-03-15 [1] CRAN (R 4.1.1)
   leaps              3.1      2020-01-16 [1] CRAN (R 4.1.1)
   lifecycle          1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
   limma              3.50.0   2021-10-26 [1] Bioconductor
   magrittr           2.0.1    2020-11-17 [1] CRAN (R 4.1.1)
   MASS               7.3-54   2021-05-03 [4] CRAN (R 4.0.5)
   Matrix             1.4-0    2021-12-08 [1] R-Forge (R 4.1.1)
   matrixStats        0.61.0   2021-09-17 [1] CRAN (R 4.1.1)
   memoise            2.0.1    2021-11-26 [1] CRAN (R 4.1.1)
   mice               3.14.0   2021-11-24 [1] CRAN (R 4.1.1)
   missMDA            1.18     2020-12-11 [1] CRAN (R 4.1.1)
   munsell            0.5.0    2018-06-12 [1] CRAN (R 4.1.1)
   mvtnorm            1.1-3    2021-10-08 [1] CRAN (R 4.1.1)
   nlme               3.1-153  2021-09-07 [1] CRAN (R 4.1.1)
   org.Hs.eg.db       3.14.0   2021-12-17 [1] Bioconductor
   pander           * 0.6.4    2021-06-13 [1] CRAN (R 4.1.1)
   patchwork          1.1.1    2020-12-17 [1] CRAN (R 4.1.1)
   pillar             1.6.4    2021-10-18 [1] CRAN (R 4.1.1)
   pkgbuild           1.3.1    2021-12-20 [1] CRAN (R 4.1.1)
   pkgconfig          2.0.3    2019-09-22 [1] CRAN (R 4.1.1)
   pkgload            1.2.4    2021-11-30 [1] CRAN (R 4.1.1)
   plyr               1.8.6    2020-03-03 [1] CRAN (R 4.1.1)
   png                0.1-7    2013-12-03 [1] CRAN (R 4.1.1)
   polyclip           1.10-0   2019-03-14 [1] CRAN (R 4.1.1)
   polylabelr         0.2.0    2020-04-19 [1] CRAN (R 4.1.1)
   prettyunits        1.1.1    2020-01-24 [1] CRAN (R 4.1.1)
   processx           3.5.2    2021-04-30 [1] CRAN (R 4.1.1)
   ps                 1.6.0    2021-02-28 [1] CRAN (R 4.1.1)
   purrr              0.3.4    2020-04-17 [1] CRAN (R 4.1.1)
   qvalue             2.26.0   2021-10-26 [1] Bioconductor
   R6                 2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
   ragg               1.2.1    2021-12-06 [1] CRAN (R 4.1.1)
   rappdirs           0.3.3    2021-01-31 [1] CRAN (R 4.1.1)
   RColorBrewer       1.1-2    2014-12-07 [1] CRAN (R 4.1.1)
   Rcpp               1.0.7    2021-07-07 [1] CRAN (R 4.1.1)
   RCurl              1.98-1.5 2021-09-17 [1] CRAN (R 4.1.1)
   reactome.db        1.77.0   2021-12-17 [1] Bioconductor
   ReactomePA         1.38.0   2021-10-26 [1] Bioconductor
   readr              2.1.1    2021-11-30 [1] CRAN (R 4.1.1)
   readxl             1.3.1    2019-03-13 [1] CRAN (R 4.1.1)
   remotes            2.4.2    2021-11-30 [1] CRAN (R 4.1.1)
   reshape2           1.4.4    2020-04-09 [1] CRAN (R 4.1.1)
   rjson              0.2.20   2018-06-08 [1] CRAN (R 4.1.1)
   rlang              0.4.12   2021-10-18 [1] CRAN (R 4.1.1)
   rmarkdown          2.11     2021-09-14 [1] CRAN (R 4.1.1)
   rprojroot          2.0.2    2020-11-15 [1] CRAN (R 4.1.1)
   RSQLite            2.2.9    2021-12-06 [1] CRAN (R 4.1.1)
   rstatix            0.7.0    2021-02-13 [1] CRAN (R 4.1.1)
   rstudioapi         0.13     2020-11-12 [1] CRAN (R 4.1.1)
   S4Vectors        * 0.32.3   2021-11-21 [1] Bioconductor
   scales             1.1.1    2020-05-11 [1] CRAN (R 4.1.1)
   scatterpie         0.1.7    2021-08-20 [1] CRAN (R 4.1.1)
   scatterplot3d      0.3-41   2018-03-14 [1] CRAN (R 4.1.1)
   sessioninfo        1.2.2    2021-12-06 [1] CRAN (R 4.1.1)
   shadowtext         0.1.0    2021-12-20 [1] CRAN (R 4.1.1)
   shape              1.4.6    2021-05-19 [1] CRAN (R 4.1.1)
   stringi            1.7.6    2021-11-29 [1] CRAN (R 4.1.1)
   stringr            1.4.0    2019-02-10 [1] CRAN (R 4.1.1)
   survival           3.2-13   2021-08-24 [4] CRAN (R 4.1.1)
   survivalAnalysis   0.2.0    2021-04-24 [1] CRAN (R 4.1.1)
   survminer          0.4.9    2021-03-09 [1] CRAN (R 4.1.1)
   survMisc           0.5.5    2018-07-05 [1] CRAN (R 4.1.1)
   systemfonts        1.0.3    2021-10-13 [1] CRAN (R 4.1.1)
   testthat           3.1.1    2021-12-03 [1] CRAN (R 4.1.1)
   textshaping        0.3.6    2021-10-13 [1] CRAN (R 4.1.1)
   tibble             3.1.6    2021-11-07 [1] CRAN (R 4.1.1)
   tidygraph          1.2.0    2020-05-12 [1] CRAN (R 4.1.1)
   tidyr              1.1.4    2021-09-27 [1] CRAN (R 4.1.1)
   tidyselect         1.1.1    2021-04-30 [1] CRAN (R 4.1.1)
   tidytidbits        0.2.3    2021-03-08 [1] CRAN (R 4.1.1)
   tidytree           0.3.6    2021-11-12 [1] CRAN (R 4.1.1)
   treeio             1.18.1   2021-11-14 [1] Bioconductor
   tweenr             1.0.2    2021-03-23 [1] CRAN (R 4.1.1)
   tzdb               0.2.0    2021-10-27 [1] CRAN (R 4.1.1)
   usethis            2.1.5    2021-12-09 [1] CRAN (R 4.1.1)
   utf8               1.2.2    2021-07-24 [1] CRAN (R 4.1.1)
   vctrs              0.3.8    2021-04-29 [1] CRAN (R 4.1.1)
   viridis            0.6.2    2021-10-13 [1] CRAN (R 4.1.1)
   viridisLite        0.4.0    2021-04-13 [1] CRAN (R 4.1.1)
   vroom              1.5.7    2021-11-30 [1] CRAN (R 4.1.1)
   withr              2.4.3    2021-11-30 [1] CRAN (R 4.1.1)
   xfun               0.29     2021-12-14 [1] CRAN (R 4.1.1)
   xml2               1.3.3    2021-11-30 [1] CRAN (R 4.1.1)
   xtable             1.8-4    2019-04-21 [1] CRAN (R 4.1.1)
   XVector            0.34.0   2021-10-26 [1] Bioconductor
   yaml               2.2.1    2020-02-01 [1] CRAN (R 4.1.1)
   yulab.utils        0.0.4    2021-10-09 [1] CRAN (R 4.1.1)
   zlibbioc           1.40.0   2021-10-26 [1] Bioconductor
   zoo                1.8-9    2021-03-09 [1] CRAN (R 4.1.1)

 [1] /home/suwu/R/x86_64-pc-linux-gnu-library/4.1
 [2] /usr/local/lib/R/site-library
 [3] /usr/lib/R/site-library
 [4] /usr/lib/R/library

────────────────────────────────────────────────────────────────────────────────
```

## Tutorials

Open the `Analysis.R` from `Script` folder from our GitHub repository.
This script contains the command lines to execute all analyses.
