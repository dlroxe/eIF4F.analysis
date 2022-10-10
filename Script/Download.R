# install.packages("RCurl")
library(RCurl)
library(R.utils)
library(utils)

### create the directory to store all downloaded datasets
data_file_directory <-
  Sys.getenv(c("DATA_FILE_DIRECTORY"), unset = "~/eIF4F.analysis/eIF4F_data")
output_directory <-
  Sys.getenv(c("OUTPUT_DIRECTORY"), unset = "~/eIF4F.analysis/eIF4F_output")

dir.create(data_file_directory)
dir.create(output_directory)

### set up download
getOption('timeout')
options(timeout = 200000)

### TCGA and GTEX DATA
#### TCGA CNV dataset (thresholded) 1
if (!file.exists(file.path(data_file_directory, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"))) {
  utils::download.file("https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz",
    destfile = file.path(data_file_directory, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"),
    method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"))
}


#### TCGA CNV dataset 2
if (!file.exists(file.path(data_file_directory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes"))) {
  utils::download.file("https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz",
    destfile = file.path(data_file_directory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"),
    method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"))
}


#### TCGA CNV ratio dataset 3
if (!file.exists(file.path(data_file_directory, "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena"))) {
  utils::download.file("https://pancanatlas.xenahubs.net/download/broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz",
                       file.path(data_file_directory, "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena.gz"))
}


#### TCGA RNA-Seq dataset 4
if (!file.exists(file.path(data_file_directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"))) {
  utils::download.file("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz",
                       file.path(data_file_directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"))
}


#### TCGA sample type annotation 5
if (!file.exists(file.path(data_file_directory, "TCGA_phenotype_denseDataOnlyDownload.tsv"))) {
  utils::download.file("https://pancanatlas.xenahubs.net/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz",
                       file.path(data_file_directory, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "TCGA_phenotype_denseDataOnlyDownload.tsv.gz"))
}


#### TCGA OS data ## 6
if (!file.exists(file.path(data_file_directory, "Survival_SupplementalTable_S1_20171025_xena_sp"))) {
  utils::download.file("https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp",
                       file.path(data_file_directory, "Survival_SupplementalTable_S1_20171025_xena_sp"),
                       method = "libcurl"
  )
}

#### TCGA and GTEX RNA-Seq dataset 7
if (!file.exists(file.path(data_file_directory, "TcgaTargetGtex_RSEM_Hugo_norm_count"))) {
  utils::download.file("https://toil.xenahubs.net/download/TcgaTargetGtex_RSEM_Hugo_norm_count.gz",
                       file.path(data_file_directory, "TcgaTargetGtex_RSEM_Hugo_norm_count.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "TcgaTargetGtex_RSEM_Hugo_norm_count.gz"))
}


#### TCGA and GTEX sample type annotation 8
if (!file.exists(file.path(data_file_directory, "TcgaTargetGTEX_phenotype.txt"))) {
  utils::download.file("https://toil.xenahubs.net/download/TcgaTargetGTEX_phenotype.txt.gz",
                       file.path(data_file_directory, "TcgaTargetGTEX_phenotype.txt.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "TcgaTargetGTEX_phenotype.txt.gz"))
}

### CPTAC DATA

#### CPTAC LUAD RNA-Seq data (Gillette et al., 2020) 9
if (!file.exists(file.path(data_file_directory, "RNA.xlsx"))) {
  utils::download.file("https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/RNA.xlsx.gz",
                       file.path(data_file_directory, "RNA.xlsx.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "RNA.xlsx.gz"))
}

#### CPTAC LUAD Proteomics (Gillette et al., 2020) 10
if (!file.exists(file.path(data_file_directory, "Protein.xlsx"))) {
  utils::download.file("https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/Protein.xlsx.gz",
                       file.path(data_file_directory, "Protein.xlsx.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "Protein.xlsx.gz"))
}

#### CPTAC LUAD Phosproteomics (Gillette et al., 2020) 11
if (!file.exists(file.path(data_file_directory, "Phos.xlsx"))) {
  utils::download.file("https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/Phos.xlsx.gz",
                       file.path(data_file_directory, "Phos.xlsx.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "Phos.xlsx.gz"))
}

#### CPTAC LUAD Sample Annotation 12
# S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx from PDC Study ID: PDC000153
if (!file.exists(file.path(data_file_directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx"))) {
  utils::download.file("https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx.gz",
                       file.path(data_file_directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx.gz"))
}

#### CPTAC LUAD Clinical Data 13
# S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx from PDC Study ID: PDC000153
if (!file.exists(file.path(data_file_directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx"))) {
  utils::download.file("https://github.com/a3609640/EIF-analysis/raw/master/LUAD%20Data/S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx.gz",
                       file.path(data_file_directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx.gz"))
}

### CCLE DATA

#### CCLE RNA-Seq data from DepMap Public 20Q4 20Q3 14
if (!file.exists(file.path(data_file_directory, "CCLE_expression_full.csv"))) {
  utils::download.file("https://ndownloader.figshare.com/files/24613349",
                       file.path(data_file_directory, "CCLE_expression_full.csv"),
                       method = "libcurl"
  )
}

#### CCLE annotation data 15
if (!file.exists(file.path(data_file_directory, "sample_info.csv"))) {
  utils::download.file("https://ndownloader.figshare.com/files/24613394",
                       file.path(data_file_directory, "sample_info.csv"),
                       method = "libcurl"
  )
}

#### CCLE proteomics data 16
if (!file.exists(file.path(data_file_directory, "protein_quant_current_normalized.csv"))) {
  utils::download.file("https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
                       file.path(data_file_directory, "protein_quant_current_normalized.csv.gz"),
                       method = "libcurl"
  )
  R.utils::gunzip(file.path(data_file_directory, "protein_quant_current_normalized.csv.gz"))
}




