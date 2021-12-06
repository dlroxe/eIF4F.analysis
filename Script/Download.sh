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
