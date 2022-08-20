## Wrapper function for directory initialization ===============================
#' @noRd
data_file_directory <- "~/Downloads/EIF_data"
output_directory <- "~/Documents/EIF_output"

#' @title Set output directories
#'
#' @description A wrapper function sets up the output directories
#'
#' @family wrapper function for setting up output directories
#'
#' @details Side effects:
#'
#' (1) output directory for CNV analysis results
#'  * `~/Documents/EIF_output/CNV`
#'
#' (2) output directory for differential gene expression and ratio analysis
#'  results
#'  * `~/Documents/EIF_output/DEG`
#'
#' (3) output directories for survival analysis results
#'  * `~/Documents/EIF_output/Survival/KM`
#'  * `~/Documents/EIF_output/Survival/CoxPH`
#'
#' (4) output directory for PCA results
#'  * `~/Documents/EIF_output/PCA/TCGA_tumor+GTEX_healthy`
#'  * `~/Documents/EIF_output/PCA/TCGA_tumor`
#'  * `~/Documents/EIF_output/PCA/GTEX_healthy`
#'  * `~/Documents/EIF_output/PCA/matched_tumor_and_healthy`
#'  * `~/Documents/EIF_output/PCA/LUAD`
#'
#' (5) output directory for analysis results on EIF4F correlating genes
#'  * `~/Documents/EIF_output/CORs`
#'
#' (6) output directory for correlation analysis between RNA and protein levels
#'  * `~/Documents/EIF_output/RNApro`
#'
#' (7) output directory for co-expression analysis among EIF4F subunits and
#'  their differential expression
#'  * `~/Documents/EIF_output/Proteomics/COEXP`
#'  * `~/Documents/EIF_output/Proteomics/DIFFEXP`
#'
#' (8) output directory for processed datasets used as global variables
#'  their differential expression
#'  * `~/Documents/EIF_output/ProcessedData`
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_dir()
#' }
#'
initialize_dir <- function() {
  dir.create(file.path(output_directory, "CNV"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "DEG"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Survival", "KM"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Survival", "CoxPH"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "TCGA_tumor+GTEX_healthy"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "TCGA_tumor"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "GTEX_healthy"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "matched_tumor_and_healthy"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "LUAD"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "CORs"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "RNApro"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Proteomics","COEXP"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Proteomics","DIFFEXP"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "ProcessedData"),
             showWarnings = FALSE, recursive = TRUE)
}


## Wrapper function for data initialization ====================================
#' @title A wrapper function reads datasets from the download data files
#'
#' @description A wrapper function for data initialization loads the download
#' data files.
#'
#' @details
#' [initialize_data()] runs five data initialization functions
#' to load the relevant data files.
#' * [initialize_cnv_data()]
#' * [initialize_RNAseq_data()]
#' * [initialize_survival_data()]
#' * [initialize_proteomics_data()]
#' * [initialize_phosphoproteomics_data()]
#'
#' data initialization functions read datasets and generates global variable
#' available to functions in the package. They are not accessible to the user
#' and will not show at the users' workspace. This package has saved them as
#' csv files under `~/Documents/EIF_output/ProcessedData` in case that
#' users would like to access the data frames for their own analyses.
#'
#' Side effects:
#'
#' (1) `TCGA_CNV_value`imports the download dataset
#'  `Gistic2_CopyNumber_Gistic2_all_data_by_genes`. It is stored
#'  as `TCGA_CNV_value.csv`
#'
#' (2) `TCGA_CNV_sampletype` imports and merges datasets from
#'  `Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes` and
#'  `TCGA_phenotype_denseDataOnlyDownload.tsv`. It is stored
#'  as `TCGA_CNV_sampletype.csv`
#'
#' (3) `TCGA_CNVratio_sampletype` imports and merges datasets from
#'  `broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena` and
#'  `TCGA_phenotype_denseDataOnlyDownload.tsv`. It is stored
#'  as `TCGA_CNVratio_sampletype.csv`
#'
#' (4) `TCGA_GTEX_RNAseq_sampletype` imports and merges datasets from
#'  `TcgaTargetGtex_RSEM_Hugo_norm_count` and `TcgaTargetGTEX_phenotype.txt`.
#'  It is stored as `TCGA_GTEX_RNAseq_sampletype.csv`.
#'
#' (5) `TCGA_RNAseq_OS_sampletype` imports and merges datasets from
#'  `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`,
#'  `Survival_SupplementalTable_S1_20171025_xena_sp` and
#'  `TCGA_phenotype_denseDataOnlyDownload.tsv`.
#'  It is stored as `TCGA_GTEX_RNAseq_sampletype.csv`.
#'
#' (6) `CCLE_RNAseq` imports the download dataset
#'  `CCLE_expression_full.csv`. It is stored as `CCLE_RNAseq.csv`.
#'
#' (7) `CCLE_Anno`imports the download dataset `sample_info.csv`. It is stored
#'  as `CCLE_Anno.csv`.
#'
#' (8) `CCLE_Proteomics` imports the download dataset
#'  `protein_quant_current_normalized.csv`. It is stored as `CCLE_Proteomics.csv`.
#'
#' (9) `CPTAC_LUAD_Proteomics` imports the download dataset `Protein.xlsx`.
#'  It is stored as `CPTAC_LUAD_Proteomics.csv`.
#'
#' (10) `CPTAC_LUAD_RNAseq`imports the download dataset `RNA.xlsx`.
#'  It is stored as `CPTAC_LUAD_RNAseq.csv`.
#'
#' (11) `CPTAC_LUAD_Phos`imports the download dataset `Phos.xlsx`.
#'  It is stored as `CPTAC_LUAD_Phos.csv`.
#'
#' (12) `CPTAC_LUAD_Clinic_Sampletype` imports and merges datasets from
#'   `S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx` and
#'   `S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx`.
#'   It is stored as `CPTAC_LUAD_Clinic_Sampletype.csv`.
#'
#' @family wrapper function for data initialization
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_data()
#' }
#'
initialize_data <- function() {
  tryCatch({
    rlang::env_binding_unlock(parent.env(environment()), nms = NULL)
    initialize_cnv_data()
    initialize_RNAseq_data()
    initialize_survival_data()
    initialize_proteomics_data()
    initialize_phosphoproteomics_data()
  },
  finally = {
    rlang::env_binding_lock(parent.env(environment()), nms = NULL)
  })
}


## Wrapper function for format initialization ==================================
#' @noRd
## due to NSE notes in R CMD check
black_bold_tahoma_7 <- black_bold_12 <- black_bold_12_45 <- NULL
black_bold_16 <- black_bold_16_right <- NULL
black_bold_16_45 <- black_bold_16_90 <- black_bold_18 <- NULL
col_vector <- NULL

#' @title Set format for font and color of plots
#'
#' @description A wrapper function set up the font type, size and color for
#'  ggplots.
#'
#' @family wrapper function for format initialization
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_format()
#' }
#'
initialize_format <- function() {
  parent_env <- parent.env(environment())

  rlang::env_binding_unlock(parent_env, nms = NULL)

  assign("black_bold_tahoma_7",
         element_text(color = "black",
                      face = "bold",
                      size = 7),
         envir = parent_env)


  assign("black_bold_12",
         element_text(color = "black",
                      face = "bold",
                      size = 12),
         envir = parent_env)

  assign("black_bold_12_45",
         element_text(
           color = "black",
           face = "bold",
           size = 12,
           angle = 45,
           hjust = 1
         ),
         envir = parent_env)

  assign("black_bold_16",
         element_text(color = "black",
                      face = "bold",
                      size = 16),
         envir = parent_env)

  assign("black_bold_16_right",
         element_text(color = "black",
                      face = "bold",
                      size = 16,
                      angle = 90),
         envir = parent_env)

  assign("black_bold_16_45",
         element_text(
           color = "black",
           face = "bold",
           size = 16,
           angle = 45,
           hjust = 1
         ),
         envir = parent_env)

  assign("black_bold_16_90",
         element_text(
           color = "black",
           face = "bold",
           size = 16,
           angle = 90,
           hjust = 1,
           vjust = 0.5
         ),
         envir = parent_env)

  assign("black_bold_18",
         element_text(
           color = "black",
           face = "bold",
           size = 18
         ),
         envir = parent_env)

  assign("col_vector",
         color(),
         envir = parent_env)

  rlang::env_binding_lock(parent.env(environment()), nms = NULL)
}


#' @title Generates a number of most distinctive colors.
#'
#' @description A helper function generates a number of most distinctive colors
#'  for plotting
#'
#' @details This function should not be used directly,
#' only inside [initialize_format()] function.
#'
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#'
#' @return a vector with 74 different colors
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' color()
#' }
#'
color <- function() {
  qual_col_pals <- RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual", ]

  return(unlist(mapply(
    RColorBrewer::brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  )))
}


