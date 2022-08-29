# Correlation between RNA and protein levels for EIF4F
# This R script contains four sections.
#
# (1) CCLE and LUAD proteomics data preparation
#
# (2) selection of RNA and protein expression data and plotting
#
# (3) composite functions to execute a pipeline of functions to select related
#  expression with supply of EIF4F gene names for correlation analysis and
#  plotting.
#
# (4) wrapper function to call all master functions with inputs
#

## Wrapper function for data initialization of RNA and proteomics datasets =====

#' @noRd
## due to NSE notes in R CMD check
CCLE_RNAseq <- CCLE_Anno <- CCLE_Proteomics <- CPTAC_LUAD_Proteomics <- NULL
CPTAC_LUAD_RNAseq <- NULL

#' @title Read all proteomics related datasets from CCLE and CPTAC LUAD
#'
#' @description
#'
#' A wrapper function reads all proteomics related datasets from CCLE and
#'  CPTAC LUAD.
#'
#' @details Side effects:
#'
#' (1) `CCLE_RNAseq`: the RNAseq data of CCLE from the download dataset
#'  `CCLE_expression_full.csv`.
#'  It is stored as `CCLE_RNAseq.csv` in
#'  `~/Documents/EIF_output/ProcessedData` folder.
#'
#' (2) `CCLE_Anno`: the annotation data of CCLE from `sample_info.csv`.
#'  It is stored as `CCLE_Anno.csv` in
#'  `~/Documents/EIF_output/ProcessedData` folder.
#'
#' (3) `CCLE_Proteomics`: the proteomics data of CCLE from
#'  `protein_quant_current_normalized.csv`.
#'  It is stored as `CCLE_Proteomics.csv` in
#'  `~/Documents/EIF_output/ProcessedData` folder.
#'
#' (4) `CPTAC_LUAD_Proteomics`: the proteomics data with annotation of CPTAC
#'  LUAD study from the download data file `Protein.xlsx`.
#'  It is stored as `CPTAC_LUAD_Proteomics.csv` in
#'  `~/Documents/EIF_output/ProcessedData` folder.
#'
#' (5) `CPTAC_LUAD_RNAseq`: the RNAseq data of CPTAC LUAD samples from
#'  the download data file `RNA.xlsx`.
#'  It is stored as `CPTAC_LUAD_RNAseq.csv` in
#'  `~/Documents/EIF_output/ProcessedData` folder.
#'
#' @family wrapper function for data initialization
#'
#' @importFrom readxl read_excel
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_proteomics_data()
#' }
#'
initialize_proteomics_data <- function() {
  assign("CCLE_RNAseq",
         fread(
           file.path(
             data_file_directory,
             "CCLE_expression_full.csv"
           ),
           data.table = FALSE
         ),
         envir = parent.env(environment()))

  readr::write_csv(CCLE_RNAseq, file.path(output_directory,
                                          "ProcessedData",
                                          "CCLE_RNAseq.csv"))

  assign("CCLE_Anno",
         fread(
           file.path(
             data_file_directory,
             "sample_info.csv"
           ),
           data.table = FALSE
         ) %>% dplyr::select(1, 2),
         envir = parent.env(environment()))

  readr::write_csv(CCLE_Anno, file.path(output_directory,
                                        "ProcessedData",
                                        "CCLE_Anno.csv"))

  assign("CCLE_Proteomics",
         fread(
           file.path(
             data_file_directory,
             "protein_quant_current_normalized.csv"
           ),
           data.table = FALSE
         ),
         envir = parent.env(environment()))

  readr::write_csv(CCLE_Proteomics, file.path(output_directory,
                                              "ProcessedData",
                                              "CCLE_Proteomics.csv"))

  assign("CPTAC_LUAD_Proteomics",
         readxl::read_excel(
           file.path(data_file_directory, "Protein.xlsx"),
           col_names = FALSE
         ) %>%
           # as.data.frame(.) %>%
           dplyr::mutate(...1 = make.unique(.data$...1)) %>% # relabel the duplicates
           tibble::column_to_rownames(var = "...1") %>%
           t() %>%
           tibble::as_tibble() %>%
           dplyr::mutate_at(dplyr::vars(-.data$Type, -.data$Sample),
                            # exclude two columns convert character to number
                            ~ suppressWarnings(as.numeric(.))),
         envir = parent.env(environment()))

  readr::write_csv(CPTAC_LUAD_Proteomics, file.path(output_directory,
                                                    "ProcessedData",
                                                    "CPTAC_LUAD_Proteomics.csv"))

  assign("CPTAC_LUAD_RNAseq",
         readxl::read_excel(
           file.path(data_file_directory, "RNA.xlsx"),
           col_names = FALSE
         ) %>%
           # as_tibble(.) %>%
           dplyr::mutate(...1 = make.unique(.data$...1)) %>% # relabel the duplicates
           tibble::column_to_rownames(var = "...1") %>%
           t() %>%
           tibble::as_tibble() %>%
           dplyr::mutate_at(dplyr::vars(-.data$Type, -.data$Sample),
                            # exclude two columns convert character to number
                            ~ suppressWarnings(as.numeric(.))),
         envir = parent.env(environment()))

  readr::write_csv(CPTAC_LUAD_RNAseq, file.path(output_directory,
                                                "ProcessedData",
                                                "CPTAC_LUAD_RNAseq.csv"))
}


## Select EIF RNA/pro data and plotting ================================

#' @title Select subset of CCLE RNAseq data
#'
#' @description An internal helper function selects the data of input gene
#'  `EIF` from the data frame `CCLE_RNAseq` - a global variable generated from
#'  [initialize_proteomics_data].
#'
#' @details The function should not be used directly, only inside
#'  [.plot_scatter_RNApro_CCLE()] function.
#'
#' @param gene_name gene name
#'
#' @return a data frame of CCLE RNAseq data from the input `gene_name` genes
#'
#' @family helper function for RNA protein correlation analysis
#'
#' @importFrom dplyr starts_with
#'
#' @examples \dontrun{
#' .get_CCLE_RNAseq_subset()
#' }
#'
#' @keywords internal
#'
.get_CCLE_RNAseq_subset <- function(gene_name) {
  return(CCLE_RNAseq %>%
    rename("DepMap_ID" = "V1") %>%
    dplyr::select("DepMap_ID", dplyr::starts_with((!!paste(gene_name, "(ENSG")))))
}

#' @title Select subset of CCLE proteomics data
#'
#' @description A helper function selects the data of input gene `gene_name`.
#' from the data frame `CCLE_Proteomics` - a global variable generated from
#'  [initialize_proteomics_data].
#'
#' @details The function should not be used directly, only inside
#'  [.plot_scatter_RNApro_CCLE()] function.
#'
#' @param gene_name gene name
#'
#' @return a data frame of CCLE proteomics data from the input `gene_name`
#'
#' @family helper function for RNA protein correlation analysis
#'
#' @importFrom dplyr across
#'
#' @importFrom tidyselect contains vars_select_helpers
#'
#' @examples \dontrun{
#' .get_CCLE_Proteomics_subset()
#' }
#'
#' @keywords internal
#'
.get_CCLE_Proteomics_subset <- function(gene_name) {
  .CCLE_Proteomics_subset <- CCLE_Proteomics %>%
    dplyr::filter(.data$Gene_Symbol == gene_name)
  if (nrow(.CCLE_Proteomics_subset) > 1) {
    df <- .CCLE_Proteomics_subset %>%
      dplyr::select(contains("_Peptides")) %>%
      dplyr::mutate(sum = rowSums(
        dplyr::across(tidyselect::vars_select_helpers$where(is.numeric))))
    name <- rownames(df[df$sum == max(df$sum), ]) # for the maximum value
    .CCLE_Proteomics_subset <- .CCLE_Proteomics_subset %>%
      dplyr::filter(row.names(.CCLE_Proteomics_subset) == name)
  } else {
    TRUE
  }

  .CCLE_Proteomics_subset <- .CCLE_Proteomics_subset %>%
    tibble::column_to_rownames(var = "Gene_Symbol") %>%
    dplyr::select(contains("TenPx"), -contains("_Peptides")) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "rn") %>%
    dplyr::mutate(
      Celline = sub("\\_.*", "", .data$rn),
      Type = sub(".*_ *(.*?) *_.*", "\\1", .data$rn)
    ) %>%
    dplyr::select(-.data$rn) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    # select(!!paste0(EIF,".pro") := EIF)
    rename(!!paste0(gene_name, ".pro") := gene_name) # rename with dynamic variables

  return(.CCLE_Proteomics_subset)
}

#' @title Scatter plots of correlation between RNAseq and proteomics data
#'
#' @description A helper function plots RNAseq and proteomics data in scatter
#'  plot and analyzes the correlation
#'
#' @details This function should not be used directly, only inside
#'  [.plot_scatter_RNApro_CCLE()] or [.plot_scatter_RNApro_LUAD()] function.
#'
#' Side effects:
#'
#' (1) scatter plots to show correlation between RNA and protein
#'  expressions
#'
#' @param df dataframe of both RNAseq and proteomics generated inside
#'
#' @param gene_name protein or gene name
#'
#' @param cohort database name. `CCLE` or `LUAD`
#'
#' @family helper function for RNA protein correlation analysis
#'
#' @importFrom ggpubr ggscatter
#'
#' @keywords internal
#'
.RNApro_scatterplot <- function(df, gene_name, cohort) {
  p1 <- ggscatter(df,
    x = paste0(gene_name, ".pro"),
    y = paste0(gene_name, ".rna"), # color = "Type",
    add = "reg.line", # conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    title = paste(gene_name, "(", cohort, ")"),
    xlab = "Protein expresion",
    ylab = "RNA expression"
  ) +
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = black_bold_12,
      axis.title.y = black_bold_12,
      axis.text.x = black_bold_12,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = black_bold_12,
      strip.background = element_rect(fill = "white")
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "RNApro"),
    filename = paste(cohort, gene_name, "cor", ".pdf"),
    plot = p1,
    #width = 3,
    #height = 3,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )

  return(NULL)
}


## Composite function to perform RNA protein correlation =======================

#' @title Correlation between CCLE RNAseq and proteomics data
#'
#' @description
#'
#' An internal composite function generates correlation scatter plot for eIF4F RNAseq and
#' proteomics data from CCLE.
#'
#' @details This function
#'
#' * merges the CCLE RNAseq values of EIF4F genes in the data frame prepared
#'  from [.get_CCLE_RNAseq_subset()], the annotation data of CCLE `CCLE_Anno`,
#'  and proteomics data of the same EIF4F protein in the data frame prepared
#'  from [.get_CCLE_Proteomics_subset()].
#'
#' * uses the combined data to calculate the correlation coefficients
#'  between protein and RNA levels, and plot the result with the function
#'  [.RNApro_scatterplot()]
#'
#' This function is not accessible to the user and will not show at the users'
#' workspace. It can only be called by the exported [EIF4F_RNA_pro_correlation()]
#' function.
#'
#' Side effects:
#'
#' (1) scatter plot to show correlation between RNA and protein of
#'  input `gene_name` expressions in CCLE samples
#'
#' @param gene_name gene name
#'
#' @family composite function to call RNA protein correlation analysis and
#'  plotting
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' lapply(c("EIF4G1", "EIF4A1", "EIF4E"), .plot_scatter_RNApro_CCLE)
#' }
#'
.plot_scatter_RNApro_CCLE <- function(gene_name) {
  .df <- .get_CCLE_RNAseq_subset(gene_name) %>%
    dplyr::full_join(CCLE_Anno, by = "DepMap_ID") %>%
    stats::na.omit() %>%
    dplyr::select(-.data$DepMap_ID) %>%
    dplyr::rename(
      "Celline" = "stripped_cell_line_name",
      !!paste0(gene_name, ".rna") := tidyselect::contains(gene_name)
    ) %>%
    dplyr::full_join(.get_CCLE_Proteomics_subset(gene_name), by = "Celline") %>%
    stats::na.omit()

  .RNApro_scatterplot(df = .df, gene_name = gene_name, cohort = "CCLE")

  return(NULL)
}

#' @title Correlation between LUAD RNAseq and proteomics data
#'
#' @description
#' A composite function generates correlation scatter plot for eIF4F RNAseq and
#' proteomics data from LUAD.
#'
#' @details This function
#'
#' * merges the LUAD RNAseq data of EIF4F genes prepared from
#'  `CPTAC_LUAD_Proteomics` and proteomics data of the EIF4F proteins
#'  prepared from `CPTAC_LUAD_RNAseq`.
#' * uses the combined data to calculate the correlation coefficients between
#'  protein and RNA levels, and plot the result with the function
#'  [.RNApro_scatterplot()]
#'
#' This function is not accessible to the user and will not show at the users'
#' workspace. It can only be called by the exported [EIF4F_RNA_pro_correlation()]
#' function.
#'
#' Side effects:
#'
#' (1) scatter plot to show correlation between RNA and protein
#'  expressions of input `gene_name` in CPTAC LUAD samples
#'
#' @param gene_name gene name
#'
#' @family composite function to call RNA protein correlation analysis and
#'  plotting
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' lapply(c("EIF4G1", "EIF4A1", "EIF4E"), .plot_scatter_RNApro_LUAD)
#' }
#'
.plot_scatter_RNApro_LUAD <- function(gene_name) {
  .EIF_LUAD_Proteomics <- CPTAC_LUAD_Proteomics %>%
    dplyr::select(dplyr::all_of(gene_name), "Type", "Sample") %>%
    dplyr::filter(.data$Type == "Tumor")

  .EIF_LUAD_RNAseq <- CPTAC_LUAD_RNAseq %>%
    dplyr::select(dplyr::all_of(gene_name), "Type", "Sample") %>%
    dplyr::filter(.data$Type == "Tumor")

  df <- merge(.EIF_LUAD_Proteomics,
    .EIF_LUAD_RNAseq,
    by = c("Sample", "Type"),
    suffixes = c(".pro", ".rna")
  )

  .RNApro_scatterplot(df = df, gene_name = gene_name, cohort = "LUAD")

  return(NULL)
}


## Wrapper function to call all composite functions with inputs ================

#' @title Perform RNA protein correlation and generate scatter plots
#'
#' @description
#'
#' A wrapper function to call all composite functions for RNA and protein
#'  correlation with inputs.
#'
#' @details
#'
#' This function run the internal composite functions [.plot_scatter_RNApro_CCLE()] and
#'  [.plot_scatter_RNApro_LUAD()] with EIF4F gene name as inputs.
#'
#' Side effects:
#'
#' (1) scatter plot to show correlation between EIF4F RNA and protein
#'  expressions in CCLE samples
#'
#' (2) scatter plot to show correlation between EIF4F RNA and protein
#'  expressions in CPTAC LUAD samples
#'
#' @family wrapper function to call all composite functions with inputs
#'
#' @export
#'
#' @examples \dontrun{
#' EIF4F_RNA_pro_correlation()
#' }
#'
EIF4F_RNA_pro_correlation <- function() {
  lapply(c("EIF4G1", "EIF4A1", "EIF4E"), .plot_scatter_RNApro_CCLE)
  lapply(c("EIF4G1", "EIF4A1", "EIF4E"), .plot_scatter_RNApro_LUAD)
}
