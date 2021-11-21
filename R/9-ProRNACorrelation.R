## prepare RNA proteomics datasets from CCLE and CPTAC LUAD ====================

CCLE_RNAseq <- CCLE_Anno <- CCLE_Proteomics <- CPTAC_LUAD_Proteomics <- CPTAC_LUAD_RNAseq <- NULL
#' Read all proteomics related datasets from CCLE and CPTAC LUAD
#'
#' @description This function reads all proteomics related datasets from CCLE and CPTAC LUAD. It generates five global variables.
#'
#' CCLE_RNAseq: the RNAseq data of CCLE from \code{CCLE_expression_full.csv}
#'
#' CCLE_Anno: the annotation data of CCLE from \code{sample_info.csv}
#'
#' CCLE_Proteomics: the proteomics data of CCLE from \code{protein_quant_current_normalized.csv}
#'
#' CPTAC_LUAD_Proteomics: the proteomics data with annotation of CPTAC LUAD from \code{Protein.xlsx}
#'
#' CPTAC_LUAD_RNAseq: the RNAseq data of CPTAC LUAD from \code{RNA.xlsx}
#'
#' @importFrom readxl read_excel
#'
#' @export
#'
#' @examples \dontrun{initialize_proteomics_data()}
#'
initialize_proteomics_data <- function() {
  CCLE_RNAseq <<- fread(
    file.path(
      data.file.directory,
      "CCLE_expression_full.csv"
    ),
    data.table = FALSE
  )

  CCLE_Anno <<- fread(
    file.path(
      data.file.directory,
      "sample_info.csv"
    ),
    data.table = FALSE
  ) %>% dplyr::select(1, 2)

  CCLE_Proteomics <<- fread(
    file.path(
      data.file.directory,
      "protein_quant_current_normalized.csv"
    ),
    data.table = FALSE
  )

  CPTAC_LUAD_Proteomics <<- read_excel(
    file.path(data.file.directory, "Protein.xlsx"),
    col_names = FALSE
  ) %>%
    # as.data.frame(.) %>%
    mutate(...1 = make.unique(.data$...1)) %>% # relabel the duplicates
    tibble::column_to_rownames(var = "...1") %>%
    t() %>%
    as_tibble() %>%
    mutate_at(vars(-.data$Type, -.data$Sample), ~ suppressWarnings(as.numeric(.))) # exclude two columns convert character to number

  CPTAC_LUAD_RNAseq <<- read_excel(
    file.path(data.file.directory, "RNA.xlsx"),
    col_names = FALSE
  ) %>%
    # as_tibble(.) %>%
    mutate(...1 = make.unique(.data$...1)) %>% # relabel the duplicates
    tibble::column_to_rownames(var = "...1") %>%
    t() %>%
    as_tibble() %>%
    mutate_at(vars(-.data$Type, -.data$Sample), ~ suppressWarnings(as.numeric(.))) # exclude two columns convert character to number
}


## Select EIF RNA/pro data and plotting ================================

#' Select subset of CCLE RNAseq data
#' @description This function selected the CCLE RNAseq data from the input \code{EIF}.
#'
#' @details The function should not be used directly, only inside \code{\link{plot_scatter_RNApro_CCLE}} function.
#' @param EIF gene name
#' @return a data frame of CCLE RNAseq data from the input \code{EIF} genes
#' @importFrom dplyr starts_with
#' @examples \dontrun{.get_CCLE_RNAseq_subset()}
#' @keywords internal
.get_CCLE_RNAseq_subset <- function(EIF) {
  .CCLE_RNAseq_subset <- CCLE_RNAseq %>%
    rename("DepMap_ID" = "V1") %>%
    dplyr::select("DepMap_ID", starts_with((!!paste(EIF, "(ENSG"))))
  return(.CCLE_RNAseq_subset)
}

#' Select subset of CCLE proteomics data
#' @description This function selected the CCLE proteomics data from the input \code{EIF}.
#'
#' @details The function should not be used directly, only inside \code{\link{plot_scatter_RNApro_CCLE}} function.
#' @param EIF protein name
#' @return a data frame of CCLE proteomics data from the input \code{EIF} genes
#' @importFrom dplyr across
#' @importFrom tidyselect contains vars_select_helpers
#' @examples \dontrun{.get_CCLE_Proteomics_subset()}
#' @keywords internal
.get_CCLE_Proteomics_subset <- function(EIF) {
  .CCLE_Proteomics_subset <- CCLE_Proteomics %>%
    dplyr::filter(.data$Gene_Symbol == EIF)
  if (nrow(.CCLE_Proteomics_subset) > 1) {
    df <- .CCLE_Proteomics_subset %>%
      dplyr::select(contains("_Peptides")) %>%
      mutate(sum = rowSums(dplyr::across(tidyselect::vars_select_helpers$where(is.numeric))))
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
    rownames_to_column(var = "rn") %>%
    mutate(
      Celline = sub("\\_.*", "", .data$rn),
      Type = sub(".*_ *(.*?) *_.*", "\\1", .data$rn)
    ) %>%
    select(-.data$rn) %>%
    mutate_if(is.character, as.factor) %>%
    # select(!!paste0(EIF,".pro") := EIF)
    rename(!!paste0(EIF, ".pro") := EIF) # rename with dynamic variables
  return(.CCLE_Proteomics_subset)
}

#' Scatter plots of correlation between RNAseq and proteomics data
#' @description This function should not be used directly,
#' only inside \code{\link{plot_scatter_RNApro_CCLE}} or \code{\link{plot_scatter_RNApro_LUAD}} function.
#' @param df dataframe of both RNAseq and proteomics generated inside
#' @param x protein or gene name
#' @param y database name. \code{CCLE} or \code{LUAD}
#' @importFrom ggpubr ggscatter
#' @keywords internal
.RNApro_scatterplot <- function(df, x, y) {
  p1 <- ggscatter(df,
    x = paste0(x, ".pro"),
    y = paste0(x, ".rna"), # color = "Type",
    add = "reg.line", # conf.int = TRUE,
    cor.coef = TRUE,
    cor.method = "pearson",
    title = paste(x, "(", y, ")"),
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
    path = file.path(output.directory, "LUAD"),
    filename = paste(y, x, "cor", ".pdf"),
    plot = p1,
    width = 3,
    height = 3,
    useDingbats = FALSE
  )
}


## master function to perform RNA protein correlation ==========================

#' Correlation between CCLE RNAseq and proteomics data
#' @description generates correlation scatter plot for eIF4F RNAseq and proteomics data from CCLE
#' @param EIF gene name
#' @details  This function merge the CCLE RNAseq values from EIF4F genes
#' in the data frame prepared from \code{\link{.get_CCLE_RNAseq_subset}} and
#' proteomics data of the same protein in the data frame prepared from \code{\link{.get_CCLE_Proteomics_subset}}.
#'
#' Then it uses the combined data to calculate the correlation coefficients
#' between protein and RNA levels, and plot the result with the function \code{\link{.RNApro_scatterplot}}
#' @return the correlation scatter plot for \code{EIF} RNA and protein values in CCLE
#' @export
#' @examples \dontrun{
#' plot_scatter_RNApro_CCLE("EIF4G1")
#' lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "PABPC1"), plot_scatter_RNApro_CCLE)
#' }
plot_scatter_RNApro_CCLE <- function(EIF) {
  .df <- .get_CCLE_RNAseq_subset(EIF) %>%
    full_join(CCLE_Anno, by = "DepMap_ID") %>%
    na.omit() %>%
    dplyr::select(-.data$DepMap_ID) %>%
    rename(
      "Celline" = "stripped_cell_line_name",
      !!paste0(EIF, ".rna") := contains(EIF)
    ) %>%
    full_join(.get_CCLE_Proteomics_subset(EIF), by = "Celline") %>%
    na.omit()
  .RNApro_scatterplot(df = .df, x = EIF, y = "CCLE")
}

#' Correlation between LUAD RNAseq and proteomics data
#' @description generates correlation scatter plot for eIF4F RNAseq and proteomics data from LUAD
#' @param EIF gene name
#' @details  This function merge the LUAD RNAseq values from EIF4F genes
#' in the data frame prepred from \code{CPTAC_LUAD_Proteomics} and
#' proteomics data of the same protein in the data frame prepared from \code{CPTAC_LUAD_RNAseq}.
#'
#' Then it uses the combined data to calculate the correlation coefficients
#' between protein and RNA levels, and plot the result with the function \code{\link{.RNApro_scatterplot}}
#' @return the correlation scatter plot for \code{EIF} RNA and protein values in LUAD
#' @export
#' @examples \dontrun{
#' plot_scatter_RNApro_LUAD("EIF4G1")
#' lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "PABPC1"), plot_scatter_RNApro_LUAD)
#' }
plot_scatter_RNApro_LUAD <- function(EIF) {
  .EIF_LUAD_Proteomics <- CPTAC_LUAD_Proteomics %>%
    dplyr::select(all_of(EIF), "Type", "Sample") %>%
    dplyr::filter(.data$Type == "Tumor")

  .EIF_LUAD_RNAseq <- CPTAC_LUAD_RNAseq %>%
    dplyr::select(all_of(EIF), "Type", "Sample") %>%
    dplyr::filter(.data$Type == "Tumor")

  df <- merge(.EIF_LUAD_Proteomics,
    .EIF_LUAD_RNAseq,
    by = c("Sample", "Type"),
    suffixes = c(".pro", ".rna")
  )

  .RNApro_scatterplot(df = df, x = EIF, y = "LUAD")
}

