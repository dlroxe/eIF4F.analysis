# CNV analyses of EIF genes in TCGA data ---------------------------------------

# This R script contains four sections:
# (1) CNV data preparation
# (2) CNV data analyses and plotting
# (3) master functions to execute a pipeline of functions to select related CNV data
# with supply of EIF4F gene names for analysis and plotting.
# (4) wrapper function to call all master functions with inputs


## prepare CNV related datasets from TCGA ======================================

# due to NSE notes in R CMD check
TCGA_CNV_value <- TCGA_CNV_sampletype <- TCGA_CNVratio_sampletype <- NULL

#' Read all CNV related datasets from TCGA
#'
#' @description This function reads all CNV related datasets from TCGA and generates three global variables.
#'
#' TCGA_CNV_value: the unthreshold CNV value data generated from \code{\link{.get_TCGA_CNV_value}}
#'
#' TCGA_CNV_sampletype: the merged dataset from .TCGA.CNV, the threshold CNV data generated from \code{\link{.get_TCGA_CNV}},
#' and .TCGA_sampletype, the annotation data from the "TCGA_phenotype_denseDataOnlyDownload.tsv"
#' dataset with selection of sample.type and primary.disease columns. “Solid Tissue Normal” samples are excluded.
#'
#' TCGA_CNVratio_sampletype: the merged dataset from .TCGA_CNV_ratio, the CNV ratio data generated from \code{\link{.get_TCGA_CNV_ratio}},
#' and .TCGA_sampletype.
#'
#' @importFrom data.table fread transpose
#' @importFrom dplyr distinct filter select select_if mutate mutate_at summarise rename group_by
#' @importFrom stats na.omit
#' @importFrom tibble remove_rownames column_to_rownames
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_cnv_data()
#' }
#'
initialize_cnv_data <- function() {
  # .TCGA_CNV <- .TCGA_CNV_ratio <- .TCGA_sampletype <- NULL
  TCGA_CNV_value <<- .get_TCGA_CNV_value()

  .TCGA_CNV <- .get_TCGA_CNV()

  .TCGA_CNV_ratio <- .get_TCGA_CNV_ratio()

  .TCGA_sampletype <- readr::read_tsv(file.path(
    data_file_directory,
    "TCGA_phenotype_denseDataOnlyDownload.tsv"
  )) %>%
    as.data.frame() %>%
    dplyr::distinct(.data$sample, .keep_all = TRUE) %>%
    na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "sample") %>%
    dplyr::select("sample_type", "_primary_disease") %>%
    dplyr::rename(
      "sample.type" = "sample_type",
      "primary.disease" = "_primary_disease"
    )

  TCGA_CNV_sampletype <<- merge(.TCGA_CNV,
    .TCGA_sampletype,
    by    = "row.names",
    all.x = TRUE
  ) %>%
    dplyr::filter(.data$sample.type != "Solid Tissue Normal") %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Row.names")

  TCGA_CNVratio_sampletype <<- merge(.TCGA_CNV_ratio,
    .TCGA_sampletype,
    by    = "row.names",
    all.x = TRUE
  ) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Row.names")
}

#' Read CNV threshold dataset from TCGA
#' @description This function reads the threshold CNV data "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
#' which groups the CNV values in tumors into five categories.
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize_cnv_data}} function.
#' @importFrom tibble as_tibble
#' @return a data frame that contains CNV threshold data for all TCGA tumor
#' @examples \dontrun{
#' .get_TCGA_CNV()
#' }
#' @keywords internal
.get_TCGA_CNV <- function() {
  .TCGA_pancancer <- data.table::fread(
    file.path(
      data_file_directory,
      "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
    ),
    data.table = FALSE
  ) %>%
    tibble::as_tibble() %>%
    # as.data.frame(.) %>%
    dplyr::distinct(.data$Sample, .keep_all = TRUE) %>%
    na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Sample")

  # transpose function from the data.table library keeps numeric values as numeric.
  .TCGA_pancancer_transpose <- data.table::transpose(.TCGA_pancancer)
  # get row and column names in order
  rownames(.TCGA_pancancer_transpose) <- colnames(.TCGA_pancancer)
  colnames(.TCGA_pancancer_transpose) <- rownames(.TCGA_pancancer)
  return(.TCGA_pancancer_transpose)
}

#' Read CNV value dataset from TCGA
#' @description This function reads the unthreshold CNV value data from TCGA "Gistic2_CopyNumber_Gistic2_all_data_by_genes".
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize_cnv_data}} function.
#' @return a data frame that contains CNV value data for each tumor
#' @examples \dontrun{
#' .get_TCGA_CNV_value()
#' }
#' @keywords internal
.get_TCGA_CNV_value <- function() {
  .TCGA_pancancer <- fread(
    file.path(
      data_file_directory,
      "Gistic2_CopyNumber_Gistic2_all_data_by_genes"
    ),
    data.table = FALSE
  ) %>%
    as_tibble() %>%
    # as.data.frame(.) %>%
    dplyr::distinct(.data$Sample, .keep_all = TRUE) %>%
    na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Sample")

  # transpose function from the data.table library keeps numeric values as numeric.
  .TCGA_pancancer_transpose <- data.table::transpose(.TCGA_pancancer)
  # get row and colnames in order
  rownames(.TCGA_pancancer_transpose) <- colnames(.TCGA_pancancer)
  colnames(.TCGA_pancancer_transpose) <- rownames(.TCGA_pancancer)
  return(.TCGA_pancancer_transpose)
}

#' Read CNV ratio dataset from TCGA
#' @description This function reads the CNV ratio data "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena",
#' which calculates the ratios of CNV value in each tumor to the average CNV value from normals in the same cancer type.
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize_cnv_data}} function.
#' @return a data frame that contains CNV ratio data for each tumor
#' @examples \dontrun{
#' .get_TCGA_CNV_ratio()
#' }
.get_TCGA_CNV_ratio <- function() {
  .TCGA_pancancer <- fread(
    file.path(
      data_file_directory,
      "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena"
    ),
    data.table = FALSE
  ) %>%
    as_tibble() %>%
    # as.data.frame(.) %>%
    dplyr::distinct(.data$sample, .keep_all = TRUE) %>%
    na.omit() %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "sample")

  # transpose function from the data.table library keeps numeric values as numeric.
  .TCGA_pancancer_transpose <- data.table::transpose(.TCGA_pancancer)
  # get row and colnames in order
  rownames(.TCGA_pancancer_transpose) <- colnames(.TCGA_pancancer)
  colnames(.TCGA_pancancer_transpose) <- rownames(.TCGA_pancancer)
  return(.TCGA_pancancer_transpose)
}

## CNV data analysis and plotting ==============================================

#' Calculates the frequency of CNV status in all TCGA cancer types combined
#' @description This function calculates the frequency of each CNV status across tumors in all TCGA cancer types for every EIF4F gene.
#' @details It should not be used directly, only inside \code{\link{.plot_bargraph_CNV_TCGA}} function.
#' @param df \code{.TCGA_CNV_sampletype_subset} generated inside \code{\link{.plot_bargraph_CNV_TCGA}}
#' @return a summary table ranking the EIF4F gene by the frequencies of duplication (CNV threshold labeled as “1” in the dataset)
#' @importFrom reshape2 dcast melt
#' @examples \dontrun{
#' .CNV_all_cancer(.TCGA_CNV_sampletype_subset)
#' }
#' @keywords internal
.CNV_all_cancer <- function(df) {
  .TCGA_CNV_anno_subset_long <- melt(df,
    id = c(
      "sample.type",
      "primary.disease"
    ),
    value.name = "CNV"
  ) %>%
    mutate_if(is.character, as.factor)

  .CNV_sum <-
    table(.TCGA_CNV_anno_subset_long[, c("CNV", "variable")]) %>%
    # tibble(.) %>%
    as.data.frame() %>%
    mutate(CNV = factor(.data$CNV, levels = c("-2", "-1", "0", "1", "2")))

  # reorder stack bars by the frequency of duplication.
  Freq.sum <- dcast(.CNV_sum, variable ~ CNV, mean)
  .CNV_sum$variable <- factor(.CNV_sum$variable,
    levels = Freq.sum[order(Freq.sum$`1`), ]$variable
  )
  return(.CNV_sum)
}

#' Stacked bar plots of the CNV status summary
#' @description This function generates the stacked bar plots to the summary of CNV statuses
#' @details This plot function uses dataset generated from \code{\link{.CNV_all_cancer}} function
#'
#' It should not be used directly, only inside \code{\link{.plot_bargraph_CNV_TCGA}}function.
#' @param data summary table from \code{.CNV_all_cancer(.TCGA_CNV_sampletype_subset)}
#' @return stacked bar plots showing the summary table of the CNV statuses
#' @examples \dontrun{
#' .CNV_sum_barplot(.CNV_all_cancer(.TCGA_CNV_sampletype_subset))
#' }
#' @keywords internal
.CNV_sum_barplot <- function(data) {
  p1 <- ggplot2::ggplot(
    data,
    aes_(
      fill = ~CNV,
      y = ~Freq,
      x = ~variable
    )
  ) +
    geom_bar(
      stat = "identity",
      position = "fill"
    ) +
    geom_col() +
    geom_text(aes(label = paste0(.data$Freq / 100, "%")),
      position = position_stack(vjust = 0.5),
      size = 4
    ) +
    # scale_y_continuous(labels = scales::percent_format())+
    labs(
      x = "Tumor types (TCGA pan cancer atlas 2018)",
      y = "All TCGA tumors combined"
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title = black_bold_16,
      axis.title.x = black_bold_16,
      axis.title.y = element_blank(),
      axis.text.x = black_bold_16,
      axis.text.y = black_bold_16,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_16,
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_16
    ) +
    guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
    scale_fill_manual(
      name = "Copy number variation",
      breaks = c("-2", "-1", "0", "1", "2"),
      labels = c(
        "Deep del\n 0",
        "Shallow del\n 1",
        "Diploid\n 2",
        "Gain\n 3",
        "Amp\n 3+"
      ),
      values = c("darkblue", "blue", "lightgreen", "red", "darkred")
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "CNV"),
    filename = "EIFCNVsum.pdf",
    plot = p1,
    width = 9,
    height = 9,
    useDingbats = FALSE
  )
}

#' Calculates the frequency of CNV status in individual TCGA cancer types
#' @description a data analysis function that calculates the frequency of CNV status for each EIF4F gene in individual TCGA cancer types.
#' @details It should not be used directly, only inside \code{\link{.plot_bargraph_CNV_TCGA}} function.
#' @param df \code{.TCGA_CNV_sampletype_subset} generated inside \code{\link{.plot_bargraph_CNV_TCGA}}
#' @param x one gene from the input argument of \code{\link{.plot_bargraph_CNV_TCGA}}
#' @return a list with the summary table of CNV in individual TCGA cancer types and gene name
#' @importFrom forcats fct_rev
#' @importFrom tidyselect all_of any_of
#' @examples \dontrun{
#' lapply(EIF, .CNV_ind_cancer, df = .TCGA_CNV_sampletype_subset)
#' }
#' @keywords internal
.CNV_ind_cancer <- function(df, x) {
  .TCGA_CNV_anno_subset_long <- df %>%
    dplyr::select(
      all_of(x),
      "sample.type",
      "primary.disease"
    ) %>%
    melt(
      id = c(
        "sample.type",
        "primary.disease"
      ),
      value.name = "CNV"
    ) %>%
    mutate_if(is.character, as.factor)

  .CNV_sum <-
    table(.TCGA_CNV_anno_subset_long[, c("CNV", "primary.disease")]) %>%
    # tibble(.) %>%
    as.data.frame() %>%
    mutate(CNV = factor(.data$CNV, levels = c("-2", "-1", "0", "1", "2"))) %>%
    mutate(primary.disease = forcats::fct_rev(.data$primary.disease))

  output <- list(.CNV_sum, x)
  return(output)
}

#' Stacked bar plots of the CNV status
#' @description This function generates stacked bar plots for CNV status of each gene in an individual cancer type.
#' @details This plot function uses dataset generated from \code{\link{.CNV_ind_cancer}} function
#' It should not be used directly, only inside \code{\link{.plot_bargraph_CNV_TCGA}} function.
#' @param df \code{.EIF_CNV_ind_cancer} generated inside \code{\link{.plot_bargraph_CNV_TCGA}}
#' @return stacked bar plots for CNV status of each gene in an individual cancer type.
#' @examples \dontrun{
#' lapply(.EIF_CNV_ind_cancer, .CNV_barplot)
#' }
#' @keywords internal
.CNV_barplot <- function(df) {
  p1 <- ggplot(
    df[[1]],
    aes_(
      fill = ~CNV,
      order = ~ as.numeric(CNV),
      y = ~Freq,
      x = ~primary.disease
    )
  ) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      x = "Tumor types (TCGA pan cancer atlas 2018)",
      y = paste0("Percentages of ", df[[2]], " CNVs")
    ) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = black_bold_12,
      axis.title.y = element_blank(),
      axis.text.x = black_bold_12,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_12,
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_12
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
    scale_fill_manual(
      name = "Copy number variation",
      breaks = c("-2", "-1", "0", "1", "2"),
      labels = c(
        "Deep del\n 0",
        "Shallow del\n 1",
        "Diploid\n 2",
        "Gain\n 3",
        "Amp\n 3+"
      ),
      values = c(
        "darkblue", "blue",
        "lightgreen", "red",
        "darkred"
      )
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "CNV"),
    filename = paste0(df[[2]], "pancancerCNV.pdf"),
    plot = p1,
    width = 7.5,
    height = 9,
    useDingbats = FALSE
  )
}

#' Correlation coefficients for CNV values
#' @description This function calculates the correlation coefficients between every two genes and plot the correlation matrix with the function.
#' @details This plot function uses dataset \code{TCGA_CNV_value} generated from \code{\link{initialize_cnv_data}} function
#' It should not be used directly, only inside \code{\link{.plot_matrix_CNVcorr_TCGA}} function.
#' @param df output of \code{TCGA_CNV_value \%>\% select(all_of(EIF))} generated inside \code{\link{.plot_matrix_CNVcorr_TCGA}}
#' @return stacked bar plots for CNV status of each gene in an individual cancer type.
#' @importFrom corrplot cor.mtest corrplot
#' @importFrom grDevices dev.off pdf
#' @importFrom stats cor cor.test setNames
#' @examples \dontrun{
#' TCGA_CNV_value %>%
#'   dplyr::select(all_of(EIF_list)) %>%
#'   .matrix_plot()
#' }
#' @keywords internal
.matrix_plot <- function(df) {
  # correlation plot
  M <- stats::cor(df)
  testRes <- corrplot::cor.mtest(df, conf.level = 0.95)

  p1 <- corrplot::corrplot(
    M,
    method      = "color",
    cl.pos      = "n", # remove color legend
    tl.cex      = 1,
    number.cex  = 1,
    addgrid.col = "gray",
    addCoef.col = "black",
    tl.col      = "black",
    type        = "lower",
    order       = "FPC",
    tl.srt      = 45,
    p.mat       = testRes$p,
    sig.level   = 0.05, # insig = "blank"
  )
  print(p1) # print correlation matrix on the screen
  # save correlation plot as a pdf file
  pdf(
    file.path(output_directory, "CNV", "EIFCNVcormatrix.pdf"),
    width = 9,
    height = 9,
    useDingbats = FALSE
  )
  corrplot::corrplot(
    M,
    method      = "color",
    cl.pos      = "n", # remove color legend
    tl.cex      = 1,
    number.cex  = 1,
    addgrid.col = "gray",
    addCoef.col = "black",
    tl.col      = "black",
    type        = "lower",
    order       = "FPC",
    tl.srt = 45,
    p.mat       = testRes$p,
    sig.level   = 0.05, # insig = "blank"
  )
  dev.off()
}

#' Calculates the frequency of CNV status in all TCGA cancer types combined
#' @description This function selects the CNV ratio data in tumors vs adjacent normals from individual TCGA cancer types for each input gene.
#' @details It should not be used directly, only inside \code{\link{.plot_boxgraph_CNVratio_TCGA}} function.
#' @param df \code{.TCGA_CNVratio_sampletype_subset} generated inside \code{\link{.plot_boxgraph_CNVratio_TCGA}}
#' @param x one gene from the input argument of \code{\link{.plot_boxgraph_CNVratio_TCGA}}
#' @return a list with the data frame of CNV ratio of the input gene in individual TCGA cancer types and gene name
#' @examples \dontrun{
#' lapply(EIF, .CNVratio_tumor, df = .TCGA_CNVratio_sampletype_subset)
#' }
#' @keywords internal
.CNVratio_tumor <- function(df, x) {
  .CNVratio_data <- df %>%
    dplyr::select(all_of(x), "sample.type", "primary.disease") %>%
    melt(
      id = c("sample.type", "primary.disease"),
      value.name = "CNV"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(primary.disease = forcats::fct_rev(.data$primary.disease))
  output <- list(.CNVratio_data, x)
  return(output)
}

#' Box plots of the CNV ratios in tumors vs adjacent normals
#' @description This function generates boxplot for CNV ratios in tumors vs adjacent normals from individual TCGA cancer types.
#' @details This plot function uses dataset \code{.TCGA_CNVratio_sampletype_subset} generated from \code{\link{.plot_boxgraph_CNVratio_TCGA}} function.
#'
#' It should not be used directly, only inside \code{\link{.plot_boxgraph_CNVratio_TCGA}} function.
#' @param df \code{.EIF_CNVratio_ind_cancer} generated inside \code{\link{.plot_boxgraph_CNVratio_TCGA}}
#' @return boxplot for CNV ratios in tumors vs adjacent normals from individual TCGA cancer types.
#' @examples \dontrun{
#' lapply(.EIF_CNVratio_ind_cancer, .CNVratio_boxplot)
#' }
#' @keywords internal
.CNVratio_boxplot <- function(df) {
  p1 <- ggplot(
    data = df[[1]],
    aes_(
      y = ~ 2**CNV,
      x = ~primary.disease,
      # x = f.ordered1,
      color = ~primary.disease
    )
  ) +
    ylim(0, 3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    stat_n_text(
      size = 5,
      fontface = "bold",
      hjust = 0
    ) +
    geom_boxplot(
      alpha = .01,
      outlier.colour = NA,
      # size     = .75,
      # width    = 1,
      position = position_dodge(width = .9)
    ) +
    labs(
      x = "primary disease",
      y = paste(df[[2]], "CNV ratio", "(tumor/normal)")
    ) +
    # scale_color_manual(values = col_vector) +
    coord_flip() +
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = black_bold_12,
      axis.title.y = element_blank(),
      axis.text.x = black_bold_12,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_12,
      legend.position = "none",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_12
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "CNV"),
    filename = paste0(df[[2]], "pancancerCNVratio.pdf"),
    plot = p1,
    width = 7,
    height = 9,
    useDingbats = FALSE
  )
}


## master function to call CNV data analysis and plotting ======================

#' Summary of CNV statuses in bar plots
#' @description Provides the summary of CNV statuses of EIF4F genes in tumors from all TCGA
#' cancer types combined and in tumors from individual TCGA cancer types
#' @param EIF_list gene names in a vector of characters
#' @details  This function first selects CNV and sample type data from only EIF4F gene
#' in the data frame \code{TCGA_CNV_sampletype} prepared from \code{\link{initialize_cnv_data}}.
#'
#' Then it uses the subset data \code{.TCGA_CNV_sampletype_subset} to perform the CNV status analysis of all inquired genes
#' with the functions \code{\link{.CNV_all_cancer}} and plot the results as a bar plot with \code{\link{.CNV_sum_barplot}}
#'
#' It also uses the same subset data to perform CNV status analysis for each gene across
#' all tumors types with the function \code{\link{.CNV_ind_cancer}} and plots the CNV results and \code{\link{.CNV_barplot}}
#'
#' It should not be used directly, only inside \code{\link{EIF4F_CNV_analysis}} function.
#' @return stacked bar plots for grouped CNV status of \code{EIF_list} in TCGA tumors
#' @examples \dontrun{
#' plot_bargraph_CNV_TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))
#' }
#' @examples \dontrun{
#' plot_bargraph_CNV_TCGA(c(
#'   "TP53", "EIF4A1", "EIF4A2", "EIF4E",
#'   "EIF4E2", "EIF4E3", "MYC", "EIF3D", "EIF4EBP1", "EIF4G1", "EIF4G2", "EIF4G3",
#'   "PABPC1", "MKNK1", "MKNK2"
#' ))
#' }
#' @keywords internal
.plot_bargraph_CNV_TCGA <- function(EIF_list) {
  .TCGA_CNV_sampletype_subset <- NULL
  # combine CNV data with annotation data
  .TCGA_CNV_sampletype_subset <- TCGA_CNV_sampletype %>%
    dplyr::select(all_of(EIF_list), "sample.type", "primary.disease")

  # stacked bar plots for CNV status in combined TCGA tumors
  .CNV_all_cancer(.TCGA_CNV_sampletype_subset) %>%
    .CNV_sum_barplot()

  # stacked bar plots for CNV status in each TCGA tumor type
  .EIF_CNV_ind_cancer <- lapply(EIF_list,
    .CNV_ind_cancer,
    df = .TCGA_CNV_sampletype_subset
  )

  lapply(.EIF_CNV_ind_cancer, .CNV_barplot)
}


#' Correlation matrix for CNV values
#' @description generates correlation matrix for eIF4F CNV in tumors from all TCGA cancer type combined
#' @param EIF_list gene names in a vector of characters
#' @details  This function first selects unthreshold CNV values from only EIF4F gene
#' in the data frame \code{TCGA_CNV_value} prepared from \code{\link{initialize_cnv_data}}.
#'
#' Then it uses the subset data to calculate the correlation coefficients
#' between every two genes and plot the correlation matrix with the function \code{\link{.matrix_plot}}
#'
#' It should not be used directly, only inside \code{\link{EIF4F_CNV_analysis}} function.
#' @return the correlation matrix plot for \code{EIF_list} CNV values in TCGA tumors
#' @examples \dontrun{
#' plot_matrix_CNVcorr_TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))
#' }
#' @keywords internal
.plot_matrix_CNVcorr_TCGA <- function(EIF_list) {
  TCGA_CNV_value %>%
    dplyr::select(all_of(EIF_list)) %>%
    .matrix_plot()
}


#' CNV ratios in tumors vs adjacent normal tissues across tumor types
#' @description This function generates boxplot for CNV ratios in tumors vs adjacent normals from individual TCGA cancer types.
#' @param EIF_list gene names in a vector of characters
#' @details  This function first selects CNV and sample type data from only EIF4F gene
#' in the data frame \code{TCGA_CNVratio_sampletype} prepared from \code{\link{initialize_cnv_data}}.
#'
#' Then it uses the subset data \code{.TCGA_CNVratio_sampletype_subset} to perform CNV status analysis for each gene across
#' all tumors types with the function \code{\link{.CNVratio_tumor}} and plots the CNV results and \code{\link{.CNVratio_boxplot}}
#'
#' It should not be used directly, only inside \code{\link{EIF4F_CNV_analysis}} function.
#' @return stacked bar plots for grouped CNV status of \code{EIF_list} in TCGA tumors
#' @examples \dontrun{
#' plot_boxgraph_CNVratio_TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))
#' }
#' @keywords internal
.plot_boxgraph_CNVratio_TCGA <- function(EIF_list) {
  .TCGA_CNVratio_sampletype_subset <- NULL
  .TCGA_CNVratio_sampletype_subset <- TCGA_CNVratio_sampletype %>%
    dplyr::select(
      all_of(EIF_list),
      "sample.type",
      "primary.disease"
    )

  .EIF_CNVratio_ind_cancer <- lapply(EIF_list,
    .CNVratio_tumor,
    df = .TCGA_CNVratio_sampletype_subset
  )

  lapply(.EIF_CNVratio_ind_cancer, .CNVratio_boxplot)
}


## wrapper function to call all master functions with inputs ===================

#' Perform all CNV related analysis and generate plots
#' @description A wrapper function to call all master functions for CNV analysis with inputs
#'
#' @details  This function run three master functions together:
#' \code{\link{.plot_bargraph_CNV_TCGA}},
#' \code{\link{.plot_matrix_CNVcorr_TCGA}} and
#' \code{\link{.plot_boxgraph_CNVratio_TCGA}} with inputs
#'
#' @return CNV analysis plots
#'
#' @export
#'
#' @examples \dontrun{
#' EIF4F_CNV_analysis()
#' }
EIF4F_CNV_analysis <- function() {
  .plot_bargraph_CNV_TCGA(c(
    "TP53", "EIF4A1", "EIF4A2",
    "EIF4E", "EIF4E2", "EIF4E3",
    "MYC", "EIF3D", "EIF4EBP1",
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4H", "EIF4B"
  ))

  .plot_matrix_CNVcorr_TCGA(c(
    "TP53", "EIF4A1", "EIF4A2",
    "EIF4E", "EIF4E2", "EIF4E3",
    "MYC", "EIF3D", "EIF4EBP1",
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4H", "EIF4B"
  ))

  .plot_boxgraph_CNVratio_TCGA(c(
    "EIF4G1", "EIF4E", "EIF4A1", "EIF4EBP1", "EIF4G2", "EIF4G3",
    "EIF4E2", "EIF4E3", "EIF4A2", "EIF3D", "EIF4H", "EIF4B"
  ))
}
