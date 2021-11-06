# prepare CNV related datasets from TCGA ---------------------------------------

#' Read all CNV related datasets from TCGA
#'
#' @description This function reads all CNV related datasets from TCGA and generates three global variables.
#'
#' TCGA.CNV.value: the CNV value data generated from \code{\link{.get.TCGA.CNV.value}}
#'
#' TCGA.CNV.sampletype: the merged dataset from TCGA.CNV, the threshold CNV data generated from \code{\link{.get.TCGA.CNV}},
#' and TCGA.sampletype, the annotation data from the "TCGA_phenotype_denseDataOnlyDownload.tsv"
#' dataset with selection of sample.type and primary.disease columns. “Solid Tissue Normal” samples are excluded.
#'
#' TCGA.CNVratio.sampletype: the merged dataset from TCGA.CNV.ratio, the CNV ratio data generated from \code{\link{.get.TCGA.CNV.ratio}},
#' and TCGA.sampletype.
#'
#' @importFrom dplyr distinct filter select select_if mutate mutate_at summarise rename group_by all_of
#' @importFrom tibble remove_rownames column_to_rownames
#' @importFrom data.table fread transpose
#'
#' @export
#'
#' @examples \dontrun{initialize.cnv.data()}
#'
initialize.cnv.data <- function() {
  TCGA.CNV.value <<- .get.TCGA.CNV.value()

  TCGA.CNV <- .get.TCGA.CNV()

  TCGA.CNV.ratio <- .get.TCGA.CNV.ratio()

  TCGA.sampletype <- readr::read_tsv(file.path(
    data.file.directory,
    "TCGA_phenotype_denseDataOnlyDownload.tsv"
  )) %>%
    {
      #as_tibble(.) %>%
        as.data.frame(.) %>%
        dplyr::distinct(., sample, .keep_all = TRUE) %>%
        na.omit(.) %>%
        tibble::remove_rownames(.) %>%
        tibble::column_to_rownames(., var = "sample") %>%
        dplyr::select("sample_type", "_primary_disease") %>%
        rename("sample.type" = "sample_type",
               "primary.disease" = "_primary_disease")
    }

  TCGA.CNV.sampletype <<- merge(TCGA.CNV,
    TCGA.sampletype,
    by    = "row.names",
    all.x = TRUE
  ) %>%
    dplyr::filter(sample.type != "Solid Tissue Normal") %>%
    {
      tibble::remove_rownames(.) %>%
        tibble::column_to_rownames(., var = "Row.names")
    }

  TCGA.CNVratio.sampletype <<- merge(TCGA.CNV.ratio,
    TCGA.sampletype,
    by    = "row.names",
    all.x = TRUE
  ) %>%
    {
      tibble::remove_rownames(.) %>%
        tibble::column_to_rownames(var = "Row.names")
    }
}

#' Read CNV threshold dataset from TCGA
#' @description This function reads the threshold CNV data "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes",
#' which groups the CNV values in tumors into five categories.
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize.cnv.data}} function.
#' @return a data frame that contains CNV threshold data for all TCGA tumor
#' @examples \dontrun{.get.TCGA.CNV()}
#' @keywords internal
.get.TCGA.CNV <- function() {
  TCGA.pancancer <- data.table::fread(
    file.path(
      data.file.directory,
      "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes"
    ),
    data.table = FALSE
  ) %>%
    as_tibble(.) %>%
    # as.data.frame(.) %>%
    dplyr::distinct(., Sample, .keep_all = TRUE) %>%
    na.omit(.) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "Sample")

  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and column names in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return(TCGA.pancancer_transpose)
}

#' Read CNV value dataset from TCGA
#' @description This function reads the unthreshold CNV value data from TCGA "Gistic2_CopyNumber_Gistic2_all_data_by_genes".
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize.cnv.data}} function.
#' @return a data frame that contains CNV value data for each tumor
#' @examples \dontrun{.get.TCGA.CNV.value()}
#' @keywords internal
.get.TCGA.CNV.value <- function() {
  TCGA.pancancer <- fread(
    file.path(
      data.file.directory,
      "Gistic2_CopyNumber_Gistic2_all_data_by_genes"
    ),
    data.table = FALSE
  ) %>%
    as_tibble(.) %>%
    # as.data.frame(.) %>%
    dplyr::distinct(., Sample, .keep_all = TRUE) %>%
    na.omit(.) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "Sample")

  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and colnames in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return(TCGA.pancancer_transpose)
}

#' Read CNV ratio dataset from TCGA
#' @description This function reads the CNV ratio data "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena",
#' which calculates the ratios of CNV value in each tumor to the average CNV value from normals in the same cancer type.
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize.cnv.data}} function.
#' @return a data frame that contains CNV ratio data for each tumor
#' @examples \dontrun{.get.TCGA.CNV.ratio()}
.get.TCGA.CNV.ratio <- function() {
  TCGA.pancancer <- fread(
    file.path(
      data.file.directory,
      "broad.mit.edu_PANCAN_Genome_Wide_SNP_6_whitelisted.gene.xena"
    ),
    data.table = FALSE
  ) %>%
    as_tibble(.) %>%
    # as.data.frame(.) %>%
    dplyr::distinct(., sample, .keep_all = TRUE) %>%
    na.omit(.) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "sample")

  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and colnames in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return(TCGA.pancancer_transpose)
}

# CNV data analysis and plotting -----------------------------------------------

#' Calculates the frequency of CNV status in all TCGA cancer types combined
#' @description This function calculates the frequency of each CNV status across tumors in all TCGA cancer types for every EIF4F gene.
#' @details It should not be used directly, only inside \code{\link{plot.bargraph.CNV.TCGA}} function.
#' @param df \code{TCGA.CNV.sampletype.subset} generated inside \code{\link{plot.bargraph.CNV.TCGA}}
#' @return a summary table ranking the EIF4F gene by the frequencies of duplication (CNV threshold labeled as “1” in the dataset)
#' @importFrom reshape2 dcast melt
#' @examples \dontrun{.CNV.all.cancer(TCGA.CNV.sampletype.subset)}
#' @keywords internal
.CNV.all.cancer <- function(df) {
  TCGA.CNV.anno.subset.long <- melt(df,
    id = c(
      "sample.type",
      "primary.disease"
    ),
    value.name = "CNV"
  ) %>%
    mutate_if(is.character, as.factor)

  CNV.sum <-
    table(TCGA.CNV.anno.subset.long[, c("CNV", "variable")]) %>%
    # tibble(.) %>%
    as.data.frame(.) %>%
    mutate(CNV = factor(CNV, levels = c("-2", "-1", "0", "1", "2")))

  # reorder stack bars by the frequency of duplication.
  Freq.sum <- dcast(CNV.sum, variable ~ CNV, mean)
  CNV.sum$variable <- factor(CNV.sum$variable,
    levels = Freq.sum[order(Freq.sum$`1`), ]$variable
  )
  return(CNV.sum)
}

#' Stacked bar plots of the CNV status summary
#' @description This function generates the stacked bar plots to the summary of CNV statuses
#' @details This plot function uses dataset generated from \code{\link{.CNV.all.cancer}} function
#'
#' It should not be used directly, only inside \code{\link{plot.bargraph.CNV.TCGA}}function.
#' @param data summary table from \code{.CNV.all.cancer(TCGA.CNV.sampletype.subset)}
#' @return stacked bar plots showing the summary table of the CNV statuses
#' @examples \dontrun{.CNV.sum.barplot(.CNV.all.cancer(TCGA.CNV.sampletype.subset))}
#' @keywords internal
.CNV.sum.barplot <- function(data) {
  p1 <- ggplot2::ggplot(
    data,
    aes(
      fill = CNV,
      y = Freq,
      x = variable
    )
  ) +
    geom_bar(
      stat = "identity",
      position = "fill"
    ) +
    geom_col() +
    geom_text(aes(label = paste0(Freq / 100, "%")),
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
    path = file.path(output.directory, "CNV"),
    filename = "EIFCNVsum.pdf",
    plot = p1,
    width = 9,
    height = 9,
    useDingbats = FALSE
  )
}

#' Calculates the frequency of CNV status in individual TCGA cancer types
#' @description a data analysis function that calculates the frequency of CNV status for each EIF4F gene in individual TCGA cancer types.
#' @details It should not be used directly, only inside \code{\link{plot.bargraph.CNV.TCGA}} function.
#' @param df \code{TCGA.CNV.sampletype.subset} generated inside \code{\link{plot.bargraph.CNV.TCGA}}
#' @param x one gene from the input argument of \code{\link{plot.bargraph.CNV.TCGA}}
#' @return a list with the summary table of CNV in individual TCGA cancer types and gene name
#' @importFrom forcats fct_rev
#' @examples \dontrun{lapply(EIF,.CNV.ind.cancer,df = TCGA.CNV.sampletype.subset)}
#' @keywords internal
.CNV.ind.cancer <- function(df, x) {
  TCGA.CNV.anno.subset.long <- df %>%
    dplyr::select(
      all_of(x),
      "sample.type",
      "primary.disease"
    ) %>%
    melt(.,
      id = c(
        "sample.type",
        "primary.disease"
      ),
      value.name = "CNV"
    ) %>%
    mutate_if(is.character, as.factor)

  CNV.sum <-
    table(TCGA.CNV.anno.subset.long[, c("CNV", "primary.disease")]) %>%
    # tibble(.) %>%
    as.data.frame(.) %>%
    mutate(CNV = factor(CNV, levels = c("-2", "-1", "0", "1", "2"))) %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))

  output <- list(CNV.sum, x)
  return(output)
}

#' Stacked bar plots of the CNV status
#' @description This function generates stacked bar plots for CNV status of each gene in an individual cancer type.
#' @details This plot function uses dataset generated from \code{\link{.CNV.ind.cancer}} function
#' It should not be used directly, only inside \code{\link{plot.bargraph.CNV.TCGA}} function.
#' @param df \code{EIF.CNV.ind.cancer} generated inside \code{\link{plot.bargraph.CNV.TCGA}}
#' @return stacked bar plots for CNV status of each gene in an individual cancer type.
#' @examples \dontrun{lapply(EIF.CNV.ind.cancer, .CNV.barplot)}
#' @keywords internal
.CNV.barplot <- function(df) {
  p1 <- ggplot(
    df[[1]],
    aes(
      fill = CNV,
      order = as.numeric(CNV),
      y = Freq,
      x = primary.disease
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
    path = file.path(output.directory, "CNV"),
    filename = paste0(df[[2]], "pancancerCNV.pdf"),
    plot = p1,
    width = 7.5,
    height = 9,
    useDingbats = FALSE
  )
}

#' Correlation coefficients for CNV values
#' @description This function calculates the correlation coefficients between every two genes and plot the correlation matrix with the function.
#' @details This plot function uses dataset \code{TCGA.CNV.value} generated from \code{\link{initialize.cnv.data}} function
#' It should not be used directly, only inside \code{\link{plot.matrix.CNVcorr.TCGA}} function.
#' @param df output of \code{TCGA.CNV.value \%>\% select(all_of(EIF))} generated inside \code{\link{plot.matrix.CNVcorr.TCGA}}
#' @return stacked bar plots for CNV status of each gene in an individual cancer type.
#' @importFrom corrplot cor.mtest corrplot
#' @importFrom grDevices dev.off
#' @examples \dontrun{TCGA.CNV.value %>% dplyr::select(all_of(EIF)) %>% .matrix.plot()}
#' @keywords internal
.matrix.plot <- function(df) {
  # correlation plot
  M <- cor(df)
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
    file.path(output.directory, "CNV", "EIFCNVcormatrix.pdf"),
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
#' @details It should not be used directly, only inside \code{\link{plot.boxgraph.CNVratio.TCGA}} function.
#' @param df \code{TCGA.CNVratio.sampletype.subset} generated inside \code{\link{plot.boxgraph.CNVratio.TCGA}}
#' @param x one gene from the input argument of \code{\link{plot.boxgraph.CNVratio.TCGA}}
#' @return a list with the data frame of CNV ratio of the input gene in individual TCGA cancer types and gene name
#' @examples \dontrun{lapply(EIF,.CNVratio.tumor,df = TCGA.CNVratio.sampletype.subset)}
#' @keywords internal
.CNVratio.tumor <- function(df, x) {
  CNVratio.data <- df %>%
    dplyr::select(all_of(x), "sample.type", "primary.disease") %>%
    melt(.,
      id = c("sample.type", "primary.disease"),
      value.name = "CNV"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  output <- list(CNVratio.data, x)
  return(output)
}

#' Box plots of the CNV ratios in tumors vs adjacent normals
#' @description This function generates boxplot for CNV ratios in tumors vs adjacent normals from individual TCGA cancer types.
#' @details This plot function uses dataset \code{TCGA.CNVratio.sampletype.subset} generated from \code{\link{plot.boxgraph.CNVratio.TCGA}} function.
#'
#' It should not be used directly, only inside \code{\link{plot.boxgraph.CNVratio.TCGA}} function.
#' @param df \code{EIF.CNVratio.ind.cancer} generated inside \code{\link{plot.boxgraph.CNVratio.TCGA}}
#' @return boxplot for CNV ratios in tumors vs adjacent normals from individual TCGA cancer types.
#' @examples \dontrun{lapply(EIF.CNVratio.ind.cancer, .CNVratio.boxplot)}
#' @keywords internal
.CNVratio.boxplot <- function(df) {
  p1 <- ggplot(
    data = df[[1]],
    aes(
      y = 2**CNV,
      x = primary.disease,
      # x = f.ordered1,
      color = primary.disease
    )
  ) +
    # ylim(0, 3) +
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
    path = file.path(output.directory, "CNV"),
    filename = paste0(df[[2]], "pancancerCNVratio.pdf"),
    plot = p1,
    width = 7,
    height = 9,
    useDingbats = FALSE
  )
}


# master function to call CNV data analysis and plotting -----------------------
#' Summary of CNV statuses in bar plots
#' @description Provides the summary of CNV statuses of EIF4F genes in tumors from all TCGA
#' cancer types combined and in tumors from individual TCGA cancer types
#' @param EIF gene names in a vector of characters
#' @details  This function first selects CNV and sample type data from only EIF4F gene
#' in the data frame \code{TCGA.CNV.sampletype} prepared from \code{\link{initialize.cnv.data}}.
#'
#' Then it uses the subset data \code{TCGA.CNV.sampletype.subset} to perform the CNV status analysis of all inquired genes
#' with the functions \code{\link{.CNV.all.cancer}} and plot the results as a bar plot with \code{\link{.CNV.sum.barplot}}
#'
#' It also uses the same subset data to perform CNV status analysis for each gene across
#' all tumors types with the function \code{\link{.CNV.ind.cancer}} and plots the CNV results and \code{\link{.CNV.barplot}}
#' @return stacked bar plots for grouped CNV status of \code{EIF} in TCGA tumors
#' @export
#' @examples \dontrun{plot.bargraph.CNV.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))}
#' @examples \dontrun{plot.bargraph.CNV.TCGA(c("TP53", "EIF4A1", "EIF4A2","EIF4E",
#' "EIF4E2", "EIF4E3", "MYC", "EIF3D", "EIF4EBP1", "EIF4G1", "EIF4G2", "EIF4G3",
#' "PABPC1", "MKNK1", "MKNK2"))}
plot.bargraph.CNV.TCGA <- function(EIF) {
  # combine CNV data with annotation data
  TCGA.CNV.sampletype.subset <- TCGA.CNV.sampletype %>%
    dplyr::select(all_of(EIF), "sample.type", "primary.disease")

  # stacked bar plots for CNV status in combined TCGA tumors
  .CNV.all.cancer(TCGA.CNV.sampletype.subset) %>%
    .CNV.sum.barplot(.)

  # stacked bar plots for CNV status in each TCGA tumor type
  EIF.CNV.ind.cancer <- lapply(EIF,
    .CNV.ind.cancer,
    df = TCGA.CNV.sampletype.subset
  )

  lapply(EIF.CNV.ind.cancer, .CNV.barplot)
}


#' Correlation matrix for CNV values
#' @description generates correlation matrix for eIF4F CNV in tumors from all TCGA cancer type combined
#' @param EIF gene names in a vector of characters
#' @details  This function first selects unthreshold CNV values from only EIF4F gene
#' in the data frame \code{TCGA.CNV.value} prepared from \code{\link{initialize.cnv.data}}.
#'
#' Then it uses the subset data to calculate the correlation coefficients
#' between every two genes and plot the correlation matrix with the function \code{\link{.matrix.plot}}
#' @return the correlation matrix plot for \code{EIF} CNV values in TCGA tumors
#' @export
#' @examples \dontrun{plot.matrix.CNVcorr.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))}
plot.matrix.CNVcorr.TCGA <- function(EIF) {
  TCGA.CNV.value %>%
    dplyr::select(all_of(EIF)) %>%
    .matrix.plot()
}


#' CNV ratios in tumors vs adjacent normal tissues across tumor types
#' @description This function generates boxplot for CNV ratios in tumors vs adjacent normals from individual TCGA cancer types.
#' @param EIF gene names in a vector of characters
#' @details  This function first selects CNV and sample type data from only EIF4F gene
#' in the data frame \code{TCGA.CNVratio.sampletype} prepared from \code{\link{initialize.cnv.data}}.
#'
#' Then it uses the subset data \code{TCGA.CNVratio.sampletype.subset} to perform CNV status analysis for each gene across
#' all tumors types with the function \code{\link{.CNVratio.tumor}} and plots the CNV results and \code{\link{.CNVratio.boxplot}}
#' @return stacked bar plots for grouped CNV status of \code{EIF} in TCGA tumors
#' @export
#' @examples \dontrun{plot.boxgraph.CNVratio.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))}
plot.boxgraph.CNVratio.TCGA <- function(EIF) {
  TCGA.CNVratio.sampletype.subset <- TCGA.CNVratio.sampletype %>%
    dplyr::select(
      all_of(EIF),
      "sample.type",
      "primary.disease"
    )

  EIF.CNVratio.ind.cancer <- lapply(EIF,
    .CNVratio.tumor,
    df = TCGA.CNVratio.sampletype.subset
  )

  lapply(EIF.CNVratio.ind.cancer, .CNVratio.boxplot)
}
