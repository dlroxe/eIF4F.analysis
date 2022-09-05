# Analyses on EIF4F correlating genes (CORs)
# This R script contains four sections.
#
# (1) identify positively or negatively correlating genes of EIF4F
#
# (2) analyze and plot the correlating genes
#
# (3) composite functions to execute a pipeline of functions to select related
#  RNASeq data for correlation analyses and plotting.
#
# (4) wrapper function that serves as the entry point for this file
#


## Identify correlating genes for EIF4F genes ==================================

#' @title Identify EIF4F correlating genes
#'
#' @description This function
#'
#' * takes the specific data frames
#'  `.TCGA_GTEX_RNAseq_sampletype_subset` and `sample_type` that are
#'  generated inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' * calculates the correlation coefficiency between each EIF4F gene
#'  and the rest of cellular mRNAs with [.correlation_coefficient()]
#'
#' * combines the correlation coefficiency data from EIF4E, EIF4A1, EIF4G1, and
#'  EIF4EBP1
#'
#' * selects positive correlating genes with [.is_significant_poscor()] and
#'  negative correlating genes with [.is_significant_negcor()]
#'
#' * summarizes the total number of posCORs or negCORs identified for
#'  each EIF4F gene with [.summarize_counts()]
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' @family helper function to identify correlating genes for EIF4F genes
#'
#' @param df the data frame `.TCGA_GTEX_RNAseq_sampletype_subset` generated
#'  inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @param sample_type sample types, either `all.tumor.type` or
#'  `c("Normal Tissue")` generated inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @return
#'
#' a list output with four elements:
#'  * `cor_value_combined` for the heatmap
#'  * `CORs_summary_tbl` for bargraph
#'  * `posCOR_EIF4F` for Venn plots
#'  * `negCOR_EIF4F` for Venn plots
#'
#' @importFrom AnnotationDbi mapIds
#'
#' @importFrom clusterProfiler compareCluster dotplot
#'
#' @importFrom purrr discard
#'
#' @importFrom stats setNames
#'
#' @examples \dontrun{
#' .EIF_correlation(
#'   df = .TCGA_GTEX_RNAseq_sampletype_subset,
#'   sample_type = all.tumor.type
#' )
#' }
#'
#' @keywords internal
#'
.EIF_correlation <- function(df, sample_type) {
  TCGA_GTEX_RNAseq_subset <- df[df$sample.type %in% sample_type, ] %>%
    stats::na.omit() # select tumor or healthy samples
  # provide a list of all cellular genes for the correlation analysis
  Gene_ID <- names(df) %>%
    purrr::discard(~ .x %in% c(
      "Row.names",
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ))

  # find all genes positively correlate with EIF4F expression
  # lapply function gives a large list, need to convert it to a dataframe
  EIF_cor_value <- function(EIF) {
    cor.data <- do.call(
      rbind.data.frame,
      lapply(Gene_ID,
        .correlation_coefficient,
        gene2 = EIF,
        df = TCGA_GTEX_RNAseq_subset
      )
    )
    rownames(cor.data) <- cor.data[, 1]
    return(cor.data)
  }

  EIF4E.cor <- EIF_cor_value("EIF4E")
  EIF4G1.cor <- EIF_cor_value("EIF4G1")
  EIF4A1.cor <- EIF_cor_value("EIF4A1")
  EIF4EBP1.cor <- EIF_cor_value("EIF4EBP1")

  cor_value_combined <- cbind(
    stats::setNames(data.frame(EIF4E.cor[, c(3, 4)]), c("EIF4E", "EIF4E.p")),
    stats::setNames(data.frame(EIF4G1.cor[, c(3, 4)]), c("EIF4G1", "EIF4G1.p")),
    stats::setNames(data.frame(EIF4A1.cor[, c(3, 4)]), c("EIF4A1", "EIF4A1.p")),
    stats::setNames(data.frame(EIF4EBP1.cor[, c(3, 4)]), c(
      "EIF4EBP1",
      "EIF4EBP1.p"
    ))
  )

  # identify posCORs for each EIF4F gene
  posCOR_EIF4F <- cbind(
    .is_significant_poscor(EIF4E.cor$estimate, EIF4E.cor$p.value),
    .is_significant_poscor(EIF4G1.cor$estimate, EIF4G1.cor$p.value),
    .is_significant_poscor(EIF4A1.cor$estimate, EIF4A1.cor$p.value),
    .is_significant_poscor(EIF4EBP1.cor$estimate, EIF4EBP1.cor$p.value)
  ) %>%
    `colnames<-`(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1"))

  # identify negCORs for each EIF4F gene
  negCOR_EIF4F <- cbind(
    .is_significant_negcor(EIF4E.cor$estimate, EIF4E.cor$p.value),
    .is_significant_negcor(EIF4G1.cor$estimate, EIF4G1.cor$p.value),
    .is_significant_negcor(EIF4A1.cor$estimate, EIF4A1.cor$p.value),
    .is_significant_negcor(EIF4EBP1.cor$estimate, EIF4EBP1.cor$p.value)
  ) %>%
    `colnames<-`(c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1"))

  CORs_summary_tbl <- cbind(
    .summarize_counts(posCOR_EIF4F, "posCORs"),
    .summarize_counts(negCOR_EIF4F, "negCORs")
  )

  ## four output values:
  ## cor_value_combined for the heatmap function
  ## CORs_summary_tbl for the bargraph function
  ## posCOR_EIF4F and negCOR_EIF4F for Venn plots
  return(list(cor_value_combined, CORs_summary_tbl, posCOR_EIF4F, negCOR_EIF4F))
}


#' @title Calculate the correlation co-efficiency between two gene expressions
#'
#' @description
#'
#' A helper function calculates the correlation co-efficiency between two gene
#'  expressions, using the RNAseq expression dataset.
#'
#' It should not be used directly, only inside [.EIF_correlation()]
#'  function.
#'
#' @family helper function to identify correlating genes for EIF4F genes
#'
#' @param gene1 gene name
#'
#' @param gene2 gene name
#'
#' @param df expression dataset
#'
#' @return a dataframe with correlation coefficiency between two gene expression
#'
#' @importFrom stats cor.test
#'
#'
#' @keywords internal
#'
.correlation_coefficient <- function(gene1, gene2, df) {
  result <- stats::cor.test(df[[gene1]], df[[gene2]], method = "pearson")

  return(data.frame(gene1,
    gene2,
    result[c("estimate", "p.value", "statistic", "method")],
    stringsAsFactors = FALSE
  ))
}


#' @title Select positive correlating genes based on coefficiency and pvalue
#'
#' @description
#'
#' A helper function selects positive correlating genes based on coefficiency
#'  and pvalue.
#'
#' It should not be used directly, only inside [.EIF_correlation()]
#'  function.
#'
#' @family helper function to identify correlating genes for EIF4F genes
#'
#' @param magnitude correlation coefficiency value
#'
#' @param pvalue p value for Pearson correlation
#'
#' @return a data frame of positive correlating gene with correlation
#'  coefficiency and p value
#'
#' @keywords internal
#'
.is_significant_poscor <- function(magnitude, pvalue) {
  return(magnitude > 0.3 & pvalue <= 0.05)
}


#' @title Select negative correlating genes based on coefficiency and pvalue
#'
#' @description
#'
#' A helper function selects negative correlating genes based on coefficiency
#'  and pvalue.
#'
#' It should not be used directly, only inside [.EIF_correlation()]
#'  function.
#'
#' @family helper function to identify correlating genes for EIF4F genes
#'
#' @param magnitude correlation coefficiency value
#'
#' @param pvalue p value for Pearson correlation
#'
#' @return a dataframe of negative correlating gene with correlation
#'  coefficiency and p value
#'
#' @keywords internal
#'
.is_significant_negcor <- function(magnitude, pvalue) {
  return(magnitude < -0.3 & pvalue <= 0.05)
}


#' @title Summary of total number of CORs identified for each EIF4F gene
#'
#' @description
#'
#' This function summarizes the total number of posCORs or negCORs identified for
#'  each EIF4F gene
#'
#' It should not be used directly, only inside [.EIF_correlation()]
#'  function.
#'
#' @family helper function to identify correlating genes for EIF4F genes
#'
#' @param identified_COR a data frame with logic statement for each gene tested
#'
#' @param COR_type posCORs or negCORs
#'
#' @return a table with counts of posCORs or negCORs for each EIF4F gene
#'
#' @importFrom dplyr bind_rows
#'
#'
#' @keywords internal
#'
.summarize_counts <- function(identified_COR, COR_type) {
  return(
    as.data.frame(identified_COR) %>%
      lapply(table) %>%
      dplyr::bind_rows(.id = "column_label") %>%
      as.data.frame() %>%
      tibble::column_to_rownames(var = "column_label") %>%
      dplyr::select("TRUE") %>%
      dplyr::rename(!!COR_type := "TRUE")
  )
}


## Correlation analysis and plotting ===========================================

#' @title Plot Venn diagrams for correlating genes
#'
#' @description
#'
#' This function draw Venn diagrams, using the correlation data generated from
#'  [.EIF_correlation()].
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' Side effects:
#'
#' (1) Venn diagrams on screen and as pdf file to show the overlap of
#'  EIF correlating genes
#'
#' @family helper function for correlation analysis plotting
#'
#' @param df correlation data
#'
#' @param tissue_type input argument of [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @param sample_type tumor or normal for the title of Venn diagram
#'
#' @param CORs_type posCOR or negCORs for the title of Venn diagram
#'
#' @importFrom eulerr euler
#'
#' @importFrom limma vennCounts vennDiagram
#'
#' @examples \dontrun{
#' .CORs_vennplot(
#' df = EIF.cor.tumor[[3]],
#' tissue_type = tissue_type,
#' sample_type = "tumor",
#' CORs_type = "posCOR")
#' }
#'
#' @keywords internal
#'
.CORs_vennplot <- function(df, tissue_type, sample_type, CORs_type) {
  b <- limma::vennCounts(df)
  colnames(b) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1", "Counts")
  limma::vennDiagram(b)
  ## eulerr generates area-proportional Euler diagrams that display set
  ## relationships (intersections, unions, and disjoints) with circles.
  pos.Venn2 <- eulerr::euler(
    c(
      EIF4E       = b[9, "Counts"], # EIF4E
      EIF4G1      = b[5, "Counts"], # EIF4G1
      EIF4A1      = b[3, "Counts"], # EIF4A1
      EIF4EBP1    = b[2, "Counts"], # EIF4EBP1
      "EIF4E&EIF4G1"   = b[13, "Counts"],
      "EIF4E&EIF4A1"   = b[11, "Counts"],
      "EIF4E&EIF4EBP1" = b[10, "Counts"],
      "EIF4G1&EIF4A1"   = b[7, "Counts"],
      "EIF4G1&EIF4EBP1"   = b[6, "Counts"],
      "EIF4A1&EIF4EBP1"   = b[4, "Counts"],
      "EIF4E&EIF4G1&EIF4A1" = b[15, "Counts"],
      "EIF4E&EIF4G1&EIF4EBP1" = b[14, "Counts"],
      "EIF4E&EIF4A1&EIF4EBP1" = b[12, "Counts"],
      "EIF4G1&EIF4A1&EIF4EBP1" = b[8, "Counts"],
      "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]
    ),
    # shape = "ellipse"
  )
  p2 <- plot(pos.Venn2,
    # key = TRUE,
    main = paste(tissue_type, sample_type, CORs_type),
    lwd = 0,
    fill = c(
      "#999999", "#009E73",
      "#56B4E9", "#E69F00"
    ),
    quantities = list(cex = 1.25),
    labels = list(
      labels = c(
        "EIF4E", "EIF4G1",
        "EIF4A1", "EIF4EBP1"
      ),
      cex = 1.25
    )
  )
  print(p2)
  ggplot2::ggsave(
    path = file.path(output_directory, "CORs"),
    filename = paste(tissue_type, sample_type, CORs_type, "Venn.pdf"),
    plot = p2,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )

  return(NULL)
}


#' @title Combine the summary of number of posCORs or negCORs from sample types
#'
#' @description
#'
#' This function combines the summary of number of posCORs or negCORs from
#'  sample types, using the CORs table generated from [.EIF_correlation()].
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' @family helper function for correlation analysis
#'
#' @param COR_tbl1 a COR summary table of specific sample type `EIF.cor.tumor[[2]]`
#'
#' @param sample_type1 the label of sample type where `COR_df1` is derived from
#'
#' @param COR_tbl2 a COR summary table of specific sample type `EIF.cor.normal[[2]]`
#'
#' @param sample_type2 the label of sample type where `COR_df2` is derived from
#'
#' @return a data frame with combined COR summary table
#'
#' @examples \dontrun{
#' .combine_CORs_summary(
#' COR_tbl1 = EIF.cor.tumor[[2]], sample_type1 = "tumor",
#' COR_tbl2 = EIF.cor.normal[[2]], sample_type2 = "normal")
#' }
#'
#' @keywords internal
#'
.combine_CORs_summary <- function(COR_tbl1, sample_type1,
                                  COR_tbl2, sample_type2) {
  EIF.cor.counts.tumor <- COR_tbl1 %>%
    tibble::add_column(label = sample_type1) %>%
    tibble::rownames_to_column(var = "gene")

  EIF.cor.counts.normal <- COR_tbl2 %>%
    tibble::add_column(label = sample_type2) %>%
    tibble::rownames_to_column(var = "gene")

  return(rbind(EIF.cor.counts.tumor,
    EIF.cor.counts.normal,
    make.row.names = F
  ) %>%
    dplyr::mutate(label = factor(.data$label,
      levels = c("tumor", "normal"),
      labels = c("Tumors", "Healthy tissues")
    )) %>%
    dplyr::mutate(gene = factor(.data$gene,
      levels = c(
        "EIF4EBP1", "EIF4A1",
        "EIF4G1", "EIF4E"
      )
    )))
}

#' @title Plot bargraph for the numbers of correlating genes
#'
#' @description
#'
#' This function draw bargraph, using the correlation data generated from
#'  [.EIF_correlation()].
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' Side effects:
#'
#' (1) bar graphs on screen and as pdf file to show the numbers of
#'  identified correlating genes for each EIF4F subunit
#'
#' @family helper function for correlation analysis plotting
#'
#' @param df combined CORs summary table
#'
#' @param tissue_type input argument of [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @param CORs_type posCORs or negCORs for the title of Venn diagram
#'
#' @param coord_flip.ylim the limit of y axis in the bar plot
#'
#' @return bargraph for the numbers of posCOR or negCORs
#'
#' @examples \dontrun{
#' .CORs_summary_bargraph(
#'  df = EIF.cor,
#'  tissue_type = tissue_type,
#'  CORs_type = "posCORs",
#'  coord_flip.ylim = 14000)
#' }
#'
#' @keywords internal
#'
.CORs_summary_bargraph <- function(df, tissue_type,
                                   CORs_type, coord_flip.ylim) {
  p1 <- ggplot(
    data = df,
    aes_string(
      x = "gene",
      y = CORs_type,
      # y = !!sym(CORs_type), # quote the passed variable CORs_type
      #color = "label",
      fill = "label"),
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge()
    ) +
    geom_text(aes_string(label = CORs_type),
      position = position_dodge(width = 0.9),
      size = 3.5
    ) +
    scale_fill_manual(values = c(
      "#CC79A7", "#0072B2", "#E69F00",
      "#009E73", "#D55E00"
    )) + # for color-blind palettes
    labs(y = paste("number of ", tissue_type, CORs_type)) +
    coord_flip(ylim = c(0, coord_flip.ylim)) +
    # Flip ordering of legend without altering ordering in plot
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_bw() +
    theme(
      plot.title = black_bold_18,
      axis.title.x = black_bold_18,
      axis.title.y = element_blank(),
      axis.text.x = black_bold_18,
      axis.text.y = black_bold_18,
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_18,
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_18
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "CORs"),
    filename = paste(tissue_type, CORs_type, "bargraph", ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )

  return(NULL)
}


#' @title Combine all the EIF4F correlating gene data from tumor and healthy
#'  tissues
#'
#' @description This function
#'
#' * combines all the EIF4F correlating genes and their correlation
#'  coefficients to a data frame, with the data generated from
#'  [.EIF_correlation()].
#' * selects significant correlating genes with [.is_significant_correlation()]
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' @family helper function for correlation analysis
#'
#' @param df1 `cor_value_combined` of `EIF.cor.tumor[[2]]`
#'
#' @param df2 `cor_value_combined` of `EIF.cor.normal[[2]]`
#'
#' @return
#'
#' data frame with posCORs or negCORs from both tumors and healthy tissue samples
#'
#' @importFrom stats setNames
#'
#' @examples \dontrun{
#' .combine_COR_list(df1 = EIF.cor.tumor[[1]], df2 = EIF.cor.normal[[1]])
#' }
#'
#' @keywords internal
#'
.combine_CORs_list <- function(df1, df2) {
  cor.data <- cbind(
    stats::setNames(
      data.frame(df1[1:8]),
      c(
        "EIF4E.tumor", "EIF4E.p.tumor",
        "EIF4G1.tumor", "EIF4G1.p.tumor",
        "EIF4A1.tumor", "EIF4A1.p.tumor",
        "EIF4EBP1.tumor", "EIF4EBP1.p.tumor"
      )
    ),
    stats::setNames(
      data.frame(df2[1:8]),
      c(
        "EIF4E.normal", "EIF4E.p.normal",
        "EIF4G1.normal", "EIF4G1.p.normal",
        "EIF4A1.normal", "EIF4A1.p.normal",
        "EIF4EBP1.normal", "EIF4EBP1.p.normal"
      )
    )
  )
  DF <- as.matrix(na.omit(cor.data[
    .is_significant_correlation(cor.data$EIF4E.tumor, cor.data$EIF4E.p.tumor)
    | .is_significant_correlation(cor.data$EIF4G1.tumor, cor.data$EIF4G1.p.tumor)
    | .is_significant_correlation(cor.data$EIF4A1.tumor, cor.data$EIF4A1.p.tumor)
    | .is_significant_correlation(cor.data$EIF4EBP1.tumor, cor.data$EIF4EBP1.p.tumor)
    | .is_significant_correlation(cor.data$EIF4E.normal, cor.data$EIF4E.p.normal)
    | .is_significant_correlation(cor.data$EIF4G1.normal, cor.data$EIF4G1.p.normal)
    | .is_significant_correlation(cor.data$EIF4A1.normal, cor.data$EIF4A1.p.normal)
    | .is_significant_correlation(cor.data$EIF4EBP1.normal, cor.data$EIF4EBP1.p.normal),
  ]))

  return(DF[, c(1, 3, 5, 7, 9, 11, 13, 15)])
}


#' @title Select both positive and negative correlating genes based on
#'  coefficiency and pvalue
#'
#' @description
#'
#' This function selects both positive and negative correlating genes based on
#'  coefficiency and pvalue.
#'
#' It should not be used directly, only inside [.combine_CORs_list()]
#'  function.
#'
#' @family helper function for correlation analysis
#'
#' @param magnitude correlation coefficiency value
#'
#' @param pvalue p value for Pearson correlation
#'
#' @return a data frame of both positive and negative correlating gene with
#'  correlation coefficiency and p value
#'
#' @keywords internal
#'
.is_significant_correlation <- function(magnitude, pvalue) {
  return(.is_significant_poscor(magnitude, pvalue)
         | .is_significant_negcor(magnitude, pvalue))
}


#' @title Visualize associations between different sources of correlation data
#'  and reveal potential pattern
#'
#' @description
#'
#' This function draw heatmap of combined COR list data generated from
#'  [.EIF_correlation()].
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' Side effects:
#'
#' (1) heatmap on screen and as pdf file to show correlation strength
#'  and clustering pattern of EIF4F correlating genes.
#'
#' @family helper function for correlation analysis plotting
#'
#' @param df combined COR data generated from [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @param tissue_type "All" for all tissue type, or "Lung" for specific tissue
#'
#' @return
#'
#' a data structure that includes clustering analysis results.
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_text draw
#'
#' @importFrom grid gpar unit
#'
#' @importFrom circlize colorRamp2
#'
#' @examples \dontrun{
#' .CORs_coeff_heatmap(df = DF, tissue_type = tissue_type)
#' }
#'
#' @keywords internal
#'
.CORs_coeff_heatmap <- function(df, tissue_type) {
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation
  ## for more options)
  ht1 <- ComplexHeatmap::Heatmap(df,
    name = paste("Correlation Coefficient Heatmap", tissue_type),
    heatmap_legend_param = list(
      labels_gp = grid::gpar(font = 15),
      legend_width = grid::unit(6, "cm"),
      direction = "horizontal"
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = ComplexHeatmap::HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      type = c(
        "tumor", "tumor", "tumor", "tumor",
        "normal", "normal", "normal", "normal"
      ),
      col = list(type = c(
        "normal" = "royalblue",
        "tumor" = "pink"
      )),
      cn = ComplexHeatmap::anno_text(gsub("\\..*", "", colnames(df)),
        location = 0,
        rot = 0,
        just = "center",
        gp = grid::gpar(
          fontsize = 15,
          fontface = "bold"
        )
      )
    ),
    # cluster_rows  = as.dendrogram(hclust(dist(DF))),
    # row_split     = 3,
    row_km = 3,
    # row_km_repeats = 100,
    row_title = "cluster_%s",
    row_title_gp = grid::gpar(
      fontsize = 15,
      fontface = "bold"
    ),
    border = TRUE,
    col = circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "#EEEEEE", "red")
    )
  )
  ht <- ComplexHeatmap::draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )

  pdf(file.path(
    path = file.path(output_directory, "CORs"),
    filename = paste(tissue_type, "tumors heatmap.pdf")
  ),
  width = 8,
  height = 10,
  useDingbats = FALSE
  )
  ht <- ComplexHeatmap::draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )
  dev.off()

  return(ht1)
}


#' @title Retrieve the gene names from each cluster in heatmap
#'
#' @description
#'
#' This function retrieve the gene names from the clusters in heatmap from
#'  [.CORs_coeff_heatmap()].
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' @family helper function for correlation analysis
#'
#' @param df1 combined COR data generated from [.combine_CORs_list()]
#'
#' @param df2 gene name in the order of heatmap generated from
#'  [.CORs_coeff_heatmap()]
#'
#' @return
#'
#' a heatmap plot and dataframe with gene names of clustering
#'
#' @importFrom ComplexHeatmap row_order
#'
#' @importFrom AnnotationDbi mapIds
#'
#' @examples \dontrun{
#' .get_cluster_genes(df1 = DF, df2 = ht1)
#' }
#'
#' @keywords internal
#'
.get_cluster_genes <- function(df1, df2) {
  cluster.geneID.list <- function(cluster_label) {
    c1 <- row.names(df1[ComplexHeatmap::row_order(df2)[[cluster_label]], ]) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      stats::setNames("gene")
    # c1$V1 <- as.character(c1$V1)
    c1$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
      keys = c1$gene,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    # c1 <- c1[!is.na(c1)]
    c1 <- stats::na.omit(c1)
    return(c1$entrez)
  }
  cluster.num <- as.character(c(1:3))
  names(cluster.num) <- paste("cluster", 1:3)

  return(lapply(cluster.num, cluster.geneID.list))
}


#' @title Plot the enriched pathways from gene lists
#'
#' @description
#'
#' This function plots pathway enrichment analysis results of the gene from each
#'  cluster in heatmap generated from [.get_cluster_genes()].
#'
#' It should not be used directly, only inside [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  function.
#'
#' Side effects:
#'
#' (1) dotplot on screen and as pdf file to show the enriched pathways
#'  in the genes from different clusters in the heatmap
#'
#' @family helper function for correlation analysis plotting
#'
#' @param df pathway enrichment analysis results generated from
#'  [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @param pathway pathway name, "GO", ""REACTOME", "KEGG"
#'
#' @param tissue_type tissue type, the same argument of
#'  [.plot_Corr_RNAseq_TCGA_GTEX()]
#'
#' @importFrom clusterProfiler dotplot
#'
#' @examples \dontrun{
#' .pathway_dotplot(df = ck.REACTOME, pathway = "REACTOME",
#'  tissue_type = tissue_type)
#' }
#'
#' @keywords internal
#'
.pathway_dotplot <- function(df, pathway, tissue_type) {
  p1 <- clusterProfiler::dotplot(df,
    title = paste("The Most Enriched", pathway, "Pathways"),
    showCategory = 8,
    font.size = 10,
    includeAll = FALSE
  ) +
    theme_bw() +
    theme(
      plot.title = black_bold_16,
      axis.title = black_bold_12,
      axis.text.x = black_bold_12, #black_bold_16
      axis.text.y = black_bold_12, #black_bold_16
      #axis.line.x = element_line(color = "black"),
      #axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.title = black_bold_12,
      legend.text = black_bold_12,
      strip.text = black_bold_12
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "CORs"),
    filename = paste(tissue_type, " tumors enriched", pathway, ".pdf"),
    plot = p1,
    width = 12, #10
    height = 15, #10
    useDingbats = FALSE
  )

  return(NULL)
}


## Composite functions to analyze the EIF4F CORs ===============================

#' @title Perform correlation analysis on RNAseq data from all tumors and
#'  healthy tissues
#'
#' @description This function
#'
#' * finds correlating genes of EIF4F from the data frame `TCGA_GTEX_RNAseq_sampletype`
#'  - a global variable generated by [initialize_RNAseq_data()].
#' * identifies positively or negatively correlating genes from tumor or
#'  healthy samples by [.EIF_correlation()],
#' * finds the overlap of CORs for EIF4F subunits with [.CORs_vennplot()],
#' * compares the number of CORs for EIF4F subunits with
#'  [.combine_CORs_summary()] and plot the results in bargraph
#'  with [.CORs_summary_bargraph()]
#' * compares the correlation strength of CORs from healthy tissues and tumors
#'  from merged data frame generated from [.combine_CORs_list()] and visualizes
#'  clustering results in heatmap with [.CORs_coeff_heatmap()]
#' * retrieve the gene names from the clusters in heatmap with
#'  [.get_cluster_genes()]
#' * performs the pathways enrichment analysis on the retrieved gene list and
#'  plot the results with [.pathway_dotplot()]
#'
#' This function is not accessible to the user and will not show at the users'
#'  workspace. It can only be called by the exported [EIF4F_Corrgene_analysis()]
#'  function.
#'
#' Side effects:
#'
#' (1) Venn diagrams on screen and as pdf file to show the overlap of
#'  EIF correlating genes identify from RNAseq of `tissue_type`
#'
#' (2) bar graphs on screen and as pdf file to show the numbers of
#'  identified correlating genes for each EIF4F subunit identify from
#'  RNAseq of `tissue_type`
#'
#' (3) heatmap on screen and as pdf file to show correlation strength
#'  and clustering pattern of EIF4F correlating genes
#'
#' (4) dotplot on screen and as pdf file to show the enriched pathways
#'  in the genes from different clusters in the heatmap
#'
#' @family composite function to call correlation analysis and plotting
#'
#' @param tissue_type all tissue/tumor types or one specific type
#'
#' @importFrom purrr map pluck discard
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' .plot_Corr_RNAseq_TCGA_GTEX(x = "All")
#' .plot_Corr_RNAseq_TCGA_GTEX(x = "Lung")
#' }
#'
.plot_Corr_RNAseq_TCGA_GTEX <- function(tissue_type) {
  .TCGA_GTEX_RNAseq_sampletype_subset <- TCGA_GTEX_RNAseq_sampletype %>%
    dplyr::filter(if (tissue_type == "All") {
      TRUE
    } else {
      .data$primary.site == tissue_type
    }) %>%
    stats::na.omit() %>%
    # mutate_if(is.character, as.factor) %>%
    dplyr::mutate_at(c(
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ), factor)

  all.tumor.type <- .TCGA_GTEX_RNAseq_sampletype_subset %>%
    dplyr::select(.data$sample.type) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    purrr::map(levels) %>%
    purrr::pluck("sample.type") %>%
    purrr::discard(~ .x %in% c(
      "Cell Line",
      "Normal Tissue",
      "Solid Tissue Normal"
    ))

  # identify CORs for EIF4F genes from normal tissues
  EIF.cor.tumor <- .EIF_correlation(
    df = .TCGA_GTEX_RNAseq_sampletype_subset,
    sample_type = all.tumor.type
  )
  .CORs_vennplot(
    df = EIF.cor.tumor[[3]],
    tissue_type = tissue_type,
    sample_type = "tumor",
    CORs_type = "posCOR"
  )
  .CORs_vennplot(
    df = EIF.cor.tumor[[4]],
    tissue_type = tissue_type,
    sample_type = "tumor",
    CORs_type = "negCOR"
  )

  # identify CORs for EIF4F genes from normal tissues
  EIF.cor.normal <- .EIF_correlation(
    df = .TCGA_GTEX_RNAseq_sampletype_subset,
    sample_type = c("Normal Tissue")
  )
  .CORs_vennplot(
    df = EIF.cor.normal[[3]],
    tissue_type = tissue_type,
    sample_type = "normal",
    CORs_type = "posCOR"
  )
  .CORs_vennplot(
    df = EIF.cor.normal[[4]],
    tissue_type = tissue_type,
    sample_type = "normal",
    CORs_type = "negCOR"
  )

  # combine summary of CORs counts for EIF4F genes in tumors and normal tissues
  # for bargraph comparison
  EIF.cor <- .combine_CORs_summary(
    COR_tbl1 = EIF.cor.tumor[[2]], sample_type1 = "tumor",
    COR_tbl2 = EIF.cor.normal[[2]], sample_type2 = "normal"
  )
  .CORs_summary_bargraph(
    df = EIF.cor,
    tissue_type = tissue_type,
    CORs_type = "posCORs",
    coord_flip.ylim = 14000
  )
  .CORs_summary_bargraph(
    df = EIF.cor,
    tissue_type = tissue_type,
    CORs_type = "negCORs",
    coord_flip.ylim = 14000
  )

  # combine CORs data with coefficiency value for EIF4F genes in tumors and
  # normal tissues for clustering and heatmap analysis
  DF <- .combine_CORs_list(
    df1 = EIF.cor.tumor[[1]],
    df2 = EIF.cor.normal[[1]]
  )
  ht1 <- .CORs_coeff_heatmap(df = DF, tissue_type = tissue_type)

  # retrieve the gene list from clusters and perform pathway enrichment analysis
  cluster.data <- .get_cluster_genes(df1 = DF, df2 = ht1)

  ck.GO <- clusterProfiler::compareCluster(
    geneClusters = cluster.data,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db"
  )

  ck.KEGG <- clusterProfiler::compareCluster(
    geneClusters = cluster.data,
    fun = "enrichKEGG"
  )

  ck.REACTOME <- clusterProfiler::compareCluster(
    geneClusters = cluster.data,
    fun = "enrichPathway",
  )

  .pathway_dotplot(df = ck.GO, pathway = "GO", tissue_type = tissue_type)
  .pathway_dotplot(df = ck.KEGG, pathway = "KEGG", tissue_type = tissue_type)
  .pathway_dotplot(df = ck.REACTOME, pathway = "REACTOME",
                   tissue_type = tissue_type)

  return(NULL)
}


## Wrapper function to call all composite functions with inputs ================

#' @title Perform differential correlation analyses and generate plots
#'
#' @description
#'
#' A wrapper function to call all composite functions for PCA with inputs
#'
#' @details
#'
#' This function run the composite functions  [.plot_Corr_RNAseq_TCGA_GTEX()]
#'  with two inputs.
#'
#'  * "All" for all cancer types of TCGA and all healthy tissues types of GTEx.
#'  * "Lung" for lung cancer types (LUAD and LUSC of TCGA) and healthy lung
#'   tissues of GTEx
#'
#' Side effects:
#'
#' (1) Venn diagrams on screen and as pdf file to show the overlap of
#'  EIF correlating genes identify from RNAseq of `tissue_type`
#'
#' (2) bar graphs on screen and as pdf file to show the numbers of
#'  identified correlating genes for each EIF4F subunit identify from
#'  RNAseq of `tissue_type`
#'
#' (3) heatmap on screen and as pdf file to show correlation strength
#'  and clustering pattern of EIF4F correlating genes.
#'
#' (4) dotplot on screen and as pdf file to show the enriched pathways
#'  in the genes from different clusters in the heatmap
#'
#' @family wrapper function to call all composite functions with inputs
#'
#' @export
#'
#' @examples \dontrun{
#' EIF4F_Corrgene_analysis()
#' }
#'
EIF4F_Corrgene_analysis <- function() {
  .plot_Corr_RNAseq_TCGA_GTEX(tissue_type = "All")
  .plot_Corr_RNAseq_TCGA_GTEX(tissue_type = "Lung")

  return(NULL)
}
