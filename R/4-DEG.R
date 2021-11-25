# Differential gene expression analysis of EIF4F genes in TCGA -----------------

# This R script contains four sections.
# (1) RNAseq data preparation
# (2) analyze differential mRNA gene expression and plotting,
# (3) master functions to execute a pipeline of functions to select related RNAseq data
# with supply of EIF4F gene names for analysis and plotting.
# (4) wrapper function to call all master functions with inputs


## prepare RNA-seq related dataset from TCGA and GTEx ==========================

# due to NSE notes in R CMD check
TCGA_GTEX_RNAseq_sampletype <- NULL

#' Read all RNA-seq related datasets from TCGA and GTEx
#'
#' @description This function reads RNA-seq related datasets from TCGA and GTEx
#'
#' .TCGA_GTEX_RNAseq: the recomputed RNAseq data from both TCGA and GTEx generated from \code{\link{.get_TCGA_GTEX_RNAseq}}
#'
#' .TCGA_GTEX_sampletype: the annotation data from the "TcgaTargetGTEX_phenotype.txt"
#' dataset with selection of sample.type and primary.disease columns.
#'
#' TCGA_GTEX_RNAseq_sampletype: the merged dataset from .TCGA_GTEX_RNAseq and .TCGA_GTEX_sampletype.
#'
#' @importFrom dplyr distinct filter select select_if mutate mutate_at summarise rename group_by all_of
#' @importFrom readr read_tsv
#' @importFrom tibble remove_rownames column_to_rownames
#' @importFrom data.table fread transpose
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_RNAseq_data()
#' }
#'
initialize_RNAseq_data <- function() {
  TCGA_GTEX_RNAseq <- .get_TCGA_GTEX_RNAseq()

  TCGA_GTEX_sampletype <- readr::read_tsv(
    file.path(
      data_file_directory,
      "TcgaTargetGTEX_phenotype.txt"
    ),
    show_col_types = FALSE
  ) %>%
    # {
    as_tibble() %>%
    dplyr::distinct(.data$sample, .keep_all = TRUE) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "sample") %>%
    dplyr::select(
      "_sample_type",
      "primary disease or tissue",
      "_primary_site",
      "_study"
    ) %>%
    dplyr::rename(
      "sample.type" = "_sample_type",
      "primary.disease" = "primary disease or tissue",
      "primary.site" = "_primary_site",
      "study" = "_study"
    )
  # }

  TCGA_GTEX_RNAseq_sampletype <<- merge(TCGA_GTEX_RNAseq,
    TCGA_GTEX_sampletype,
    by    = "row.names",
    all.x = TRUE
  ) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "Row.names")
}

#' Read recomputed RNAseq data from both TCGA and GTEx
#' @description This function reads recomputed RNAseq data from both TCGA and GTEx
#' "TcgaTargetGtex_RSEM_Hugo_norm_count".
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize_RNAseq_data}} function.
#' @return a data frame that contains the recomputed RNAseq data from both TCGA and GTEx
#' @examples \dontrun{
#' .get_TCGA_GTEX_RNAseq()
#' }
#' @keywords internal
.get_TCGA_GTEX_RNAseq <- function() {
  .TCGA_pancancer <- data.table::fread(
    file.path(
      data_file_directory,
      "TcgaTargetGtex_RSEM_Hugo_norm_count"
    ),
    data.table = FALSE
  ) %>%
    as_tibble() %>%
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


## Differential expression analysis and plotting ===============================

#' Compares expressions among genes in tumor samples
#' @description This function selects the RNAseq data from TCGA samples,
#' excludes “Solid Tissue Normal” and ranks the genes by their mRNA expression.
#'
#' It should not be used directly, only inside \code{\link{.plot_boxgraph_RNAseq_TCGA}} function.
#' @param df \code{.TCGA_GTEX_RNAseq_sampletype_subset} generated inside \code{\link{.plot_boxgraph_RNAseq_TCGA}}
#' @return a data frame ranking genes by their mRNA expressions in lung adenocarcinoma
#' @importFrom dplyr pull
#' @examples \dontrun{
#' .RNAseq_all_gene(.TCGA_GTEX_RNAseq_sampletype_subset)
#' }
#' @keywords internal
#'
.RNAseq_all_gene <- function(df) {
  df <- df %>%
    dplyr::filter(.data$study == "TCGA" & .data$sample.type != "Solid Tissue Normal")
  # order the EIF genes according to the order of the expression means
  order <- df %>%
    # dplyr::filter out category equal to 'Lung Adenocarcinoma'
    dplyr::filter(.data$primary.disease == "Lung Adenocarcinoma") %>%
    # use the same groups as in the ggplot
    group_by(.data$variable) %>%
    # calculate the means
    summarise(mean_RNAseq = mean(.data$RNAseq)) %>%
    mutate(variable = fct_reorder(.data$variable, .data$mean_RNAseq)) %>%
    pull(.data$variable) %>%
    levels()
  output <- df %>%
    mutate(variable = factor(.data$variable, levels = rev(order)))
  return(output)
}

#' Grouped box plots of RNA expression across tumors
#' @description This function should not be used directly, only inside \code{\link{.plot_boxgraph_RNAseq_TCGA}}function.
#' @param data output dataset generated from \code{\link{.RNAseq_all_gene}} function
#' @keywords internal
.RNAseq_grouped_boxplot <- function(df) {
  p1 <- ggplot(
    data = df,
    aes_(
      x = ~primary.disease,
      y = ~ 2**RNAseq
    )
  ) +
    scale_y_continuous(
      trans = log2_trans(),
      # limits = c(2**4, 2**17),
      labels = label_comma()
    ) +
    stat_n_text(
      size = 5,
      fontface = "bold",
      angle = 90,
      hjust = 0
    ) +
    geom_boxplot(aes_(
      # colour = factor(variable),
      colour = ~variable,
    ),
    outlier.shape = NA,
    position = position_dodge(width = 1)
    ) +
    labs(
      x = "primary disease",
      y = paste("Normalized expression (RNA-Seq counts)")
    ) +
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = element_blank(),
      axis.title.y = black_bold_12,
      axis.text.x = black_bold_12_45,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_12,
      legend.position = c(0, 0.95),
      legend.direction = "horizontal",
      legend.justification = c(0, 1),
      legend.box = "horizontal",
      strip.text = black_bold_12
    )
  print(p1)

  ggplot2::ggsave(
    path = file.path(output_directory, "DEG"),
    filename = "RNAseqGroupedBoxplot.pdf",
    plot = p1,
    width = 16.5,
    height = 8,
    useDingbats = FALSE
  )
}


#' Analyzes differential gene expression in tumors vs adjacent normal tissues
#' @description This function selects the RNAseq data from TCGA samples including all
#' tumors and solid tissue normal samples for comparison.
#'
#' It should not be used directly, only inside \code{\link{.plot_boxgraph_RNAseq_TCGA}} function.
#' @param df \code{.TCGA_GTEX_RNAseq_sampletype_subset} generated inside \code{\link{.plot_boxgraph_RNAseq_TCGA}}
#' @param x one gene from the input argument of \code{\link{.plot_boxgraph_RNAseq_TCGA}}
#' @return a data frame of differential gene expression in tumors vs adjacent normal tissues from individual TCGA cancer types
#' @importFrom dplyr case_when
#' @importFrom forcats fct_rev
#' @keywords internal
#'
.RNAseq_ind_gene <- function(df, x) {
  df <- df %>%
    dplyr::filter(.data$study == "TCGA") %>%
    droplevels() %>%
    mutate(sample.type = case_when(
      sample.type != "Solid Tissue Normal" ~ "Tumor",
      sample.type == "Solid Tissue Normal" ~ "Normal"
    )) %>%
    dplyr::filter(.data$variable == x) %>%
    mutate(primary.disease = forcats::fct_rev(.data$primary.disease))
  output <- list(df, x)
  return(output)
}

#' Box plots of differential gene expression across tumors
#' @description This function should not be used directly, only inside \code{\link{.plot_boxgraph_RNAseq_TCGA}}function.
#' @param data output dataset generated from \code{\link{.RNAseq_ind_gene}} function
#' @importFrom scales log2_trans label_comma
#' @keywords internal
.RNAseq_boxplot <- function(df) {
  p1 <- ggplot(
    data = df[[1]],
    aes_(
      x = ~primary.disease,
      y = ~ 2**RNAseq,
      color = ~sample.type
    )
  ) +
    scale_y_continuous(
      trans = scales::log2_trans(),
      # limits = c(2**11, 2**17),# for 4g
      # limits = c(2**7, 2**14),# for eif4E
      labels = scales::label_comma()
    ) +
    stat_n_text(
      size = 5,
      fontface = "bold",
      hjust = 0
    ) +
    geom_boxplot(
      # alpha = .1,
      # fill = sample.type,
      outlier.shape = NA,
      # size = .75,
      # width = 1,
      position = "dodge"
    ) +
    scale_color_manual(
      values = c(
        "Tumor" = "#CC79A7",
        "Normal" = "#0072B2"
      ),
      breaks = c("Tumor", "Normal"),
      labels = c("Tumor\n", "Normal\n")
    ) +
    labs(
      x = "primary disease",
      y = paste(df[[2]], "expression (RNA-Seq counts)")
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
    )
  # geom_signif(comparisons=list(c("Tumor", "Normal")))
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "DEG"),
    filename = paste0(df[[2]], "tumorvsnormal.pdf"),
    plot = p1,
    width = 7.5,
    height = 9,
    useDingbats = FALSE
  )
}


#' Analyzes differential gene expression in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined
#' @description This function selects the RNAseq data from TCGA samples including
#' tumor samples that are labeled as “metastatic”, “primary tumor”, and “solid tissue normal” for comparison..
#'
#' It should not be used directly, only inside \code{\link{.plot_boxgraph_RNAseq_TCGA}} function.
#' @param df \code{.TCGA_GTEX_RNAseq_sampletype_subset} generated inside \code{\link{.plot_boxgraph_RNAseq_TCGA}}
#' @return a data frame of differential gene expression in tumors vs adjacent normal tissues in individual TCGA cancer types
#' @keywords internal
#'
.RNAseq_tumortype <- function(df) {
  df1 <- df %>%
    dplyr::filter(.data$study == "TCGA") %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::filter(.data$sample.type %in% c(
      "Metastatic",
      "Primary Tumor",
      "Solid Tissue Normal"
    ))
  return(df1)
}

#' Violin plots of differential gene expression/ratio in primary, metastatic tumors vs adjacent normal tissues
#' @description This function should not be used directly, only inside \code{\link{.plot_boxgraph_RNAseq_TCGA}} or \code{\link{.plot_boxgraph_RNAratio_TCGA}} function.
#' @param data output dataset generated from \code{\link{.RNAseq_tumortype}} or \code{\link{.RNAratio_tumortype}} function
#' @importFrom EnvStats stat_n_text
#' @keywords internal
.violinplot <- function(df, y.axis.title, y.axis.break, yintercept) {
  p1 <- ggplot(
    data = df,
    # data = EIF.TCGA.RNAseq.anno.subset.long,
    aes_(
      x = ~sample.type,
      y = ~ 2**RNAseq,
      color = ~sample.type,
      fill = ~sample.type
    )
  ) +
    EnvStats::stat_n_text(
      size = 6,
      fontface = "bold",
      angle = 90,
      hjust = 0
    ) +
    ggplot2::facet_grid(. ~ variable,
      scales = "free",
      space  = "free"
    ) +
    # facet_wrap(~ variable, ncol = 6) +
    geom_violin(trim = TRUE) +
    geom_boxplot(
      alpha = .01,
      width = .25,
      color = "black",
      position = position_dodge(width = .9)
    ) +
    labs(
      x = "sample type",
      y = y.axis.title
      # y = "normalized RNA counts"
    ) +
    scale_x_discrete(labels = c(
      "Metastatic Tumor",
      "Primary Tumor",
      # "Normal Tissue",
      "Adjacent Normal"
    )) +
    scale_y_continuous(
      trans = scales::log2_trans(),
      labels = scales::label_comma(),
      breaks = y.axis.break,
      # breaks = c(1, 128, 2048, 32768),
      # labels = c("1","8","64","512"),
      position = "left"
    ) +
    geom_hline(yintercept = yintercept, linetype = "dashed") +
    scale_color_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + # for color-blind palettes
    scale_fill_manual(values = c("#56B4E9", "#009E73", "#D55E00")) + # for color-blind palettes
    theme_bw() +
    theme(
      plot.title = black_bold_16,
      axis.title.x = element_blank(),
      axis.title.y = black_bold_16_right,
      axis.text.x = black_bold_16_90,
      axis.text.y = black_bold_16_90,
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = black_bold_16
    ) +
    ggpubr::stat_compare_means(
      comparisons = list(
        c("Metastatic", "Solid Tissue Normal"),
        c("Primary Tumor", "Solid Tissue Normal"),
        c("Metastatic", "Primary Tumor")
      ),
      method = "t.test",
      label = "p.signif",
      size = 6
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "DEG"),
    filename = "EIF4Fviolin.pdf",
    plot = p1,
    width = 18,
    height = 9,
    useDingbats = FALSE
  )
}


#' Calculate RNA ratios in tumors vs adjacent normal tissues
#' @description This function calculates the RNA ratio data from TCGA samples including all
#' tumors and solid tissue normal samples for comparison.
#'
#' It should not be used directly, only inside \code{\link{.plot_boxgraph_RNAratio_TCGA}} function.
#' @param EIF4E gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4E2 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4E3 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4EBP1 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4G1 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4G2 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4G3 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF3D gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4A1 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @param EIF4A2 gene name, same input from \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @importFrom rlang :=
#' @return a data frame of differential RNA ratios in tumors vs adjacent normal tissues from individual TCGA cancer types
#' @keywords internal
#'
.RNAratio_calculation <- function(EIF4E, EIF4E2, EIF4E3, EIF4EBP1,
                                  EIF4G1, EIF4G2, EIF4G3, EIF3D,
                                  EIF4A1, EIF4A2) {
  .genes_names <- c(
    EIF4E, EIF4E2, EIF4E3, EIF4EBP1,
    EIF4G1, EIF4G2, EIF4G3, EIF3D,
    EIF4A1, EIF4A2
  )
  .RNAratio_data <- TCGA_GTEX_RNAseq_sampletype %>%
    dplyr::select(
      all_of(.genes_names),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    as_tibble() %>%
    filter(EIF4E != 0 & !is.na(.data$primary.site)) %>%
    # calculate the ratio of mRNA counts
    dplyr::mutate(
      (!!paste0(EIF4E, "+", EIF4EBP1)) := log2(2**(!!as.name(EIF4E)) + 2**(!!as.name(EIF4EBP1)) - 1),
      (!!paste0(EIF4A1, ":", "\n", EIF4E)) := (!!as.name(EIF4A1)) - (!!as.name(EIF4E)),
      (!!paste0(EIF4A1, ":", "\n", EIF4E2)) := (!!as.name(EIF4A1)) - (!!as.name(EIF4E2)),
      (!!paste0(EIF4A2, ":", "\n", EIF4E)) := (!!as.name(EIF4A2)) - (!!as.name(EIF4E)),
      (!!paste0(EIF4A2, ":", "\n", EIF4E2)) := (!!as.name(EIF4A2)) - (!!as.name(EIF4E2)),
      (!!paste0(EIF4G1, ":", "\n", EIF4E)) := (!!as.name(EIF4G1)) - (!!as.name(EIF4E)),
      (!!paste0(EIF4G1, ":", "\n", EIF4E2)) := (!!as.name(EIF4G1)) - (!!as.name(EIF4E2)),
      (!!paste0(EIF4G3, ":", "\n", EIF4E)) := (!!as.name(EIF4G3)) - (!!as.name(EIF4E)),
      (!!paste0(EIF4G3, ":", "\n", EIF4E2)) := (!!as.name(EIF4G3)) - (!!as.name(EIF4E2)),
      (!!paste0(EIF4A1, ":", "\n", EIF4G1)) := (!!as.name(EIF4A1)) - (!!as.name(EIF4G1)),
      (!!paste0(EIF4A2, ":", "\n", EIF4G1)) := (!!as.name(EIF4A2)) - (!!as.name(EIF4G1)),
      (!!paste0(EIF4A1, ":", "\n", EIF4G2)) := (!!as.name(EIF4A1)) - (!!as.name(EIF4G2)),
      (!!paste0(EIF4A2, ":", "\n", EIF4G2)) := (!!as.name(EIF4A2)) - (!!as.name(EIF4G2)),
      (!!paste0(EIF4E, ":", "\n", EIF4EBP1)) := (!!as.name(EIF4E)) - (!!as.name(EIF4EBP1)),
      (!!paste0(EIF4E2, ":", "\n", EIF4E)) := (!!as.name(EIF4E2)) - (!!as.name(EIF4E)),
      (!!paste0(EIF4G2, ":", "\n", EIF4G1)) := (!!as.name(EIF4G2)) - (!!as.name(EIF4G1)),
      (!!paste0(EIF4G1, ":", "\n", EIF4G3)) := (!!as.name(EIF4G1)) - (!!as.name(EIF4G3)),
      (!!paste0(EIF4A1, ":", "\n", EIF4A2)) := (!!as.name(EIF4A1)) - (!!as.name(EIF4A2)),
      (!!paste0(EIF4G1, ":", "\n", EIF4E, "+", EIF4EBP1)) := (!!as.name(EIF4G1)) - log2(2**(!!as.name(EIF4E)) + 2**(!!as.name(EIF4EBP1)) - 1),
      (!!paste0(EIF4A1, ":", "\n", EIF4E, "+", EIF4EBP1)) := (!!as.name(EIF4A1)) - log2(2**(!!as.name(EIF4E)) + 2**(!!as.name(EIF4EBP1)) - 1)
    ) %>%
    dplyr::select(
      (!!paste0(EIF4A1, ":", "\n", EIF4E)), (!!paste0(EIF4A1, ":", "\n", EIF4E2)),
      (!!paste0(EIF4A2, ":", "\n", EIF4E)), (!!paste0(EIF4A2, ":", "\n", EIF4E2)),
      (!!paste0(EIF4G1, ":", "\n", EIF4E)), (!!paste0(EIF4G1, ":", "\n", EIF4E2)),
      (!!paste0(EIF4G3, ":", "\n", EIF4E)), (!!paste0(EIF4G3, ":", "\n", EIF4E2)),
      (!!paste0(EIF4A1, ":", "\n", EIF4G1)), (!!paste0(EIF4A2, ":", "\n", EIF4G1)),
      (!!paste0(EIF4A1, ":", "\n", EIF4G2)), (!!paste0(EIF4A2, ":", "\n", EIF4G2)),
      (!!paste0(EIF4E, ":", "\n", EIF4EBP1)), (!!paste0(EIF4E2, ":", "\n", EIF4E)),
      (!!paste0(EIF4G2, ":", "\n", EIF4G1)), (!!paste0(EIF4G1, ":", "\n", EIF4G3)),
      (!!paste0(EIF4A1, ":", "\n", EIF4A2)),
      (!!paste0(EIF4G1, ":", "\n", EIF4E, "+", EIF4EBP1)),
      (!!paste0(EIF4A1, ":", "\n", EIF4E, "+", EIF4EBP1)),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    na.omit()
  return(.RNAratio_data)
}

#' Select RNA ratio data for plotting
#' @description This function should not be used directly, only inside \code{\link{.plot_boxgraph_RNAratio_TCGA}}function.
#' @param df output dataset generated from \code{\link{.RNAratio_calculation}} function
#' @param x input argument for selected variables generated from \code{\link{.plot_boxgraph_RNAratio_TCGA}} function
#' @keywords internal
.RNAratio_selection <- function(df, x) {
  .RNAratio_data <- df %>%
    dplyr::filter(.data$study == "TCGA") %>%
    droplevels() %>%
    dplyr::mutate(sample.type = case_when(
      sample.type != "Solid Tissue Normal" ~ "Tumor",
      sample.type == "Solid Tissue Normal" ~ "NAT"
    )) %>%
    dplyr::select(
      all_of(x),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    melt(
      id = c(
        "sample.type",
        "primary.disease",
        "primary.site",
        "study"
      ),
      value.name = "RNAseq"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(primary.disease = forcats::fct_rev(.data$primary.disease))
  return(.RNAratio_data)
}



#' Box plots of differential RNA ratios across tumors
#' @description This function should not be used directly, only inside \code{\link{.plot_boxgraph_RNAratio_TCGA}}function.
#' @param data output dataset generated from \code{\link{.RNAratio_selection}} function
#' @keywords internal
.RNAratio_boxplot <- function(df, dashline, ylimit, filename) {
  p1 <- ggplot(
    data = df,
    aes_(
      x = ~primary.disease,
      # x = f.ordered1,
      y = ~ 2**RNAseq,
      # fill  = variable,
      color = ~sample.type
    )
  ) +
    geom_boxplot(
      # alpha = .1,
      # fill = sample.type,
      outlier.shape = NA,
      # size = .75,
      # width = 1,
      position = "dodge"
    ) +
    geom_hline(
      # data = df,
      aes(yintercept = dashline),
      linetype = "dashed"
    ) +
    scale_color_manual(
      values = c("Tumor" = "#CC79A7", "NAT" = "#0072B2"),
      breaks = c("Tumor", "NAT")
    ) +
    # for color-blind palettes
    facet_wrap(~variable, scales = "free_x") +
    ggplot2::facet_grid(~variable, scales = "free_x", space = "free") +
    guides(colour = guide_legend(nrow = 1)) +
    labs(
      x = "primary disease",
      y = "Ratio of RNA counts"
    ) +
    coord_flip(ylim = ylimit) +
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = black_bold_12,
      axis.title.y = element_blank(),
      axis.text.x = black_bold_12,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      legend.text = black_bold_12,
      # strip.background = element_blank(),
      strip.text.x = black_bold_12
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output_directory, "DEG"),
    filename = filename,
    plot = p1,
    width = 18,
    height = 8,
    useDingbats = FALSE
  )
}

#' Analyzes differential RNA ratios in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined
#' @description This function selects the RNA ratio data from TCGA samples including
#' tumor samples that are labeled as “metastatic”, “primary tumor”, and “solid tissue normal” for comparison..
#'
#' It should not be used directly, only inside \code{\link{.plot_boxgraph_RNAratio_TCGA}} function.
#' @param df \code{.RNAratio_tumortype(.RNAratio_data)} generated inside \code{\link{.plot_boxgraph_RNAratio_TCGA}}
#' @return a data frame of differential RNA ratios in tumors vs adjacent normal tissues from individual TCGA cancer types
#' @keywords internal
#'
.RNAratio_tumortype <- function(df, x) {
  .RNAratio_data <- df %>%
    dplyr::filter(.data$sample.type %in% c(
      "Metastatic",
      "Primary Tumor",
      "Solid Tissue Normal"
    )) %>%
    droplevels() %>%
    dplyr::select(
      all_of(x),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    melt(
      id = c(
        "sample.type",
        "primary.disease",
        "primary.site",
        "study"
      ),
      value.name = "RNAseq"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(primary.disease = forcats::fct_rev(.data$primary.disease))
  return(.RNAratio_data)
}


## master functions to call DEG gene analysis and plotting =====================

#' Compare the expression of EIF4F genes
#' @description This function generates a summary box plot to compare the expression of all EIF4F genes across TCGA cancer types,
#' box plots for differential gene expression in tumors tumors vs adjacent normal tissues for each  gene.
#' and a violin plot to compare differential gene expression in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined.
#' @param EIF.list gene names in a vector of characters
#' @details  This function first selects RNASeq and sample type data from input gene
#' in the data frame \code{TCGA_GTEX_RNAseq_sampletype} prepared from \code{\link{initialize_RNAseq_data}}.
#'
#' With the subset data \code{.TCGA_GTEX_RNAseq_sampletype_subset}, it compares the mRNA expression of all inquired genes
#' with the functions \code{\link{.RNAseq_all_gene}} and plot the results as a bar plot with \code{\link{.RNAseq_grouped_boxplot}}
#'
#' It also performs differential expression analysis for each gene across
#' all tumors types with the function \code{\link{.RNAseq_ind_gene}} and plots with \code{\link{.RNAseq_boxplot}}.
#'
#' Finally, it performs differential expression analysis for each gene
#' in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined with
#' the function \code{\link{.RNAseq_tumortype}} and plots with \code{\link{.violinplot}}.
#' It should not be used directly, only inside \code{\link{EIF4F_DEG_analysis}} function.
#' @return box plots for differential gene expression of \code{EIF.list} in TCGA tumors
#' @keywords internal
#' @examples \dontrun{
#' .plot_boxgraph_RNAseq_TCGA(c(
#'   "EIF4G1", "EIF4G2", "EIF4G3",
#'   "PABPC1", "EIF4A1", "EIF4A2", "EIF4B", "EIF4H", "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1", "EIF3D"
#' ))
#' }
.plot_boxgraph_RNAseq_TCGA <- function(EIF.list) {
  .TCGA_GTEX_RNAseq_sampletype_subset <- TCGA_GTEX_RNAseq_sampletype %>%
    dplyr::select(
      all_of(EIF.list),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    as_tibble() %>%
    melt(
      id = c(
        "sample.type",
        "primary.disease",
        "primary.site",
        "study"
      ),
      value.name = "RNAseq"
    ) %>%
    filter(!is.na(.data$primary.site)) %>%
    # na.omit(.$primary.site) %>%
    # filter(RNAseq != 0) %>%
    mutate_if(is.character, as.factor)

  # boxplot to compare relative abundance of genes across tumors
  .RNAseq_all_gene(.TCGA_GTEX_RNAseq_sampletype_subset) %>%
    .RNAseq_grouped_boxplot()

  # boxplot to compare RNA-seq of one gene in tumor vs adjacent normal
  RNAseq.ind.gene.df <- lapply(EIF.list, .RNAseq_ind_gene,
    df = .TCGA_GTEX_RNAseq_sampletype_subset
  )
  lapply(RNAseq.ind.gene.df, .RNAseq_boxplot)

  # violin plot to compare  expression in primary, metastatic tumors vs NATs
  .RNAseq_tumortype(.TCGA_GTEX_RNAseq_sampletype_subset) %>%
    .violinplot(
      y.axis.title = "normalized RNA counts",
      y.axis.break = c(128, 2048, 32768),
      yintercept = NULL
    )
}

#' Compare the RNA ratios between EIF4F genes
#' @description This function generates box plot to compare the RNA ratios between EIF4F genes in tumors tumors vs adjacent normal tissues across TCGA cancer types.
#' and a violin plot to compare RNA ratios in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined.
#' @param EIF4E gene name as a string
#' @param EIF4E2 gene name as a string
#' @param EIF4E3 gene name as a string
#' @param EIF4EBP1 gene name as a string
#' @param EIF4G1 gene name as a string
#' @param EIF4G2 gene name as a string
#' @param EIF4G3 gene name as a string
#' @param EIF3D gene name as a string
#' @param EIF4A1 gene name as a string
#' @param EIF4A2 gene name as a string
#' @details  This function analyzes the ratios of mRNA levels between input genes within each TCGA sample,
#' in the data frame \code{TCGA_GTEX_RNAseq_sampletype} prepared from \code{\link{initialize_RNAseq_data}}.
#'
#' It calculates RNA ratio across individual tumor types from all cancers
#' with the function \code{\link{.RNAratio_calculation}}, select the ratio data with \code{\link{.RNAratio_selection}}
#' and plots with \code{\link{.RNAratio_boxplot}}.
#'
#' Finally, it compares RNA ratio in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined with
#' the function \code{\link{.RNAratio_tumortype}} and plots with \code{\link{.violinplot}}.
#' It should not be used directly, only inside \code{\link{EIF4F_DEG_analysis}} function.
#' @return box plots for RNA ratios among input argument in TCGA tumors
#' @examples \dontrun{
#' .plot_boxgraph_RNAratio_TCGA(
#'   EIF4E = "EIF4E", EIF4E2 = "EIF4E2",
#'   EIF4E3 = "EIF4E3", EIF4EBP1 = "EIF4EBP1", EIF4G1 = "EIF4G1", EIF4G2 = "EIF4G2",
#'   EIF4G3 = "EIF4G3", EIF3D = "EIF3D", EIF4A1 = "EIF4A1", EIF4A2 = "EIF4A2"
#' )
#' }
.plot_boxgraph_RNAratio_TCGA <- function(EIF4E, EIF4E2, EIF4E3, EIF4EBP1,
                                         EIF4G1, EIF4G2, EIF4G3, EIF3D,
                                         EIF4A1, EIF4A2) {
  RNAratio.data <- .RNAratio_calculation(
    EIF4E, EIF4E2, EIF4E3, EIF4EBP1,
    EIF4G1, EIF4G2, EIF4G3, EIF3D,
    EIF4A1, EIF4A2
  )
  .RNAratio_boxplot(
    df = .RNAratio_selection(RNAratio.data, c(
      (paste0(EIF4G1, ":", "\n", EIF4E)), (paste0(EIF4A1, ":", "\n", EIF4E)),
      (paste0(EIF4A2, ":", "\n", EIF4E)), (paste0(EIF4G3, ":", "\n", EIF4E)),
      (paste0(EIF4G3, ":", "\n", EIF4E2)), (paste0(EIF4G1, ":", "\n", EIF4G3))
    )),
    dashline = 1,
    ylimit = c(0, 25),
    filename = "RNAratio1.pdf"
  )

  .RNAratio_boxplot(
    df = .RNAratio_selection(RNAratio.data, c(
      (paste0(EIF4G2, ":", "\n", EIF4G1)), (paste0(EIF4E2, ":", "\n", EIF4E)),
      (paste0(EIF4A1, ":", "\n", EIF4A2)), (paste0(EIF4E, ":", "\n", EIF4EBP1)),
      (paste0(EIF4G1, ":", "\n", EIF4E, "+", EIF4EBP1)),
      (paste0(EIF4A1, ":", "\n", EIF4E, "+", EIF4EBP1))
    )),
    dashline = 4,
    ylimit = c(0, 25),
    filename = "RNAratio2.pdf"
  )

  .RNAratio_boxplot(
    df = .RNAratio_selection(RNAratio.data, c(
      (paste0(EIF4G3, ":", "\n", EIF4E)), (paste0(EIF4G3, ":", "\n", EIF4E2)),
      (paste0(EIF4G2, ":", "\n", EIF4G1)), (paste0(EIF4E2, ":", "\n", EIF4E)),
      (paste0(EIF4A1, ":", "\n", EIF4A2)), (paste0(EIF4E, ":", "\n", EIF4EBP1))
    )),
    dashline = 1,
    ylimit = c(0, 5),
    filename = "RNAratio3.pdf"
  )

  .RNAratio_tumortype(RNAratio.data, c(
    (paste0(EIF4G1, ":", "\n", EIF4E)), (paste0(EIF4A1, ":", "\n", EIF4E)),
    (paste0(EIF4A2, ":", "\n", EIF4E)), (paste0(EIF4G3, ":", "\n", EIF4E)),
    (paste0(EIF4G3, ":", "\n", EIF4E2)), (paste0(EIF4G1, ":", "\n", EIF4G3)),
    (paste0(EIF4A1, ":", "\n", EIF4A2)), (paste0(EIF4E, ":", "\n", EIF4EBP1)),
    (paste0(EIF4G1, ":", "\n", EIF4E, "+", EIF4EBP1)),
    (paste0(EIF4A1, ":", "\n", EIF4E, "+", EIF4EBP1))
    # "EIF4G1:\nEIF4E", "EIF4A1:\nEIF4E",
    # "EIF4A2:\nEIF4E", "EIF4G3:\nEIF4E",
    # "EIF4G3:\nEIF4E2", "EIF4G1:\nEIF4G3",
    # "EIF4G2:\nEIF4G1", "EIF4E2:\nEIF4E",
    # "EIF4A1:\nEIF4A2", "EIF4E:\nEIF4EBP1",
    # "EIF4G1:\nEIF4E+EIF4EBP1",
    # "EIF4A1:\nEIF4E+EIF4EBP1"
  )) %>%
    .violinplot(
      y.axis.title = "Ratio of RNA counts",
      # y.axis.break = c(1, 128, 2048, 32768)
      y.axis.break = c(0.125, 1, 4, 8, 64, 512),
      yintercept = c(1, 4)
    )
}

## wrapper function to call all master functions with inputs ===================

#' Perform all CNV related analysis and generate plots
#' @description A wrapper function to call all master functions for CNV analysis with inputs
#'
#' @details  This function run three master functions together:
#' \code{\link{.plot_boxgraph_RNAseq_TCGA}} and
#' \code{\link{.plot_boxgraph_RNAratio_TCGA}} with inputs
#'
#' @return DEG analysis plots
#'
#' @export
#'
#' @examples \dontrun{
#' EIF4F_DEG_analysis()
#' }
EIF4F_DEG_analysis <- function() {
  .plot_boxgraph_RNAseq_TCGA(c(
    "EIF4G1", "EIF4G2", "EIF4G3", "PABPC1", "EIF4A1", "EIF4A2", "EIF4B", "EIF4H",
    "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1", "EIF3D"
  ))

  .plot_boxgraph_RNAratio_TCGA(
    EIF4E = "EIF4E", EIF4E2 = "EIF4E2", EIF4E3 = "EIF4E3", EIF4EBP1 = "EIF4EBP1",
    EIF4G1 = "EIF4G1", EIF4G2 = "EIF4G2", EIF4G3 = "EIF4G3", EIF3D = "EIF3D",
    EIF4A1 = "EIF4A1", EIF4A2 = "EIF4A2"
  )
}
