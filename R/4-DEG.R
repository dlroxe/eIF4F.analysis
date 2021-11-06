# prepare RNA-seq related dataset from TCGA and GTEx----------------------------
#' Read all RNA-seq related datasets from TCGA and GTEx
#'
#' @description This function reads RNA-seq related datasets from TCGA and GTEx
#'
#' TCGA.GTEX.RNAseq: the recomputed RNAseq data from both TCGA and GTEx generated from \code{\link{.get.TCGA.GTEX.RNAseq}}
#'
#' TCGA.GTEX.sampletype: the annotation data from the "TcgaTargetGTEX_phenotype.txt"
#' dataset with selection of sample.type and primary.disease columns.
#'
#' TCGA.GTEX.RNAseq.sampletype: the merged dataset from TCGA.GTEX.RNAseq and TCGA.GTEX.sampletype.
#'
#' @importFrom dplyr distinct filter select select_if mutate mutate_at summarise rename group_by all_of
#' @importFrom readr read_tsv
#' @importFrom tibble remove_rownames column_to_rownames
#' @importFrom data.table fread transpose
#'
#' @export
#'
#' @examples \dontrun{initialize.RNAseq.data()}
#'
initialize.RNAseq.data <- function() {
  TCGA.GTEX.RNAseq <- .get.TCGA.GTEX.RNAseq()

  TCGA.GTEX.sampletype <- readr::read_tsv(
    file.path(
      data.file.directory,
      "TcgaTargetGTEX_phenotype.txt"
    ),
    show_col_types = FALSE
  ) %>%
    {
      as_tibble(.) %>%
        dplyr::distinct(., sample, .keep_all = TRUE) %>%
        tibble::remove_rownames(.) %>%
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
    }

  TCGA.GTEX.RNAseq.sampletype <<- merge(TCGA.GTEX.RNAseq,
    TCGA.GTEX.sampletype,
    by    = "row.names",
    all.x = TRUE
  ) %>%
    {
      tibble::remove_rownames(.) %>%
        tibble::column_to_rownames(var = "Row.names")
    }
}

#' Read recomputed RNAseq data from both TCGA and GTEx
#' @description This function reads recomputed RNAseq data from both TCGA and GTEx
#' "TcgaTargetGtex_RSEM_Hugo_norm_count".
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize.RNAseq.data}} function.
#' @return a data frame that contains the recomputed RNAseq data from both TCGA and GTEx
#' @examples \dontrun{.get.TCGA.GTEX.RNAseq()}
#' @keywords internal
.get.TCGA.GTEX.RNAseq <- function() {
  TCGA.pancancer <- data.table::fread(
    file.path(
      data.file.directory,
      "TcgaTargetGtex_RSEM_Hugo_norm_count"
    ),
    data.table = FALSE
  ) %>%
    as_tibble(.) %>%
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


# Differential expression analysis and plotting --------------------------------

#' Compares expressions among genes in tumor samples
#' @description This function selects the RNAseq data from TCGA samples,
#' excludes “Solid Tissue Normal” and ranks the genes by their mRNA expression.
#'
#' It should not be used directly, only inside \code{\link{plot.boxgraph.RNAseq.TCGA}} function.
#' @param df \code{TCGA.GTEX.RNAseq.sampletype.subset} generated inside \code{\link{plot.boxgraph.RNAseq.TCGA}}
#' @return a data frame ranking genes by their mRNA expressions in lung adenocarcinoma
#' @examples \dontrun{.RNAseq.all.gene(TCGA.GTEX.RNAseq.sampletype.subset)}
#' @keywords internal
#'
.RNAseq.all.gene <- function(df) {
  df <- df %>%
    dplyr::filter(study == "TCGA" & sample.type != "Solid Tissue Normal")
  # order the EIF genes according to the order of the expression means
  order <- df %>%
    # dplyr::filter out category equal to 'Lung Adenocarcinoma'
    dplyr::filter(primary.disease == "Lung Adenocarcinoma") %>%
    # use the same groups as in the ggplot
    group_by(variable) %>%
    # calculate the means
    summarise(mean_RNAseq = mean(RNAseq)) %>%
    mutate(variable = fct_reorder(variable, mean_RNAseq)) %>%
    .$variable %>%
    levels()
  output <- df %>%
    mutate(variable = factor(variable, levels = rev(order)))
  return(output)
}

#' Grouped box plots of RNA expression across tumors
#' @description This function should not be used directly, only inside \code{\link{plot.boxgraph.RNAseq.TCGA}}function.
#' @param data output dataset generated from \code{\link{.RNAseq.all.gene}} function
#' @keywords internal
.RNAseq.grouped.boxplot <- function(df) {
  p1 <- ggplot(
    data = df,
    aes(
      x = primary.disease,
      y = 2**RNAseq
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
    geom_boxplot(aes(
      colour = factor(variable),
      # fill = factor(variable)
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
    path = file.path(output.directory, "Expression"),
    filename = "tumorexpression.pdf",
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
#' It should not be used directly, only inside \code{\link{plot.boxgraph.RNAseq.TCGA}} function.
#' @param df \code{TCGA.GTEX.RNAseq.sampletype.subset} generated inside \code{\link{plot.boxgraph.RNAseq.TCGA}}
#' @param x one gene from the input argument of \code{\link{plot.boxgraph.RNAseq.TCGA}}
#' @return a data frame of differential gene expression in tumors vs adjacent normal tissues from individual TCGA cancer types
#' @importFrom dplyr case_when
#' @importFrom forcats fct_rev
#' @keywords internal
#'
.RNAseq.ind.gene <- function(df, x) {
  df <- df %>%
    dplyr::filter(study == "TCGA") %>%
    droplevels() %>%
    mutate(sample.type = case_when(
      sample.type != "Solid Tissue Normal" ~ "Tumor",
      sample.type == "Solid Tissue Normal" ~ "Normal"
    )) %>%
    dplyr::filter(variable == x) %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  output <- list(df, x)
  return(output)
}

#' Box plots of differential gene expression across tumors
#' @description This function should not be used directly, only inside \code{\link{plot.boxgraph.RNAseq.TCGA}}function.
#' @param data output dataset generated from \code{\link{.RNAseq.ind.gene}} function
#' @keywords internal
.RNAseq.boxplot <- function(df) {
  p1 <- ggplot(
    data = df[[1]],
    aes(
      x = primary.disease,
      y = 2**RNAseq,
      color = sample.type
    )
  ) +
    scale_y_continuous(
      trans = log2_trans(),
      # limits = c(2**11, 2**17),# for 4g
      # limits = c(2**7, 2**14),# for eif4E
      labels = label_comma()
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
    path = file.path(output.directory, "Expression"),
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
#' It should not be used directly, only inside \code{\link{plot.boxgraph.RNAseq.TCGA}} function.
#' @param df \code{TCGA.GTEX.RNAseq.sampletype.subset} generated inside \code{\link{plot.boxgraph.RNAseq.TCGA}}
#' @return a data frame of differential gene expression in tumors vs adjacent normal tissues in individual TCGA cancer types
#' @keywords internal
#'
.RNAseq.tumortype <- function(df) {
  TCGA.GTEX.RNAseq.sampletype.subset <- df %>%
    dplyr::filter(study == "TCGA") %>%
    mutate_if(is.character, as.factor) %>%
    dplyr::filter(sample.type %in% c(
      "Metastatic",
      "Primary Tumor",
      "Solid Tissue Normal"
    ))
  return(TCGA.GTEX.RNAseq.sampletype.subset)
}

#' Violin plots of differential gene expression/ratio in primary, metastatic tumors vs adjacent normal tissues
#' @description This function should not be used directly, only inside \code{\link{plot.boxgraph.RNAseq.TCGA}} or \code{\link{plot.boxgraph.RNAratio.TCGA}} function.
#' @param data output dataset generated from \code{\link{.RNAseq.tumortype}} or \code{\link{.RNAratio.tumortype}} function
#' @keywords internal
.violinplot <- function(df) {
  p1 <- ggplot(
    data = df,
    # data = EIF.TCGA.RNAseq.anno.subset.long,
    aes(
      x = sample.type,
      y = 2**RNAseq,
      color = sample.type,
      fill = sample.type
    )
  ) +
    stat_n_text(
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
      y = "normalized RNA counts"
    ) +
    scale_x_discrete(labels = c(
      "Metastatic Tumor",
      "Primary Tumor",
      # "Normal Tissue",
      "Adjacent Normal"
    )) +
    scale_y_continuous(
      trans = log2_trans(),
      labels = label_comma(),
      breaks = c(1, 128, 2048, 32768),
      # labels = c("1","8","64","512"),
      position = "left"
    ) +
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
    stat_compare_means(
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
    path = file.path(output.directory, "Expression"),
    filename = "EIFexpressionviolin.pdf",
    plot = p1,
    width = 18,
    height = 9,
    useDingbats = FALSE
  )
}


#' Analyzes RNA ratios in tumors vs adjacent normal tissues
#' @description This function selects the RNA ratio data from TCGA samples including all
#' tumors and solid tissue normal samples for comparison.
#'
#' It should not be used directly, only inside \code{\link{plot.boxgraph.RNAratio.TCGA}} function.
#' @param df \code{.RNAratio.ind(RNAratio.data)} generated inside \code{\link{plot.boxgraph.RNAratio.TCGA}}
#' @return a data frame of differential RNA ratios in tumors vs adjacent normal tissues from individual TCGA cancer types
#' @keywords internal
#'
.RNAratio.ind <- function(df) {
  RNAratio.data <- df %>%
    dplyr::filter(study == "TCGA") %>%
    droplevels() %>%
    mutate(sample.type = case_when(
      sample.type != "Solid Tissue Normal" ~ "Tumor",
      sample.type == "Solid Tissue Normal" ~ "NAT"
    )) %>%
    melt(.,
      id = c(
        "sample.type",
        "primary.disease",
        "primary.site",
        "study"
      ),
      value.name = "RNAseq"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  return(RNAratio.data)
}

#' Box plots of differential RNA ratios across tumors
#' @description This function should not be used directly, only inside \code{\link{plot.boxgraph.RNAratio.TCGA}}function.
#' @param data output dataset generated from \code{\link{.RNAratio.ind}} function
#' @keywords internal
.RNAratio.boxplot <- function(df, dashline, ylimit, filename) {
  p1 <- ggplot(
    data = df,
    aes(
      x = primary.disease,
      # x = f.ordered1,
      y = 2**RNAseq,
      # fill  = variable,
      color = sample.type
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
    path = file.path(output.directory, "Expression"),
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
#' It should not be used directly, only inside \code{\link{plot.boxgraph.RNAratio.TCGA}} function.
#' @param df \code{.RNAratio.tumortype(RNAratio.data)} generated inside \code{\link{plot.boxgraph.RNAratio.TCGA}}
#' @return a data frame of differential RNA ratios in tumors vs adjacent normal tissues from individual TCGA cancer types
#' @keywords internal
#'
.RNAratio.tumortype <- function(df) {
  RNAratio.data <- df %>%
    dplyr::filter(sample.type %in% c(
      "Metastatic",
      "Primary Tumor",
      "Solid Tissue Normal"
    )) %>%
    droplevels() %>%
    melt(.,
         id = c(
           "sample.type",
           "primary.disease",
           "primary.site",
           "study"
         ),
         value.name = "RNAseq"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    mutate(primary.disease = forcats::fct_rev(primary.disease))
  return(RNAratio.data)
}

# master functions to call DEG gene analysis and plotting ----------------------
#' Compare the expression of EIF4F genes
#' @description This function generates a summary box plot to compare the expression of all EIF4F genes across TCGA cancer types,
#' box plots for differential gene expression in tumors tumors vs adjacent normal tissues for each  gene.
#' and a violin plot to compare differential gene expression in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined.
#' @param EIF gene names in a vector of characters
#' @details  This function first selects RNASeq and sample type data from input gene
#' in the data frame \code{TCGA.GTEX.RNAseq.sampletype} prepared from \code{\link{initialize.RNAseq.data}}.
#'
#' With the subset data \code{TCGA.GTEX.RNAseq.sampletype.subset}, it compares the mRNA expression of all inquired genes
#' with the functions \code{\link{.RNAseq.all.gene}} and plot the results as a bar plot with \code{\link{.RNAseq.grouped.boxplot}}
#'
#' It also performs differential expression analysis for each gene across
#' all tumors types with the function \code{\link{.RNAseq.ind.gene}} and plots with \code{\link{.RNAseq.boxplot}}.
#'
#' Finally, it performs differential expression analysis for each gene
#' in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined with
#' the function \code{\link{.RNAseq.tumortype}} and plots with \code{\link{.violinplot}}.
#' @return box plots for diferenital gene expression of \code{EIF} in TCGA tumors
#' @export
#' @examples \dontrun{plot.boxgraph.RNAseq.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))}
#' @examples \dontrun{plot.boxgraph.RNAseq.TCGA(c("EIF4G1", "EIF4G2", "EIF4G3",
#' "PABPC1", "EIF4A1", "EIF4A2", "EIF4B", "EIF4H", "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1", "EIF3D" ))}
plot.boxgraph.RNAseq.TCGA <- function(EIF) {
  TCGA.GTEX.RNAseq.sampletype.subset <- TCGA.GTEX.RNAseq.sampletype %>%
    dplyr::select(
      all_of(EIF),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    as_tibble(.) %>%
    melt(.,
      id = c(
        "sample.type",
        "primary.disease",
        "primary.site",
        "study"
      ),
      value.name = "RNAseq"
    ) %>%
    na.omit(.$primary.site) %>%
    # filter(RNAseq != 0) %>%
    mutate_if(is.character, as.factor)

  # boxplot to compare relative abundance of genes across tumors
  .RNAseq.all.gene(TCGA.GTEX.RNAseq.sampletype.subset) %>%
    .RNAseq.grouped.boxplot()

  # boxplot to compare RNA-seq of one gene in tumor vs adjacent normal
  RNAseq.ind.gene.df <- lapply(EIF, .RNAseq.ind.gene,
                               df = TCGA.GTEX.RNAseq.sampletype.subset)
  lapply(RNAseq.ind.gene.df, .RNAseq.boxplot)

  # violin plot to compare  expression in primary, metastatic tumors vs NATs
  .RNAseq.tumortype(TCGA.GTEX.RNAseq.sampletype.subset) %>% .violinplot()
}

#' Compare the RNA ratios between EIF4F genes
#' @description This function generates box plot to compare the RNA ratios between EIF4F genes in tumors tumors vs adjacent normal tissues across TCGA cancer types.
#' and a violin plot to compare RNA ratios in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined.
#' @param x gene name as a string
#' @param y gene name as a string
#' @param z gene name as a string
#' @details  This function first calculate the ratios of mRNA levels between input genes within each TCGA sample,
#' in the data frame \code{TCGA.GTEX.RNAseq.sampletype} prepared from \code{\link{initialize.RNAseq.data}}.
#'
#' With the RNA ratio data \code{RNAratio.data}, it compares RNA ratio across
#' all tumors types with the function \code{\link{.RNAratio.ind}} and plots with \code{\link{.RNAratio.boxplot}}.
#'
#' Finally, it compares RNA ratio in primary, metastatic tumors vs adjacent normal tissues from all TCGA cancer types combined with
#' the function \code{\link{.RNAratio.tumortype}} and plots with \code{\link{.violinplot}}.
#' @return box plots for RNA ratios between \code{x,y,z} in TCGA tumors
#' @export
#' @examples \dontrun{plot.boxgraph.RNAseq.TCGA("EIF4G1", "EIF4A1","EIF4E")}
#' @examples \dontrun{plot.boxgraph.RNAseq.TCGA("EIF4G1", "EIF4E", "EIF4EBP1")}
plot.boxgraph.RNAratio.TCGA <- function(x,y,z) {
  genes.names <- c(x, y, z)
  RNAratio.data <- TCGA.GTEX.RNAseq.sampletype %>%
    dplyr::select(
      all_of(genes.names),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    as_tibble(.) %>%
    dplyr::filter(x != 0 & !is.na(primary.site)) %>%
    # calculate the ratio of mRNA counts
    mutate(
      (!!paste0(x,":", "\n", y)) := (!!as.name(x))-(!!as.name(y)),
      (!!paste0(x,":", "\n", z)) := (!!as.name(x))-(!!as.name(z)),
      (!!paste0(y,":", "\n", z)) := (!!as.name(y))-(!!as.name(z))
    ) %>%
    dplyr::select(
      (!!paste0(x,":", "\n", y)),
      (!!paste0(x,":", "\n", z)),
      (!!paste0(y,":", "\n", z)),
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ) %>%
    mutate_if(is.character, as.factor) %>%
    na.omit(.)

  .RNAratio.ind(RNAratio.data)%>%
    .RNAratio.boxplot(
      df = .,
      dashline = 1,
      ylimit = c(0, 25),
      filename = "RNAratio1.pdf"
    )

  .RNAratio.tumortype(RNAratio.data) %>%
    .violinplot()
}


