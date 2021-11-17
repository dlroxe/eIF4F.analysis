# prepare RNA proteomics datasets from CCLE and CPTAC LUAD
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

### select EIF expression plotting ---------------------------------------------
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


## RNA protein correlation CCLE-------------------------------------------------
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

## RNA protein correlation LUAD-------------------------------------------------
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


## Heatmap of correlation analysis----------------------------------------------
#' Identify EIF4F correlating genes
#' @description This function calculate the correlation efficiency between each EIF4F gene and the rest of cellular mRNA.
#' and identify the significantly correlating genes.
#' @details It should not be used directly, only inside \code{\link{plot_Corr_RNAseq_TCGA_GTEX}} function.
#' @param df the dataframe \code{.TCGA_GTEX_RNAseq_sampletype_subset} generated inside \code{\link{plot_Corr_RNAseq_TCGA_GTEX}}
#' @param y sample types, either \code{all.tumor.type} or \code{c("Normal Tissue")}
#' generated inside \code{\link{plot_Corr_RNAseq_TCGA_GTEX}}
#' @return a list output with four elements: cor.data for the heatmap function and CORs.counts for
#' bargraph function, c4.posCOR, c4.negCOR for Venn plots
#' @importFrom AnnotationDbi mapIds
#' @importFrom clusterProfiler compareCluster dotplot
#' @importFrom purrr discard
#' @importFrom stats setNames
#' @examples \dontrun{.EIF_correlation(df = .TCGA_GTEX_RNAseq_sampletype_subset,
#' y = all.tumor.type)}
#' @keywords internal
.EIF_correlation <- function(df, y) {
  TCGA.GTEX.tumor <- df[
    df$sample.type %in% y,
  ] %>% na.omit() # select tumor or healthy samples

  Gene.ID <- names(df) %>%
    purrr::discard( ~ .x %in% c(
      "Row.names",
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
  ))

  correlation.coefficient <- function(x, y) {
    result <- cor.test(TCGA.GTEX.tumor[[x]],
      TCGA.GTEX.tumor[[y]],
      method = "pearson"
    )
    res <- data.frame(x,
      y,
      result[c(
        "estimate",
        "p.value",
        "statistic",
        "method"
      )],
      stringsAsFactors = FALSE
    )
  }
  # find all genes positively correlate with EIF4F expression
  # lapply function gives a large list, need to convert it to a dataframe
  EIF.cor.list <- function(EIF) {
    cor.data <- do.call(
      rbind.data.frame,
      lapply(Gene.ID,
        correlation.coefficient,
        y = EIF
      )
    )
    rownames(cor.data) <- cor.data[, 1]
    # cor.data1 <- cor.data[cor.data[, "p.value"] <= 0.05,]
    return(cor.data)
  }


  EIF4E.cor <- EIF.cor.list("EIF4E")
  EIF4G1.cor <- EIF.cor.list("EIF4G1")
  EIF4A1.cor <- EIF.cor.list("EIF4A1")
  EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")


  cor.data <- cbind(
    setNames(data.frame(EIF4E.cor[, c(3, 4)]), c("EIF4E", "EIF4E.p")),
    setNames(data.frame(EIF4G1.cor[, c(3, 4)]), c("EIF4G1", "EIF4G1.p")),
    setNames(data.frame(EIF4A1.cor[, c(3, 4)]), c("EIF4A1", "EIF4A1.p")),
    setNames(data.frame(EIF4EBP1.cor[, c(3, 4)]), c("EIF4EBP1", "EIF4EBP1.p"))
  )


  c4.posCOR <- cbind(
    EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
    EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
    EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
    EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05
  )


  c4.negCOR <- cbind(
    EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
    EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
    EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
    EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05
  )

  #df <- as.data.frame(unclass(summary(c4)))
  #filter_at(vars(starts_with("Sepal"


  count.CORs <- function() {
    c4 <- c4.posCOR
    colnames(c4) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    df <- as.data.frame(summary(c4))
    #df <- as.data.frame(unclass(summary(c4)))
    df1 <- df[df$Freq %like% "TRUE", ]
    df1$Var1 <- NULL
    df1$Var2 <- gsub(" ", "", df1$Var2)
    row.names(df1) <- df1$Var2
    df1$Var2 <- NULL
    df1$Freq <- gsub("TRUE :", "", df1$Freq)
    df1$Freq <- as.numeric(df1$Freq)
    colnames(df1) <- "posCORs"

    c5 <- c4.negCOR
    colnames(c5) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    dt <- as.data.frame(summary(c5))
    dt1 <- dt[dt$Freq %like% "TRUE", ]
    dt1$Var1 <- NULL
    dt1$Var2 <- gsub(" ", "", dt1$Var2)
    row.names(dt1) <- dt1$Var2
    dt1$Var2 <- NULL
    dt1$Freq <- gsub("TRUE :", "", dt1$Freq)
    dt1$Freq <- as.numeric(dt1$Freq)
    colnames(dt1) <- "negCORs"
    df2 <- cbind(df1, dt1)
    return(df2)
  }

  CORs.counts <- count.CORs()
  ## output four variables: cor.data for the heatmap function and CORs.counts for
  ## bargraph function, c4.posCOR, c4.negCOR for Venn plots
  output <- list(cor.data, CORs.counts, c4.posCOR, c4.negCOR)
  return(output)
}

#' Plot Venn diagrams for correlating genes
#' @description This function draw Venn diagrams, using the correlation data
#' generated from \code{\link{.EIF_correlation}}.
#'
#' @param df correlation data
#' @param x input argument of \code{\link{plot_Corr_RNAseq_TCGA_GTEX}}
#' @param z tumor or normal
#' @param CORs posCOR or negCORs
#' @return vennDiagram for posCOR or negCORs
#' @importFrom eulerr euler
#' @importFrom limma vennCounts vennDiagram
#' @examples \dontrun{.Cors_vennplot(df = EIF.cor.tumor[[3]], x = x, z = "tumor", CORs = "posCOR")}
#' @keywords internal
#'
.Cors_vennplot <- function(df, x, z, CORs) {
  b <- limma::vennCounts(df)
  colnames(b) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1", "Counts")
  limma::vennDiagram(b)
  ## eulerr generates area-proportional Euler diagrams that display set
  ## relationships (intersections, unions, and disjoints) with circles or ellipses.
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
    main = paste(x, z, CORs),
    # main = paste(z, Cor),
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
    path = file.path(output.directory, "Heatmap"),
    filename = paste("all", x, z, CORs, "Venn.pdf"),
    plot = p2,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
}

#' Count the numbers of correlating genes
#' @description This function counts the numbers of correlating genes, using the data
#' generated from \code{\link{.EIF_correlation}}.
#'
#' @param df1 \code{EIF.cor.tumor[[2]]}
#' @param df2 \code{EIF.cor.normal[[2]]}
#' @return bargrah for posCOR or negCORs
#' @examples \dontrun{.count_CORs_tumor_normal(df1 = EIF.cor.tumor[[2]],
#' df2 = EIF.cor.normal[[2]])}
#' @keywords internal
#'
.count_CORs_tumor_normal <- function(df1, df2) {
  EIF.cor.counts.tumor <- df1 %>%
    tibble::add_column(label = "tumor") %>%
    tibble::rownames_to_column(var = "gene")

  EIF.cor.counts.normal <- df2 %>%
    tibble::add_column(label = "normal") %>%
    tibble::rownames_to_column(var = "gene")

  EIF.cor <- rbind(EIF.cor.counts.tumor,
    EIF.cor.counts.normal,
    make.row.names = F
  ) %>%
    mutate(label = factor(.data$label,
      levels = c("tumor", "normal"),
      labels = c("Tumors", "Healthy tissues")
    )) %>%
    mutate(gene = factor(.data$gene,
      levels = c("EIF4EBP1", "EIF4A1", "EIF4G1", "EIF4E")
    ))
}

.Cors_bargraph <- function(df, x, CORs, coord_flip.ylim) {
  p1 <- ggplot(
    data = df,
    aes_string(
      x = "gene",
      y = CORs,
      #y = !!sym(CORs), # quote the passed variable CORs
      fill = "label"
    ), color = "label"
  ) +
    geom_bar(
      stat = "identity",
      position = position_dodge()
    ) +
    geom_text(aes_string(label = CORs),
      position = position_dodge(width = 0.9),
      size = 3.5
    ) +
    scale_fill_manual(values = c(
      "#CC79A7", "#0072B2", "#E69F00",
      "#009E73", "#D55E00"
    )) + # for color-blind palettes
    labs(y = paste("number of ", x, CORs)) +
    coord_flip(ylim = c(0, coord_flip.ylim)) +
    guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
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
    path = file.path(output.directory, "Heatmap"),
    filename = paste0("all ", x, CORs, ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
}


.EIF_cor_normal_tumor <- function(df1, df2) {
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
    cor.data$EIF4E.tumor > 0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4E.tumor < -0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor > 0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor < -0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor > 0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor < -0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor > 0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor < -0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4E.normal > 0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4E.normal < -0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4G1.normal > 0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4G1.normal < -0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal > 0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal < -0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal > 0.3 & cor.data$EIF4EBP1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal < -0.3 & cor.data$EIF4EBP1.p.normal <= 0.05,
  ]))
  DF <- DF[, c(1, 3, 5, 7, 9, 11, 13, 15)]
  return(DF)
}

.Cors_heatmap <- function(df, x) {
  # DF_scaled = t(scale(t(DF)))
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options
  ht1 <- ComplexHeatmap::Heatmap(df,
    name = paste("Correlation Coefficient Heatmap", x),
    heatmap_legend_param = list(
      labels_gp = gpar(font = 15),
      legend_width = unit(6, "cm"),
      direction = "horizontal"
    ),
    show_row_names = FALSE,
    show_column_names = FALSE,
    bottom_annotation = HeatmapAnnotation(
      annotation_legend_param = list(direction = "horizontal"),
      type = c(
        "tumor", "tumor", "tumor", "tumor",
        "normal", "normal", "normal", "normal"
      ),
      col = list(type = c(
        "normal" = "royalblue",
        "tumor" = "pink"
      )),
      cn = anno_text(gsub("\\..*", "", colnames(df)),
        location = 0,
        rot = 0,
        just = "center",
        gp = gpar(
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
    row_title_gp = gpar(
      fontsize = 15,
      fontface = "bold"
    ),
    border = TRUE,
    col = circlize::colorRamp2(
      c(-1, 0, 1),
      c("blue", "#EEEEEE", "red")
    )
  )
  ht <- draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )

  pdf(file.path(
    path = file.path(output.directory, "Heatmap"),
    filename = paste(x, " tumors heatmap.pdf")
  ),
  width = 8,
  height = 8,
  useDingbats = FALSE
  )
  ht <- draw(ht1,
    merge_legends = TRUE,
    heatmap_legend_side = "top",
    annotation_legend_side = "top"
  )
  dev.off()
  return(ht1)
}


.get_cluster_pathway_data <- function(df1, df2) {
  cluster.geneID.list <- function(x) {
    c1 <- t(t(row.names(df1[ComplexHeatmap::row_order(df2)[[x]], ])))
    c1 <- as.data.frame(c1)
    c1$V1 <- as.character(c1$V1)
    c1$entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,
      keys = c1$V1,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    # c1 <- c1[!is.na(c1)]
    c1 <- na.omit(c1)
    return(c1$entrez)
  }
  cluster.num <- as.character(c(1:3))
  names(cluster.num) <- paste("cluster", 1:3)
  cluster.data <- lapply(cluster.num, cluster.geneID.list)
  return(cluster.data)
}

.pathway_dotplot <- function(df, p, x) {
  p1 <- clusterProfiler::dotplot(df,
    title = paste("The Most Enriched", p, "Pathways"),
    showCategory = 8,
    font.size = 18,
    includeAll = FALSE
  ) +
    theme_bw() +
    theme(
      plot.title = black_bold_16,
      axis.title = black_bold_16,
      axis.text.x = black_bold_16,
      axis.text.y = black_bold_16,
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.title = black_bold_16,
      legend.text = black_bold_16,
      strip.text = black_bold_16
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output.directory, "Heatmap"),
    filename = paste(x, " tumors", p, ".pdf"),
    plot = p1,
    width = 10,
    height = 8,
    useDingbats = FALSE
  )
}



# master functions to find the overlapping CORs --------------------------------
#' Perform correlation analysis on RNAseq data from all tumors and healthy tissues
#' @description This function find correlating genes of EIF4F from
#' \code{TCGA_GTEX_RNAseq_sampletype} RNAseq data generated by \code{\link{initialize_RNAseq_data}}.
#'
#' It identifies positively or negatively correlating genes from tumor or healthy samples,
#' find the overlap with Vennplot, compare their number with bargraph,
#' visualize the correlation strength with heatmap, and found enriched pathways within correlating genes
#'
#' @param x All or Lung
#' @return Venn digram, bargraph, heatmap, and dotplot
#' @importFrom purrr map pluck discard
#' @export
#' @examples \dontrun{plot_Corr_RNAseq_TCGA_GTEX(x = "All")
#' plot_Corr_RNAseq_TCGA_GTEX(x = "Lung")}
#'
plot_Corr_RNAseq_TCGA_GTEX <- function(x) {
  .TCGA_GTEX_RNAseq_sampletype_subset <- TCGA_GTEX_RNAseq_sampletype %>%
    dplyr::filter(if (x != "All") .data$primary.site == x else TRUE) %>%
    na.omit() %>%
    # mutate_if(is.character, as.factor) %>%
    mutate_at(c(
      "sample.type",
      "primary.disease",
      "primary.site",
      "study"
    ), factor)

  all.tumor.type <- .TCGA_GTEX_RNAseq_sampletype_subset %>%
    dplyr::select(.data$sample.type) %>%
    mutate_if(is.character, as.factor) %>%
    #{
    #  levels(.$sample.type)
    #} %>%
    purrr::map(levels) %>%
    purrr::pluck("sample.type")  %>%
    purrr::discard( ~ .x %in% c("Cell Line", "Normal Tissue", "Solid Tissue Normal"))
   #.[!. %in% c("Cell Line", "Normal Tissue", "Solid Tissue Normal")]


  EIF.cor.tumor <- .EIF_correlation(
    df = .TCGA_GTEX_RNAseq_sampletype_subset,
    y = all.tumor.type
  )
  .Cors_vennplot(df = EIF.cor.tumor[[3]], x = x, z = "tumor", CORs = "posCOR")
  .Cors_vennplot(df = EIF.cor.tumor[[4]], x = x, z = "tumor", CORs = "negCOR")


  EIF.cor.normal <- .EIF_correlation(
    df = .TCGA_GTEX_RNAseq_sampletype_subset,
    y = c("Normal Tissue")
  )
  .Cors_vennplot(df = EIF.cor.normal[[3]], x = x, z = "normal", CORs = "posCOR")
  .Cors_vennplot(df = EIF.cor.normal[[4]], x = x, z = "normal", CORs = "negCOR")


  EIF.cor <- .count_CORs_tumor_normal(
    df1 = EIF.cor.tumor[[2]],
    df2 = EIF.cor.normal[[2]]
  )
  .Cors_bargraph(
    df = EIF.cor, x = x,
    CORs = "posCORs",
    coord_flip.ylim = 14000
  )
  .Cors_bargraph(
    df = EIF.cor, x = x,
    CORs = "negCORs",
    coord_flip.ylim = 14000
  )


  DF <- .EIF_cor_normal_tumor(
    df1 = EIF.cor.tumor[[1]],
    df2 = EIF.cor.normal[[1]]
  )
  ht1 <- .Cors_heatmap(df = DF, x = x)


  cluster.data <- .get_cluster_pathway_data(df1 = DF, df2 = ht1)
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
    fun = "enrichPathway"
  )
  .pathway_dotplot(df = ck.GO, p = "GO", x = x)
  .pathway_dotplot(df = ck.KEGG, p = "KEGG", x = x)
  .pathway_dotplot(df = ck.REACTOME, p = "REACTOME", x = x)
}
