# prepare TCGA survival and RNA-seq dataset-------------------------------------
#' Read RNA-seq and survival datasets from TCGA
#'
#' @description This function reads RNA-seq and survival datasets from TCGA
#'
#' TCGA.RNAseq: the RNAseq data from TCGA generated from \code{\link{.get.TCGA.RNAseq}}
#'
#' TCGA.OS: the annotation data from the "Survival_SupplementalTable_S1_20171025_xena_sp"
#' dataset with selection of OS and OS.time columns.
#'
#' TCGA.sampletype: the annotation data from the "TCGA_phenotype_denseDataOnlyDownload.tsv"
#' dataset with selection of sample.type and primary.disease columns.
#'
#' TCGA.RNAseq.OS.sampletype: the merged dataset from TCGA.RNAseq, TCGA.OS and TCGA.sampletype
#'
#' @importFrom purrr reduce
#'
#' @export
#'
#' @examples \dontrun{initialize.survival.data()}
#'
initialize.survival.data <- function() {
  TCGA.RNAseq <- .get.TCGA.RNAseq()
  ## get OS data ##
  TCGA.OS <- data.table::fread(
    file.path(
      data.file.directory,
      "Survival_SupplementalTable_S1_20171025_xena_sp"
    ),
    data.table = FALSE
  ) %>%
    {
      dplyr::distinct(., sample, .keep_all = TRUE) %>%
        # remove_rownames() %>%
        # column_to_rownames(var = 'sample') %>%
        dplyr::select("sample", "OS", "OS.time") %>%
        rename(rn = sample)
    }
  ## get sample type data ##
  TCGA.sampletype <- readr::read_tsv(
    file.path(
      data.file.directory,
      "TCGA_phenotype_denseDataOnlyDownload.tsv"
    ),
    show_col_types = FALSE
  ) %>%
    {
      as_tibble(.) %>%
        dplyr::distinct(., sample, .keep_all = TRUE) %>%
        dplyr::select(
          "sample",
          "sample_type",
          "_primary_disease"
        ) %>%
        rename(
          rn = sample,
          sample.type = sample_type,
          primary.disease = `_primary_disease`
        )
    }
  ## combine OS, sample type and RNAseq data ##
  TCGA.RNAseq.OS.sampletype <<- list(TCGA.RNAseq, TCGA.OS, TCGA.sampletype) %>%
    reduce(full_join, by = "rn") %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "rn") %>%
    dplyr::filter(sample.type != "Solid Tissue Normal")
}

#' Read the RNAseq data from TCGA
#' @description This function reads the RNAseq data from TCGA
#' "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena".
#' @details The function also removes possible duplicated tumor samples and samples with NAs in the dataset.
#'
#' It should not be used directly, only inside \code{\link{initialize.survival.data}} function.
#' @return a data frame that contains the RNAseq data from TCGA
#' @examples \dontrun{.get.TCGA.RNAseq()}
#' @keywords internal
.get.TCGA.RNAseq <- function() {
  TCGA.pancancer <- fread(
    file.path(
      data.file.directory,
      "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
    ),
    data.table = FALSE
  )
  TCGA.RNAseq <- TCGA.pancancer[
    !duplicated(TCGA.pancancer$sample),
    !duplicated(colnames(TCGA.pancancer))
  ] %>%
    as_tibble(.) %>%
    # na.omit(.) %>%
    tibble::remove_rownames(.) %>%
    tibble::column_to_rownames(var = "sample")
  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.RNAseq_transpose <- data.table::transpose(TCGA.RNAseq)
  # get row and colnames in order
  # rownames(TCGA.RNAseq_transpose) <- colnames(TCGA.RNAseq)
  colnames(TCGA.RNAseq_transpose) <- rownames(TCGA.RNAseq)
  # TCGA.RNAseq_transpose$rn <- rownames(TCGA.RNAseq_transpose)
  TCGA.RNAseq_transpose$rn <- colnames(TCGA.RNAseq)
  return(TCGA.RNAseq_transpose)
}


# Survival analysis and plotting -----------------------------------------------
##  KM survival analyses
#' Kaplan Meier survival analyses of gene expression
#' @description This function correlate the gene expression within tumor samples with patient overall survival time from TCGA tumors.
#' @details It should not be used directly, only inside \code{\link{plot.km.RNAseq.TCGA}} function.
#' @param gene gene name, passed \code{EIF} argument from \code{\link{plot.km.RNAseq.TCGA}}
#' @param data \code{df} generated inside \code{\link{plot.km.RNAseq.TCGA}}
#' @param cutoff percentage of gene expression
#' @param tumor all tumor types or specific type
#' @return a KM plot
#' @importFrom reshape2 dcast melt
#' @importFrom survival survfit survdiff
#' @importFrom stats pchisq
#' @examples \dontrun{.KM.curve(gene = EIF, data = df, cutoff = cutoff, tumor = tumor)}
#' @keywords internal
.KM.curve <- function(gene, data, cutoff, tumor) {
  km <- survival::survfit(SurvObj ~ data$Group, data = data, conf.type = "log-log")
  stats <- survival::survdiff(SurvObj ~ data$Group, data = data, rho = 0) # rho = 0 log-rank
  p.val <- (1 - pchisq(stats$chisq, length(stats$n) - 1)) %>%
    signif(., 3)

  KM <- ggplot2::autoplot(
    km,
    censor = FALSE,
    xlab = "Days",
    ylab = "Survival Probability",
    # main = paste0("All TCGA cancer studies (", nrow(df), " cases)"),
    # xlim = c(0, 4100),
    color = strata
  ) +
    {
      if (tumor == "All") {
        labs(title = paste0("All TCGA cancer studies (", nrow(data), " cases)"))
      } else {
        labs(title = paste0(tumor, "studies (", nrow(data), " cases)"))
      }
    } +
    theme_bw() +
    theme(
      plot.title = black_bold_16,
      axis.title = black_bold_16,
      axis.text = black_bold_16,
      # axis.line.x = element_line(color = "black"),
      # axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      strip.text = black_bold_16,
      legend.text = black_bold_16,
      legend.title = black_bold_16,
      legend.position = c(0.9, 0.98),
      legend.justification = c(1, 1)
    ) +
    guides(fill = "none") +
    scale_color_manual(
      values = c("red", "blue"),
      name = paste(gene, "mRNA expression"),
      breaks = c("Bottom %", "Top %"),
      labels = c(
        paste("Bottom ", percent(cutoff), ", n =", round((nrow(data)) * cutoff, digits = 0)),
        paste("Top ", percent(cutoff), ", n =", round((nrow(data)) * cutoff, digits = 0))
      )
    ) +
    annotate(
      "text",
      x        = 10000,
      y        = 0.75,
      label    = paste("log-rank test \n p.val = ", p.val),
      size     = 6.5,
      hjust    = 1,
      fontface = "bold"
    )

  print(KM)
  ggplot2::ggsave(
    path = file.path(output.directory, "KM"),
    filename = paste(gene, " all tumors KM.pdf"),
    plot = KM,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )
}

## Cox regression model and forest plot
#' Univariable Cox-PH analyses of gene expression
#' @description This function generates univariable regression model of
#' the gene expression within tumor samples and patient overall survival time TCGA.
#' @details It should not be used directly, only inside \code{\link{plot.coxph.RNAseq.TCGA}} function.
#' @param gene gene names, passed \code{EIF} argument from \code{\link{plot.coxph.RNAseq.TCGA}}
#' @param data \code{df1} generated inside \code{\link{plot.coxph.RNAseq.TCGA}}
#' @param covariate_names gene names from the input arguement of \code{\link{plot.coxph.RNAseq.TCGA}}
#' @return a table of univariable Cox-PH
#' @importFrom dplyr across arrange bind_rows desc full_join slice vars
#' @importFrom stats as.formula
#' @importFrom survival coxph cox.zph
#' @importFrom survivalAnalysis analyse_multivariate
#' @importFrom purrr map
#' @examples \dontrun{.univariable.analysis(df = df1, covariate_names = EIF)}
#' @keywords internal
.univariable.analysis <- function(df, covariate_names) {
  # Multiple Univariate Analyses
  res.cox <- map(covariate_names, function(gene) {
    survivalAnalysis::analyse_multivariate(df,
      vars(OS.time, OS),
      covariates = list(gene)
    ) %>%
      with(summaryAsFrame)
  }) %>%
    bind_rows()

  # To test for the proportional-hazards (PH) assumption
  test.ph <- map(covariate_names, function(x) {
    coxph(as.formula(paste("Surv(OS.time, OS)~", x)),
      data = df
    ) %>%
      cox.zph() %>%
      print() %>%
      as.data.frame(.) %>%
      slice(1)
  }) %>%
    bind_rows() %>%
    dplyr::select("p") %>%
    rename("pinteraction" = "p") %>%
    tibble::rownames_to_column()

  data1 <- full_join(res.cox, test.ph, by = c("factor.id" = "rowname")) %>%
    # as.data.frame(.) %>%
    mutate(across(7:11, round, 3)) %>%
    mutate(across(4:6, round, 2)) %>%
    mutate(np = nrow(df)) %>%
    mutate(HRCI = paste0(HR, " (", Lower_CI, "-", Upper_CI, ")")) %>%
    mutate(p = case_when(
      p < 0.001 ~ "<0.001",
      # p > 0.05 ~ paste(p,"*"),
      TRUE ~ as.character(p)
    )) %>%
    mutate(pinteraction = case_when(
      pinteraction < 0.001 ~ "<0.001",
      pinteraction > 0.05 ~ paste(pinteraction, "*"),
      TRUE ~ as.character(pinteraction)
    )) %>%
    arrange(desc(HR))

  return(data1)
}


#' Multivariable Cox-PH analyses of gene expression
#' @description This function generates univariable regression model of
#' the gene expression within tumor samples and patient overall survival time TCGA.
#' @details It should not be used directly, only inside \code{\link{plot.coxph.RNAseq.TCGA}} function.
#' @param gene gene names, passed \code{EIF} argument from \code{\link{plot.coxph.RNAseq.TCGA}}
#' @param data \code{df1} generated inside \code{\link{plot.coxph.RNAseq.TCGA}}
#' @param covariate_names gene names from the input argument of \code{\link{plot.coxph.RNAseq.TCGA}}
#' @return a table of multivariable Cox-PH
#' @importFrom dplyr across arrange bind_rows desc full_join slice vars
#' @importFrom stats as.formula
#' @importFrom survival coxph cox.zph
#' @importFrom survivalAnalysis analyse_multivariate
#' @importFrom purrr map
#' @examples \dontrun{.univariable.analysis(df = df1, covariate_names = EIF)}
#' @keywords internal
.multivariable.analysis <- function(df, covariate_names) {
  res.cox <- survivalAnalysis::analyse_multivariate(
    df,
    vars(OS.time, OS),
    covariates = covariate_names
  ) %>%
    with(summaryAsFrame) #  to extract an element from a list

  # To test for the proportional-hazards (PH) assumption
  test.ph <- coxph(as.formula(paste(
    "Surv(OS.time, OS)~",
    paste(covariate_names, collapse = "+")
  )),
  data = df
  ) %>%
    cox.zph() %>%
    print() %>%
    as.data.frame(.) %>% # do not use as_tibble here, cause errors
    dplyr::select("p") %>%
    rename("pinteraction" = "p") %>%
    rownames_to_column() %>%
    dplyr::filter(rowname != "GLOBAL") # remove the global test result for graph

  data1 <- full_join(res.cox, test.ph, by = c("factor.id" = "rowname")) %>%
    mutate(across(7:11, round, 3)) %>%
    mutate(across(4:6, round, 2)) %>%
    mutate(np = nrow(df)) %>%
    mutate(HRCI = paste0(HR, " (", Lower_CI, "-", Upper_CI, ")")) %>%
    mutate(p = case_when(
      p < 0.001 ~ "<0.001",
      # p > 0.05 ~ paste(p,"*"),
      TRUE ~ as.character(p)
    )) %>%
    mutate(pinteraction = case_when(
      pinteraction < 0.001 ~ "<0.001",
      pinteraction > 0.05 ~ paste(pinteraction, "*"),
      TRUE ~ as.character(pinteraction)
    )) %>%
    arrange(desc(HR))

  return(data1)
}


#' Forest plots of COX-PH results
#' @description This function should not be used directly, only inside \code{\link{plot.coxph.RNAseq.TCGA}} function.
#' @param data output dataset generated from \code{\link{.univariable.analysis}} or \code{\link{.multivariable.analysis}} function
#' @param output.file output file name
#' @param plot.title title name of the forest graph
#' @param x.tics tics for x axis
#' @param x.range range for x axis
#' @importFrom forestplot forestplot fpTxtGp fpColors
#' @keywords internal
.forest.graph <- function(data, output.file, plot.title, x.tics, x.range) {
  tabletext1 <- cbind(
    c("Gene", data$factor.id),
    c("No. of\nPatients", data$np),
    c("Hazard Ratio\n(95% CI)", data$HRCI),
    c("P Value", data$p),
    c("P Value for\nInteraction", data$pinteraction)
  )

  # print the forest plot on screen
  p <- forestplot::forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "<---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(box = "black", lines = "black"),
    xticks = x.tics,
    # xlog = 0,
    clip = x.range,
    zero = 1,
    cex = 1.2,
    lineheight = "auto", # height of the graph
    boxsize = 0.2,
    colgap = unit(6, "mm"), # the gap between column
    lwd.ci = 2,
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
    new_page = getOption("forestplot_new_page", TRUE)
  )
  print(p)
  # print the forest plot as a pdf file
  pdf(
    file = output.file,
    width = 14,
    height = 12,
    onefile = F
  )
  p <- forestplot::forestplot(
    labeltext = tabletext1,
    graph.pos = 3, graphwidth = unit(12, "cm"),
    hrzl_lines = list(
      "1" = gpar(lwd = 1, col = "black"),
      "2" = gpar(lwd = 1, col = "black")
    ),
    mean = c(NA, data$HR),
    lower = c(NA, data$Lower_CI),
    upper = c(NA, data$Upper_CI),
    title = plot.title,
    xlab = "<---Good prognosis---    ---Poor prognosis--->",
    txt_gp = fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = fpColors(box = "black", lines = "black"),
    xticks = x.tics,
    clip = x.range,
    zero = 1,
    cex = 1.2,
    lineheight = "auto", # height of the graph
    boxsize = 0.2,
    colgap = unit(6, "mm"), # the gap between column
    lwd.ci = 2,
    ci.vertices = FALSE,
    ci.vertices.height = 0.02,
    new_page = getOption("forestplot_new_page", TRUE)
  )
  print(p)
  dev.off()
}


# master functions to call Survival analysis and plotting ----------------------
#' Survival analyses of TCGA patients with expression of \code{EIF} in their tumors by Kaplan Meier method
#' @description This function generates a Kaplan Meier plot to compare the expression of one gene in TCGA cancer types.
#' @details  This function first selects RNAseq of the query gene,
#' survival data and cancer types from the dataset \code{TCGA.RNAseq.OS.sampletype} prepared from \code{\link{initialize.survival.data}}.
#'
#' With the subset data \code{df}, it compares the survival data from patients with top or bottom percents of gene expression,
#' and plot the results as a KM curve plot with \code{\link{.KM.curve}}
#' @param EIF gene name
#' @param cutoff percentage of gene expression for patient stratification
#' @param tumor all tumor types or specific type
#' @return KM curve plots for TCGA patients with expression of \code{EIF} in their tumors
#' @importFrom survival Surv
#' @importFrom stats quantile
#' @export
#' @examples \dontrun{plot.km.EIF.tumor(EIF = "EIF4E", cutoff = 0.2, tumor = "lung adenocarcinoma")}
#' @examples \dontrun{plot.km.EIF.tumor(EIF = "EIF4G1", cutoff = 0.3, tumor = "All")}
plot.km.RNAseq.TCGA <- function(EIF, cutoff, tumor) {
  df <- TCGA.RNAseq.OS.sampletype %>%
    dplyr::select(
      all_of(EIF),
      "OS",
      "OS.time",
      "sample.type",
      "primary.disease"
    ) %>%
    # drop_na(EIF) %>%
    dplyr::filter(if (tumor != "All") primary.disease == tumor else TRUE) %>%
    drop_na(.) %>%
    # na.omit(.) %>%
    as_tibble(.) %>%
    rename(RNAseq = EIF) %>%
    mutate(Group = case_when(
      RNAseq < quantile(RNAseq, cutoff) ~ "Bottom %",
      RNAseq > quantile(RNAseq, (1 - cutoff)) ~ "Top %"
    )) %>%
    mutate(SurvObj = Surv(OS.time, OS == 1))
  .KM.curve(gene = EIF, data = df, cutoff = cutoff, tumor = tumor)
}

#' Survival analyses of TCGA patients with expression of \code{EIF} in their tumors by Cox-PH method
#' @description This function makes regression models between the gene expression within tumor samples and patient overall survival time.
#' @details  This function first selects RNAseq of the query gene,
#' survival data and cancer types from the dataset \code{TCGA.RNAseq.OS.sampletype} prepared from \code{\link{initialize.survival.data}}.
#'
#' With the subset data \code{df1}, it makes univariable regression models with \code{\link{.univariable.analysis}}
#' and multivariable regression models with \code{\link{.multivariable.analysis}}.
#'
#' It plots the results as a forest graph with \code{\link{.forest.graph}}
#' @param EIF gene name
#' @param tumor all tumor types or specific type
#' @return forest graph showing the relation between survival of TCGA patients and expression of \code{EIF} in their tumors
#' @export
#' @examples \dontrun{plot.coxph.RNAseq.TCGA(c("EIF4E", "EIF4E2", "EIF4E3",
#' "EIF4G1", "EIF4G2", "EIF4G3", "EIF4A1", "EIF4A2", "EIF3D", "EIF3E", "EIF4EBP1",
#' "EIF4EBP2", "MKNK1", "MKNK2", "EIF4B", "EIF4H", "MTOR", "MYC"), "All")}
#' @examples \dontrun{plot.coxph.RNAseq.TCGA(c("EIF4E", "EIF4E2", "EIF4E3",
#' "EIF4G1", "EIF4G2", "EIF4G3", "EIF4A1", "EIF4A2", "EIF3D", "EIF3E", "EIF4EBP1",
#' "EIF4EBP2", "MKNK1", "MKNK2", "EIF4B", "EIF4H", "MTOR", "MYC"), "lung adenocarcinoma")}
plot.coxph.RNAseq.TCGA <- function(EIF, tumor) {
  df1 <- TCGA.RNAseq.OS.sampletype %>%
    dplyr::filter(sample.type != "Solid Tissue Normal") %>%
    dplyr::select(
      all_of(EIF),
      "OS",
      "OS.time",
      "sample.type",
      "primary.disease"
    ) %>%
    # drop_na(EIF) %>%
    dplyr::filter(if (tumor != "All") primary.disease == tumor else TRUE) %>%
    drop_na(.) %>%
    # na.omit(.) %>%
    as.data.frame(.)

  univariable.result <- .univariable.analysis(
    df = df1,
    covariate_names = EIF
  )
  .forest.graph(
    data = univariable.result,
    output.file = if (tumor == "All") {
      file.path(output.directory, "Cox", "EIFUniCox.pdf")
    } else {
      paste0(output.directory, "/Cox/", tumor, "EIFUniCox.pdf")
    },
    plot.title = if (tumor == "All") {
      "Univariable Cox proportional-hazards regression analysis (all tumor types)"
    } else {
      paste("Univariable Cox proportional-hazards regression analysis", tumor)
    },
    x.tics = if (tumor == "All") {
      c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)
    } else {
      c(0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8)
    },
    x.range = if (tumor == "All") {
      c(0.6, 1.8)
    } else {
      c(0.4, 2.8)
    }
  )


  multivariable.result <- .multivariable.analysis(
    df = df1,
    covariate_names = EIF
  )
  .forest.graph(
    data = multivariable.result,
    output.file = if (tumor == "All") {
      file.path(output.directory, "Cox", "EIFmultiCox.pdf")
    } else {
      paste0(output.directory, "/Cox/", tumor, "EIFmultiCox.pdf")
    },
    plot.title = if (tumor == "All") {
      "Multivariable Cox proportional-hazards regression analysis (all tumor types)"
    } else {
      paste("Multivariable Cox proportional-hazards regression analysis", tumor)
    },
    x.tics = if (tumor == "All") {
      c(0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8)
    } else {
      c(0.4, 0.8, 1.2, 1.6, 2.4, 2.8, 3.2)
    },
    x.range = if (tumor == "All") {
      c(0.6, 1.8)
    } else {
      c(0.4, 3.2)
    }
  )
}
