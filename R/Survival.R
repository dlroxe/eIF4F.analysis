# Survival analyses on EIF4F gene expression in TCGA tumors
# This R script contains four sections.
#
# (1) RNAseq and clinical data preparation
#
# (2) functions to analyze KM and cox analyses as well as plotting
#
# (3) composite functions to execute a pipeline of functions to select related
#  RNAseq and clinical data for survival analysis and plot with supply of
#  EIF4F gene names as values of the arguments.
#
# (4) wrapper function to call all composite functions with inputs
#

## Wrapper function for data initialization of survival and RNA-seq ============

#' @noRd
# due to NSE notes in R CMD check
TCGA_RNAseq_OS_sampletype <- NULL

#' @title Read RNA-seq and survival datasets from TCGA
#'
#' @description
#'
#' A wrapper function to read RNA-seq and survival datasets from TCGA.
#'
#' @details
#'
#' Side effects:
#'
#' (1) `TCGA_RNAseq_OS_sampletype` a merged dataset from
#'  `.TCGA_RNAseq`, `.TCGA_OS` and `.TCGA_sampletype`.
#'
#'  * `.TCGA_RNAseq`: the RNAseq data from TCGA generated from
#'  [.get_TCGA_RNAseq()], which imports the dataset
#'  `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`.
#'
#'  * `.TCGA.OS`: the clinical data from
#'  `Survival_SupplementalTable_S1_20171025_xena_sp`. We select three columns:
#'   `OS` for overall survival status, `OS.time` for overall survival time and
#'   `sample` for sample ID of each patient.
#'
#'  * `.TCGA_sampletype`: the annotation data from the
#'  `TCGA_phenotype_denseDataOnlyDownload.tsv` dataset. We select two columns
#'  `sample.type` that annotates malignant tissues, and `primary.disease`
#'  that annotates cancer types for each sample.
#'
#' Only malignant tissue (Solid normal tissues are excluded) are selected for
#'  survival analysis
#'
#' `TCGA_RNAseq_OS_sampletype` was stored as `TCGA_RNAseq_OS_sampletype.csv` in
#'  `~/Documents/EIF_output/ProcessedData` folder.
#'
#' @family wrapper function for data initialization
#'
#' @importFrom purrr reduce
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_survival_data()
#' }
#'
initialize_survival_data <- function() {
  rlang::env_binding_unlock(parent.env(environment()), nms = NULL)

  if (!file.exists(file.path(
    output_directory,
    "ProcessedData",
    "TCGA_RNAseq_OS_sampletype.csv"
  ))) {
    TCGA_RNAseq <- .get_TCGA_RNAseq()
    ## get OS data ##
    TCGA_OS <- data.table::fread(
      file.path(
        data_file_directory,
        "Survival_SupplementalTable_S1_20171025_xena_sp"
      ),
      data.table = FALSE
    ) %>%
      dplyr::distinct(.data$sample, .keep_all = TRUE) %>%
      # remove_rownames() %>%
      # column_to_rownames(var = 'sample') %>%
      dplyr::select("sample", "OS", "OS.time") %>%
      dplyr::rename(rn = .data$sample)

    ## get sample type data ##
    TCGA_sampletype <- readr::read_tsv(
      file.path(
        data_file_directory,
        "TCGA_phenotype_denseDataOnlyDownload.tsv"
      ),
      show_col_types = FALSE
    ) %>%
      tibble::as_tibble() %>%
      dplyr::distinct(.data$sample, .keep_all = TRUE) %>%
      dplyr::select(
        "sample",
        "sample_type",
        "_primary_disease"
      ) %>%
      dplyr::rename(
        rn = .data$sample,
        sample.type = .data$sample_type,
        primary.disease = .data$`_primary_disease`
      )

    ## combine OS, sample type and RNAseq data ##
    assign("TCGA_RNAseq_OS_sampletype",
      list(TCGA_RNAseq, TCGA_OS, TCGA_sampletype) %>%
        purrr::reduce(full_join, by = "rn") %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = "rn") %>%
        dplyr::filter(.data$sample.type != "Solid Tissue Normal"),
      envir = parent.env(environment())
    )

    data.table::fwrite(TCGA_RNAseq_OS_sampletype, file.path(
      output_directory,
      "ProcessedData",
      "TCGA_RNAseq_OS_sampletype.csv"
    ),row.names = TRUE)

  } else {
    assign("TCGA_RNAseq_OS_sampletype",
      data.table::fread(file.path(
        output_directory,
        "ProcessedData",
        "TCGA_RNAseq_OS_sampletype.csv"),data.table = FALSE
      ) %>%
        tibble::as_tibble() %>%
        #tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = "V1"),
      envir = parent.env(environment())
    )
  }
  rlang::env_binding_lock(parent.env(environment()), nms = NULL)

  return(NULL)
}

#' @title Read the RNAseq data from TCGA
#'
#' @description
#'
#' A helper function reads the RNAseq data from TCGA
#' `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena`.
#'
#' @details
#'
#' This function removes possible duplicated tumor samples and samples with
#'  NAs in the dataset.
#' It should not be used directly, only inside [initialize_survival_data()]
#'  function.
#'
#' @family helper function for data initialization
#'
#' @return
#'
#' a data frame that contains the RNAseq data from TCGA
#'
#' @examples \dontrun{
#' .get_TCGA_RNAseq()
#' }
#'
#' @keywords internal
#'
.get_TCGA_RNAseq <- function() {
  .TCGA_pancancer <- fread(
    file.path(
      data_file_directory,
      "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
    ),
    data.table = FALSE
  )
  .TCGA_RNAseq <- .TCGA_pancancer[
    !duplicated(.TCGA_pancancer$sample),
    !duplicated(colnames(.TCGA_pancancer))
  ] %>%
    tibble::as_tibble() %>%
    # na.omit(.) %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(var = "sample")
  # transpose function from the data.table library keeps numeric values as numeric.
  .TCGA_RNAseq_transpose <- data.table::transpose(.TCGA_RNAseq)
  # get row and colnames in order
  colnames(.TCGA_RNAseq_transpose) <- rownames(.TCGA_RNAseq)
  .TCGA_RNAseq_transpose$rn <- colnames(.TCGA_RNAseq)

  return(.TCGA_RNAseq_transpose)
}


## Survival analysis and plotting ==============================================


#' @title Kaplan Meier survival analyses of gene expression
#'
#' @description
#'
#' A helper function correlates the gene expression within tumor samples with
#'  patient overall survival time from TCGA study groups
#'
#' @details
#'
#' It should not be used directly, only inside [.plot_KM_RNAseq_TCGA()]
#'  function.
#'
#' Side effects:
#' (1) KM curve plots for TCGA patients with gene expression in their
#'  tumors
#'
#' @family helper function for survival analysis
#'
#' @param gene gene name, passed `EIF` argument from [.plot_KM_RNAseq_TCGA()]
#'
#' @param data `df` generated inside [.plot_KM_RNAseq_TCGA()]
#'
#' @param cutoff percentage of gene expression
#'
#' @param tumor all tumor types or specific type
#'
#' @return a KM plot
#'
#' @importFrom reshape2 dcast melt
#'
#' @importFrom scales percent
#'
#' @importFrom survival survfit survdiff
#'
#' @importFrom stats pchisq
#'
#' @examples \dontrun{
#' .KM_curve(gene = EIF, data = df, cutoff = cutoff, tumor = tumor)
#' }
#'
#' @keywords internal
#'
.KM_curve <- function(gene, data, cutoff, tumor) {
  km <- survival::survfit(SurvObj ~ data$Group,
                          data = data,
                          conf.type = "log-log")
  stats <- survival::survdiff(SurvObj ~ data$Group,
                              data = data,
                              rho = 0) # rho = 0 log-rank
  p.val <- (1 - stats::pchisq(stats$chisq, length(stats$n) - 1)) %>%
    signif(3)

  KM <- ggplot2::autoplot(
    km,
    censor = FALSE,
    xlab = "Days",
    ylab = "Survival Probability",
    # main = paste0("All TCGA cancer studies (", nrow(df), " cases)"),
    # xlim = c(0, 4100),
    color = "strata"
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
        paste("Bottom ",
              scales::percent(cutoff),
              ", n =",
              round((nrow(data)) * cutoff, digits = 0)),
        paste("Top ",
              scales::percent(cutoff),
              ", n =",
              round((nrow(data)) * cutoff, digits = 0))
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
    path = file.path(output_directory, "Survival", "KM"),
    filename = paste(gene, cutoff, "cutoff", tumor, "tumors KM.pdf"),
    plot = KM,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )

  return(NULL)
}


#' @title Univariable Cox-PH analyses of gene expression
#'
#' @description
#'
#' A helper function generates univariable regression model of the gene
#'  expression within tumor samples and patient overall survival time TCGA.
#'
#' @details
#'
#' It should not be used directly, only inside [.plot_CoxPH_RNAseq_TCGA()]
#'  function.
#'
#' @family helper function for survival analysis
#'
#' @param gene gene names, passed `gene_list` argument from
#' [.plot_CoxPH_RNAseq_TCGA()]
#'
#' @param data `df1` generated inside [.plot_CoxPH_RNAseq_TCGA()]
#'
#' @param covariate_names gene names from the input argument
#' of [.plot_CoxPH_RNAseq_TCGA()]
#'
#' @return a table of univariable Cox-PH
#'
#' @importFrom dplyr across arrange bind_rows desc full_join slice vars
#'
#' @importFrom stats as.formula
#'
#' @importFrom survival coxph cox.zph
#'
#' @importFrom survivalAnalysis analyse_multivariate
#'
#' @importFrom purrr map
#'
#' @examples \dontrun{
#' .univariable_analysis(df = df1, covariate_names = EIF)
#' }
#'
#' @keywords internal
#'
.univariable_analysis <- function(df, covariate_names) {
  # Multiple Univariate Analyses
  res.cox <- map(covariate_names, function(gene) {
    survivalAnalysis::analyse_multivariate(df,
      dplyr::vars("OS.time", "OS"),
      covariates = list(gene)
    ) %>%
      purrr::pluck("summaryAsFrame") # extracts a named element from a list
  }) %>%
    dplyr::bind_rows()

  # To test for the proportional-hazards (PH) assumption
  test.ph <- map(covariate_names, function(x) {
    survival::coxph(as.formula(paste("Surv(OS.time, OS)~", x)),
      data = df
    ) %>%
      survival::cox.zph() %>%
      print() %>%
      as.data.frame() %>%
      dplyr::slice(1)
  }) %>%
    dplyr::bind_rows() %>%
    dplyr::select("p") %>%
    dplyr::rename("pinteraction" = "p") %>%
    tibble::rownames_to_column()

  return(dplyr::full_join(res.cox,
                          test.ph,
                          by = c("factor.id" = "rowname")) %>%
           # as.data.frame(.) %>%
           dplyr::mutate(dplyr::across(7:11, round, 3)) %>%
           dplyr::mutate(dplyr::across(4:6, round, 2)) %>%
           dplyr::mutate(np = nrow(df)) %>%
           dplyr::mutate(HRCI = paste0(.data$HR,
                                       " (",
                                       .data$Lower_CI,
                                       "-",
                                       .data$Upper_CI, ")")) %>%
           dplyr::mutate(p = dplyr::case_when(
             p < 0.001 ~ "<0.001",
             # p > 0.05 ~ paste(p,"*"),
             TRUE ~ as.character(p)
           )) %>%
           dplyr::mutate(pinteraction = dplyr::case_when(
             pinteraction < 0.001 ~ "<0.001",
             pinteraction > 0.05 ~ paste(pinteraction, "*"),
             TRUE ~ as.character(pinteraction)
           )) %>%
           dplyr::arrange(dplyr::desc(.data$HR)))
}


#' @title Multivariable Cox-PH analyses of gene expression
#'
#' @description
#'
#' This function generates univariable regression model of the gene expression
#'  within tumor samples and patient overall survival time TCGA.
#'
#' @details It should not be used directly, only inside
#' [.plot_CoxPH_RNAseq_TCGA()] function.
#'
#' @family helper function for survival analysis
#'
#' @param gene gene names, passed `gene_list` argument from
#' [.plot_CoxPH_RNAseq_TCGA()]
#'
#' @param data `df1` generated inside [.plot_CoxPH_RNAseq_TCGA()]
#'
#' @param covariate_names gene names from the input argument of
#' [.plot_CoxPH_RNAseq_TCGA()]
#'
#' @return a table of multivariable Cox-PH
#'
#' @importFrom dplyr across arrange bind_rows desc full_join slice vars
#'
#' @importFrom stats as.formula
#'
#' @importFrom survival coxph cox.zph
#'
#' @importFrom survivalAnalysis analyse_multivariate
#'
#' @importFrom purrr map pluck
#'
#' @examples \dontrun{
#' .univariable_analysis(df = df1, covariate_names = EIF)
#' }
#'
#' @keywords internal
#'
.multivariable_analysis <- function(df, covariate_names) {
  res.cox <- survivalAnalysis::analyse_multivariate(
    df,
    dplyr::vars("OS.time", "OS"),
    covariates = covariate_names
  ) %>%
    purrr::pluck("summaryAsFrame") #  to extract an element from a list

  # To test for the proportional-hazards (PH) assumption
  test.ph <- survival::coxph(stats::as.formula(paste(
    "Surv(OS.time, OS)~",
    paste(covariate_names, collapse = "+")
  )),
  data = df
  ) %>%
    survival::cox.zph() %>%
    print() %>%
    as.data.frame() %>% # do not use as_tibble here, cause errors
    dplyr::select("p") %>%
    dplyr::rename("pinteraction" = "p") %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(.data$rowname != "GLOBAL") # remove the global test result for graph

  return(dplyr::full_join(res.cox,
                          test.ph,
                          by = c("factor.id" = "rowname")) %>%
           dplyr::mutate(dplyr::across(7:11, round, 3)) %>%
           dplyr::mutate(dplyr::across(4:6, round, 2)) %>%
           dplyr::mutate(np = nrow(df)) %>%
           dplyr::mutate(HRCI = paste0(.data$HR,
                                       " (",
                                       .data$Lower_CI,
                                       "-",
                                       .data$Upper_CI, ")")) %>%
           dplyr::mutate(p = dplyr::case_when(
             p < 0.001 ~ "<0.001",
             # p > 0.05 ~ paste(p,"*"),
             TRUE ~ as.character(p)
           )) %>%
           dplyr::mutate(pinteraction = dplyr::case_when(
             pinteraction < 0.001 ~ "<0.001",
             pinteraction > 0.05 ~ paste(pinteraction, "*"),
             TRUE ~ as.character(pinteraction)
           )) %>%
           dplyr::arrange(dplyr::desc(.data$HR)))
}


#' @title Forest plots of COX-PH results
#'
#' @description
#'
#' This function should not be used directly,
#' only inside [.plot_CoxPH_RNAseq_TCGA()] function.
#'
#' side effects: forest graph showing the relation between survival of TCGA
#'  patients and EIF expression in their tumors from regression model
#'
#' @family helper function for survival analysis plotting
#'
#' @param data output dataset generated from [.univariable_analysis()] or
#' [.multivariable_analysis()] function
#'
#' @param output.file output file name
#'
#' @param plot.title title name of the forest graph
#'
#' @param x.tics tics for x axis
#'
#' @param x.range range for x axis
#'
#' @importFrom forestplot forestplot fpTxtGp fpColors
#'
#' @keywords internal
#'
.forest_graph <- function(data, output.file, plot.title, x.tics, x.range) {
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
    txt_gp = forestplot::fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = forestplot::fpColors(box = "black", lines = "black"),
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
    txt_gp = forestplot::fpTxtGp(
      label = gpar(cex = 1.2),
      ticks = gpar(cex = 1.2),
      xlab = gpar(cex = 1.2),
      title = gpar(cex = 1.2)
    ),
    col = forestplot::fpColors(box = "black", lines = "black"),
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

  return(NULL)
}


## Composite functions to call Survival analysis and plotting ==================

#' @title Survival analyses of TCGA patients with expression in their tumors by
#'  Kaplan Meier method
#'
#' @description
#'
#' This function generates a Kaplan Meier plot to compare the expression of
#'  one gene in TCGA cancer types.
#'
#' @details This function
#'
#'  * selects RNAseq of the query gene, survival data and cancer types from
#'   the dataset `TCGA_RNAseq_OS_sampletype` prepared from
#'   [initialize_survival_data()].
#'  * compares the survival data from patients with top or bottom percents of
#'   gene expression, and plot the results as a KM curve plot by [.KM_curve()].
#'
#' This function is not accessible to the user and will not show at the users'
#'  workspace. It can only be called by the exported [EIF4F_Survival_analysis()]
#'  function.
#'
#' Side effects:
#'
#' (1) KM curve plots on screen and as pdf files to correlate patient survival
#'  probability with expression of `gene_name` in their tumors
#'
#' @family composite function to call survival analysis and plotting
#'
#' @param gene_name gene name
#'
#' @param cutoff percentage of gene expression for patient stratification
#'
#' @param tumor all tumor types or specific type
#
#' @importFrom survival Surv
#'
#' @importFrom stats quantile
#'
#' @examples \dontrun{
#' plot.km.EIF.tumor(gene_name = "EIF4E", cutoff = 0.2,
#' tumor = "lung adenocarcinoma")
#' }
#'
#' @keywords internal
#'
.plot_KM_RNAseq_TCGA <- function(gene_name, cutoff, tumor) {
  df <- TCGA_RNAseq_OS_sampletype %>%
    dplyr::select(
      dplyr::all_of(gene_name),
      "OS",
      "OS.time",
      "sample.type",
      "primary.disease"
    ) %>%
    # drop_na(EIF) %>%
    dplyr::filter(if (tumor == "All") TRUE
                  else .data$primary.disease == tumor) %>%
    tidyr::drop_na() %>%
    # na.omit(.) %>%
    tibble::as_tibble() %>%
    dplyr::rename(RNAseq = gene_name) %>%
    dplyr::mutate(Group = dplyr::case_when(
      RNAseq < stats::quantile(RNAseq, cutoff) ~ "Bottom %",
      RNAseq > stats::quantile(RNAseq, (1 - cutoff)) ~ "Top %"
    )) %>%
    dplyr::mutate(SurvObj = survival::Surv(.data$OS.time, .data$OS == 1))
  .KM_curve(gene = gene_name, data = df, cutoff = cutoff, tumor = tumor)

  return(NULL)
}


#' @title Survival analyses of TCGA patients with expression in their tumors by
#'  Cox-PH method
#'
#' @description
#'
#' This function makes regression models between the gene expression within
#'  tumor samples and patient overall survival time.
#'
#' @details  This function
#'
#' * selects RNAseq of the query gene, survival data and cancer types from
#'  the dataset `TCGA_RNAseq_OS_sampletype` prepared from
#'  [initialize_survival_data()].
#' * makes univariable regression models with [.univariable_analysis()]
#'  and multivariable regression models with [.multivariable_analysis()].
#' * plots the results as a forest graph with [.forest_graph()]
#'
#' This function is not accessible to the user and will not show at the users'
#'  workspace. It can only be called by the exported [EIF4F_Survival_analysis()]
#'  function.
#'
#' Side effects:
#'
#' (1) forest graphs on screen and as pdf file to show the relation between
#'  survival of TCGA patients and expression of `gene_list` in their tumors
#'  by univariable and multivariable regression models
#'
#' @family composite function to call survival analysis and plotting
#'
#' @param gene_list gene names in a vector of characters
#'
#' @param tumor all tumor types or specific type
#'
#' @importFrom tidyr drop_na
#'
#' @examples \dontrun{
#' .plot_CoxPH_RNAseq_TCGA(c(
#'   "EIF4E", "EIF4E2", "EIF4E3",
#'   "EIF4G1", "EIF4G2", "EIF4G3", "EIF4A1", "EIF4A2", "EIF3D", "EIF3E",
#'   "EIF4EBP1", "EIF4EBP2", "MKNK1", "MKNK2", "EIF4B", "EIF4H", "MTOR", "MYC"
#' ), "All")
#' }
#'
#' @keywords internal
#'
.plot_CoxPH_RNAseq_TCGA <- function(gene_list, tumor) {
  df1 <- TCGA_RNAseq_OS_sampletype %>%
    dplyr::filter(.data$sample.type != "Solid Tissue Normal") %>%
    dplyr::select(
      dplyr::all_of(gene_list),
      "OS",
      "OS.time",
      "sample.type",
      "primary.disease"
    ) %>%
    # drop_na(EIF) %>%
    dplyr::filter(if (tumor != "All") .data$primary.disease == tumor
                  else TRUE) %>%
    tidyr::drop_na() %>%
    # na.omit(.) %>%
    as.data.frame()

  # perform Cox-PH univariable analysis and plot the results
  univariable.result <- .univariable_analysis(
    df = df1,
    covariate_names = gene_list
  )
  .forest_graph(
    data = univariable.result,
    output.file = if (tumor == "All") {
      file.path(output_directory, "Survival", "CoxPH", "All uniCox.pdf")
    } else {
      # paste0(output_directory, "/Survival/", tumor, "EIFUniCox.pdf")
      file.path(output_directory, "Survival", "CoxPH",
                paste(tumor, "uniCox.pdf"))
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

  # perform Cox-PH multivariable analysis and plot the results
  multivariable.result <- .multivariable_analysis(
    df = df1,
    covariate_names = gene_list
  )
  .forest_graph(
    data = multivariable.result,
    output.file = if (tumor == "All") {
      file.path(output_directory, "Survival", "CoxPH", "all multiCox.pdf")
    } else {
      # paste0(output_directory, "/Survival/", tumor, "EIFmultiCox.pdf")
      file.path(output_directory, "Survival", "CoxPH",
                paste(tumor, "multiCox.pdf"))
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

  return(NULL)
}


## Wrapper function to call all composite functions with inputs ================

#' @title Perform all related survival analysis and generate plots
#'
#' @description
#'
#' A wrapper function to call all composite functions for survival analysis with
#'  inputs
#'
#' @details
#'
#' This function runs two internal composite functions together with inputs:
#'
#' * [.plot_KM_RNAseq_TCGA()]
#' * [.plot_CoxPH_RNAseq_TCGA()]
#'
#' Side effects:
#'
#' (1) KM curve plots on screen and as pdf files to correlate patient survival
#'  probability with gene expression in their tumors
#'
#' (2) forest graphs on screen and as pdf file to show the relation between
#'  survival of TCGA patients and gene expression in their tumors
#'  by univariable and multivariable regression models
#'
#' @family wrapper function to call all composite functions with inputs
#'
#' @export
#'
#' @examples \dontrun{
#' EIF4F_Survival_analysis()
#' }
#'
EIF4F_Survival_analysis <- function() {
  lapply(c(
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4A1", "EIF4A2",
    "EIF4E", "EIF4E2", "EIF4E3",
    "EIF3D", "EIF3E",
    "EIF4EBP1", "EIF4EBP2",
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ), .plot_KM_RNAseq_TCGA,
  cutoff = 0.2, tumor = "All"
  )

  lapply(c(
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4A1", "EIF4A2",
    "EIF4E", "EIF4E2", "EIF4E3",
    "EIF3D", "EIF3E",
    "EIF4EBP1", "EIF4EBP2",
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ),
  .plot_KM_RNAseq_TCGA,
  cutoff = 0.3, tumor = "All"
  )

  lapply(c(
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4A1", "EIF4A2",
    "EIF4E", "EIF4E2", "EIF4E3",
    "EIF3D", "EIF3E", "EIF4EBP1", "EIF4EBP2",
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ),
  .plot_KM_RNAseq_TCGA,
  cutoff = 0.2,
  tumor = "lung adenocarcinoma"
  )

  lapply(c(
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4A1", "EIF4A2",
    "EIF4E", "EIF4E2", "EIF4E3",
    "EIF3D", "EIF3E", "EIF4EBP1", "EIF4EBP2",
    "EIF4H", "EIF4B", "MYC",
    "PABPC1", "MKNK1", "MKNK2"
  ),
  .plot_KM_RNAseq_TCGA,
  cutoff = 0.3,
  tumor = "lung adenocarcinoma"
  )

  .plot_CoxPH_RNAseq_TCGA(c(
    "EIF4E", "EIF4E2", "EIF4E3",
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", # "RPS6KB1",
    "MYC"
  ), "All")

  .plot_CoxPH_RNAseq_TCGA(c(
    "EIF4E", "EIF4E2", "EIF4E3",
    "EIF4G1", "EIF4G2", "EIF4G3",
    "EIF4A1", "EIF4A2", "EIF3D",
    "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
    "MKNK1", "MKNK2", "EIF4B", "EIF4H",
    "MTOR", # "RPS6KB1",
    "MYC"
  ), "lung adenocarcinoma")

  return(NULL)
}
