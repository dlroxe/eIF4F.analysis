# Co-expression among EIF4F subunits and differential expression
# This R script contains four sections.
#
# (1) LUAD phospho-proteomics data preparation
#
# (2) selection of RNA and protein expression data and plotting
#
# (3) composite functions to execute a pipeline of functions to select related
#  expression with supply of EIF4F gene names for coexpression and differential
#  expression analyses.
#
# (4) wrapper function to call all composite functions with inputs
#

## Wrapper function for data initialization of phosphoproteomics datasets ======

#' @noRd
## due to NSE notes in R CMD check
CPTAC_LUAD_Phos <- CPTAC_LUAD_Clinic_Sampletype <- NULL


#' @title Read all phosphoproteomics related datasets from CPTAC LUAD
#'
#' @description A wrapper function reads all phosphoproteomics related datasets
#'  from CPTAC LUAD.
#'
#' @details Side effects:
#'
#' (1) `CPTAC_LUAD_Phos`: the phosphoproteomics data of CPTAC LUAD from the
#'  download data `Phos.xlsx`
#'
#' (2) `CPTAC_LUAD_Clinic_Sampletype`: a merged dataset from two data frames.
#'  It provides the annotation data about sample types (tumor or healthy tissues)
#'  and tumors stages of each sample
#'   * `.CPTAC_LUAD_Clinic`, the CPTAC clinical data from the download dataset
#'   `S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx`
#'   * `.CPTAC_LUAD_sampletype`, the CPTAC annotation data from the download
#'   dataset `S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx`
#'
#' @family wrapper function for data initialization
#'
#' @importFrom dplyr case_when
#'
#' @importFrom readxl read_excel
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_phosphoproteomics_data()
#' }
#'
initialize_phosphoproteomics_data <- function() {
  CPTAC_LUAD_Phos <<- readxl::read_excel(file.path(data_file_directory,
                                                   "Phos.xlsx"),
    col_names = FALSE
  )

  .CPTAC_LUAD_Clinic <- readxl::read_excel(file.path(
    data_file_directory,
    "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx"
  ),
  sheet = 2
  )

  .CPTAC_LUAD_sampletype <- readxl::read_excel(
    file.path(
      data_file_directory,
      "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx"
    )
  ) %>%
    as.data.frame() %>%
    dplyr::select("Aliquot (Specimen Label)", "Type",
                  "Participant ID (case_id)") %>%
    dplyr::distinct(.data$`Aliquot (Specimen Label)`, .keep_all = TRUE) # %>%
  # remove_rownames() %>%
  # column_to_rownames(var = "Aliquot (Specimen Label)")

  CPTAC_LUAD_Clinic_Sampletype <<- merge(.CPTAC_LUAD_Clinic,
    .CPTAC_LUAD_sampletype,
    by.x = "case_id",
    by.y = "Participant ID (case_id)"
  ) %>%
    dplyr::select("tumor_stage_pathological", "Aliquot (Specimen Label)",
                  "Type") %>%
    dplyr::rename("Sample" = "Aliquot (Specimen Label)") %>%
    dplyr::mutate(tumor_stage_pathological = dplyr::case_when(
      Type == "Normal" ~ "Normal",
      tumor_stage_pathological %in% c("Stage I", "Stage IA", "Stage IB") ~
      "Stage I",
      tumor_stage_pathological %in% c("Stage II", "Stage IIA", "Stage IIB") ~
      "Stage II",
      tumor_stage_pathological %in% c("Stage III", "Stage IIIA", "Stage IIIB") ~
      "Stage III",
      tumor_stage_pathological %in% c("Stage IV") ~ "Stage IV"
    ))
}


## Select RNA and protein expression data and plotting ===================

#' @title Scatter plots of protein coexpression
#'
#' @description A helper function for protein coexpression analysis
#'
#' This function should not be used directly, only inside
#'  [.plot_scatterplot_protein_CPTAC()] function.
#'
#' Side effects:
#' (1) scatter plot to show correlation between two protein expressions
#'
#' @param df a subset LUAD dataset generated from
#'  [.plot_scatterplot_protein_CPTAC()]
#'
#' @param protein01 protein name
#'
#' @param protein02 protein name
#'
#' @param color color scheme for the scatter plot
#'
#' @importFrom ggpubr ggscatter
#'
#' @family helper function for protein coexpression analysis
#'
#' @keywords internal
#'
.protein_scatterplot <- function(df, protein01, protein02, color) {
  p1 <- ggpubr::ggscatter(df,
    x = protein01,
    y = protein02, # color = "black",
    add = "reg.line", # conf.int = TRUE,
    add.params = list(color = "black", fill = "lightgray"),
    cor.coef = TRUE,
    cor.method = "pearson",
    color = color,
    xlab = paste(protein01, "protein expression)"),
    ylab = paste(protein02, "protein expression)")
  ) +
    # scale_y_continuous(breaks= scales::pretty_breaks())+
    # scale_y_continuous(
    #  breaks = get_breaks(by = 1, from = -1),
    #  limits = c(-1, 2)) + # for 3G
    # Add correlation coefficient
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = black_bold_12,
      axis.title.y = black_bold_12,
      axis.text.x = black_bold_12,
      axis.text.y = black_bold_12,
      # axis.line.x      = element_line(color = "black"),
      # axis.line.y      = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = black_bold_12,
      strip.background = element_rect(fill = "white")
    ) #+
  print(p1)

  ggplot2::ggsave(
    path = file.path(output_directory, "Proteomics"),
    filename = paste(protein01, protein02, "cor.pdf"),
    plot = p1,
    #width = 3,
    #height = 3,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )

  return(NULL)
}

#' @title Select subset of CPTAC proteomics data
#'
#' @description
#'
#' This function selected the CPTAC LUAD proteomics data from the input proteins
#'  `protein_list`.
#'
#' @details
#'
#' The function should not be used directly, only inside
#' [.plot_boxgraph_protein_CPTAC()] function.
#'
#' @param protein_list protein name, passed `protein_list` argument from
#'  [.plot_boxgraph_protein_CPTAC()]
#'
#' @return a data frame of CPTAC LUAD  data from the input `protein_list` genes
#'
#' @family helper function for protein coexpression analysis
#'
#' @examples \dontrun{
#' .get_CPTAC_LUAD_Proteomics_subset(protein_list)
#' }
#'
#' @keywords internal
#'
.get_CPTAC_LUAD_Proteomics_subset <- function(protein_list) {
  return(CPTAC_LUAD_Proteomics %>%
    #  dplyr::select_if(names(.) %in% c(x, "Sample"))
    dplyr::select(dplyr::any_of(protein_list), "Sample"))
}

#' @title Select subset of CPTAC phosproteomics data
#'
#' @description
#'
#' This function selected the CPTAC LUAD phosphoproteomics data from the
#'  input `protein_list`.
#'
#' @details
#'
#' The function should not be used directly, only inside
#'  [.plot_boxgraph_protein_CPTAC()] function.
#'
#' @param protein_list protein name, passed `protein_list` argument from
#'  [.plot_boxgraph_protein_CPTAC()]
#'
#' @return a data frame of CPTAC LUAD  data from the input `protein_list` genes
#'
#' @family helper function for protein coexpression analysis
#'
#' @importFrom dplyr funs
#'
#' @examples \dontrun{
#' .get_CCLE_RNAseq_subset()
#' }
#'
#' @keywords internal
#'
.get_CPTAC_LUAD_Phosproteomics_subset <- function(protein_list) {
  return(CPTAC_LUAD_Phos %>%
    dplyr::filter(.data$...1 %in% c(protein_list, "Sample")) %>%
    # as.data.frame(.) %>%
    dplyr::mutate(phosname = paste(.data$...1, .data$...2)) %>%
    tibble::column_to_rownames(var = "phosname") %>%
    dplyr::select(-c(.data$...1, .data$...2)) %>%
    t() %>%
    tibble::as_tibble() %>%
    dplyr::mutate_at(dplyr::vars(-.data$`Sample na`), dplyr::funs(as.numeric)) %>%
    dplyr::rename("Sample" = "Sample na"))
}

#' @title Boxplots of phosphor-proteomics data across clinic stages
#'
#' @description A helper function draws boxplots for protein differential
#'  expression across clinical stages
#'
#' @details
#' This function should not be used directly,
#'  only inside [.plot_boxgraph_protein_CPTAC()].
#'
#' Side effects:
#' (1) boxplots to show the different expression of protein across
#'  different clinical tumor stages
#'
#' @param df a data frame of the combined proteomics and phosphorylation data
#'  generated inside [.plot_boxgraph_protein_CPTAC()]
#'
#' @param protein_name protein name
#'
#' @family helper function for protein coexpression analysis
#'
#' @importFrom dplyr ungroup
#'
#' @importFrom ggpubr stat_compare_means
#'
#' @importFrom stats median
#'
#' @importFrom stringr str_remove
#'
#' @keywords internal
#'
.protein_boxplot <- function(df, protein_name) {
  hline <- dplyr::summarise(dplyr::group_by(df, .data$Gene, .data$Type),
    MD = 2**stats::median(.data$normalize)
  ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$Gene == protein_name & .data$Type == "Tumor") %>%
    dplyr::select(.data$MD) %>%
    as.numeric() %>%
    round(digits = 3)
  p2 <- ggplot(
    data = df[df$Gene == protein_name, ],
    aes_(
      x = ~tumor_stage_pathological,
      y = ~ 2**normalize
    )
  ) +
    geom_boxplot(
      data = df[df$Gene == protein_name, ],
      aes_(fill = ~Gene),
      # alpha         = 0,
      # size     = .75,
      # width    = 1,
      outlier.shape = NA,
      position = position_dodge(width = .9)
    ) +
    geom_hline(
      yintercept = hline,
      colour = "dark red",
      linetype = "dashed"
    ) +
    annotate("text",
      label = hline,
      x = "Stage IV",
      y = hline,
      vjust = -1,
      size = 4,
      colour = "dark red"
    ) +
    stat_n_text(
      size = 4,
      fontface = "bold", y.pos = 0,
      angle = 0,
      hjust = 0.5
    ) +
    labs(x = NULL, y = "Normalized peptide ratio") +
    # coord_cartesian(ylim = c(0, 3)) +
    # scale_y_continuous(breaks = seq(0, 3, by = 1)) +
    # coord_cartesian(ylim = c(-0.5, 7.5)) +
    # scale_y_continuous(breaks = seq(0, 7.5, by = 1)) +
    # coord_cartesian(ylim = c(-1, 15)) +
    # scale_y_continuous(breaks = seq(-1, 15, by = 2)) +
    scale_x_discrete(limits = c(
      "Normal",
      "Stage I",
      "Stage II",
      "Stage III",
      "Stage IV"
    )) +
    scale_fill_discrete(drop = F) +
    # ggplot2::facet_grid(. ~ Gene) +
    facet_wrap(~Gene) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.x = black_bold_12,
      axis.title.y = black_bold_12,
      axis.text.x = black_bold_12_45,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = black_bold_12,
      strip.background = element_rect(fill = "white")
    ) +
    ggpubr::stat_compare_means(
      comparisons = list(
        c("Normal", "Stage I"),
        c("Normal", "Stage II"),
        c("Normal", "Stage III")
      ),
      method = "t.test",
      label = "p.signif",
      # label.y = c(2.4, 2.7, 3),
      # label.y = c(5, 6, 7),
      # label.y = c(12.4, 13.6, 14.8),
      size = 4
    )
  print(p2)
  ggplot2::ggsave(
    path = file.path(output_directory, "Proteomics"),
    filename = paste0(stringr::str_remove(protein_name, ":"), "pro.pdf"),
    plot = p2,
    #width = 3,
    #height = 3,
    width = 6,
    height = 6,
    useDingbats = FALSE
  )

  return(NULL)
}


## Composite function for coexpression and differential expression analysis ====

#' @title protein-protein coexpression in CPTAC LUAD
#'
#' @description A composite function for coexpression and differential
#'  expression analysis
#'
#' @details This function should not be used directly, only inside
#'  [EIF4F_Proteomics_analysis()] function.
#'
#' Side effects:
#'
#' (1) scatter plots to show correlation between two protein
#'  expressions
#'
#' @family composite function to call protein coexpression and differential
#'  expression analysis
#'
#' @importFrom ggpubr ggscatter
#'
#' @keywords internal
#'
.plot_scatterplot_protein_CPTAC <- function() {
  LUAD.Pro <- CPTAC_LUAD_Proteomics[CPTAC_LUAD_Proteomics$Type %in% "Tumor", ]

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "EIF4G1",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "EIF4A1",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "EIF4E",
                       color = "dark blue")

  ## cell division
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "CKAP2",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "CKAP2",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "CKAP2",
                       color = "dark blue")

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "CCNA2",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "CCNA2",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "CCNA2",
                       color = "dark blue")

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "ERCC6L",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "ERCC6L",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "ERCC6L",
                       color = "dark blue")

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "MCM7",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "MCM7",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "MCM7",
                       color = "dark blue")

  ## translation
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "RPS2",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "RPS2",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "RPS2",
                       color = "dark blue")

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "EIF3B",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "EIF3B",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "EIF3B",
                       color = "dark blue")

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "EIF3G",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "EIF3G",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "EIF3G",
                       color = "dark blue")

  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4G1", protein02 = "EIF2S3",
                       color = "dark green")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4E", protein02 = "EIF2S3",
                       color = "dark red")
  .protein_scatterplot(df = LUAD.Pro, protein01 = "EIF4A1", protein02 = "EIF2S3",
                       color = "dark blue")

  return(NULL)
}

#' @title Comparison of protein and phosphorylation levels among different
#'  stages of LUAD
#'
#' @description A composite function analyzes differential expression of eIF4F protein and
#'  phosphorylation in CPTAC LUAD data and plots results as boxplots
#'
#' @details This function
#'
#' * merges the LUAD proteomics values from EIF4F genes
#'  in the data frame prepared from [.get_CPTAC_LUAD_Proteomics_subset()] and
#'  phosphoproteomics data of the same protein in the data frame prepared
#'  from [.get_CPTAC_LUAD_Phosproteomics_subset()].
#'
#' * uses the combined data to compare the abundance across different
#'  tumor stages, and plot the results with the function [.protein_boxplot()]
#'
#' This function should not be used directly, only inside
#'  [EIF4F_Proteomics_analysis()] function.
#'
#' Side effects:
#' (1) boxplot to show the different expression of (phospho)protein
#'  across different clinical tumor stages
#'
#' @param protein_list protein names in a vector of string
#'
#' @family composite function to call protein coexpression and differential
#'  expression analysis
#'
#' @importFrom dplyr group_by right_join summarise
#'
#' @importFrom tidyr pivot_longer
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' .plot_boxgraph_protein_CPTAC(c(
#'   "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
#'   "AKT1", "MTOR", "EIF4B", "EIF4H", "MKNK1", "MKNK2"
#' ))
#' }
#'
.plot_boxgraph_protein_CPTAC <- function(protein_list) {
  .CPTAC_LUAD_Proteomics_subset <- .get_CPTAC_LUAD_Proteomics_subset(protein_list)
  .CPTAC_LUAD_Phos_subset <- .get_CPTAC_LUAD_Phosproteomics_subset(protein_list)
  EIF.LUAD.Phos.Proteomics.Sampletype <- list(
    .CPTAC_LUAD_Proteomics_subset,
    .CPTAC_LUAD_Phos_subset,
    CPTAC_LUAD_Clinic_Sampletype
  ) %>%
    purrr::reduce(full_join, by = "Sample") %>%
    tidyr::pivot_longer(
      cols = -c("Sample", "tumor_stage_pathological", "Type"),
      names_to = "Gene",
      values_to = "value",
      values_drop_na = TRUE
    ) %>%
    dplyr::mutate_if(is.character, as.factor) %>%
    na.omit()

  Normalization <- EIF.LUAD.Phos.Proteomics.Sampletype %>%
    dplyr::filter(.data$Type == "Normal") %>%
    dplyr::group_by(.data$Gene) %>%
    # group_by(Gene) %>%
    dplyr::summarise(NAT.mean = stats::median(.data$value)) %>%
    # right_join is possible with the dev dplyr
    dplyr::right_join(EIF.LUAD.Phos.Proteomics.Sampletype, by = "Gene") %>%
    # group_by(Gene) %>%
    dplyr::mutate(normalize = .data$value - .data$NAT.mean) # %>%

  lapply(sort(unique(Normalization$Gene)),
    .protein_boxplot,
    df = Normalization
  )

  return(NULL)
}


## Wrapper function to call all composite functions with inputs ================

#' @title Analyze co-expression among EIF4F subunits and differential expression
#'
#' @description A wrapper function to call two composite functions for protein
#'  coexpression and differential expression.
#'
#' @details
#'
#' This function run the composite functions [.plot_scatterplot_protein_CPTAC()]
#' and [.plot_boxgraph_protein_CPTAC()] with EIF4F gene name as inputs.
#'
#' Side effects:
#'
#' (1) scatter plot to show correlation between two protein expressions in CPTAC
#'  LUAD samples
#'
#' (2) boxplot to show the different expression of protein across
#'  different clinical tumor stages
#'
#' @family wrapper function to call all composite functions with inputs
#'
#' @export
#'
#' @examples \dontrun{
#' EIF4F_Proteomics_analysis()
#' }
#'
EIF4F_Proteomics_analysis <- function() {
  .plot_scatterplot_protein_CPTAC()
  .plot_boxgraph_protein_CPTAC(c(
    "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "AKT1", "MTOR", "EIF4B", "EIF4H",
    "MKNK1", "MKNK2"
  ))

  return(NULL)
}
