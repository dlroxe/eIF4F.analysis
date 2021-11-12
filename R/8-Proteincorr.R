# prepare phosphoproteomics datasets from CPTAC LUAD
CPTAC_LUAD_Phos <- CPTAC_LUAD_Clinic_Sampletype <- NULL
initialize_phosphoproteomics_data <- function() {
  CPTAC_LUAD_Phos <<- read_excel(file.path(data.file.directory, "Phos.xlsx"),
    col_names = FALSE
  )


  CPTAC_LUAD_Clinic <- read_excel(file.path(
    data.file.directory,
    "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx"
  ),
  sheet = 2
  )

  CPTAC_LUAD_sampletype <- read_excel(
    file.path(
      data.file.directory,
      "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx"
    )
  ) %>%
    as.data.frame() %>%
        dplyr::select("Aliquot (Specimen Label)", "Type", "Participant ID (case_id)") %>%
        dplyr::distinct(.data$`Aliquot (Specimen Label)`, .keep_all = TRUE) # %>%
      # remove_rownames() %>%
      # column_to_rownames(var = "Aliquot (Specimen Label)")

  CPTAC_LUAD_Clinic_Sampletype <<- merge(CPTAC_LUAD_Clinic,
    CPTAC_LUAD_sampletype,
    by.x = "case_id",
    by.y = "Participant ID (case_id)"
  ) %>%
    dplyr::select("tumor_stage_pathological", "Aliquot (Specimen Label)", "Type") %>%
    rename("Sample" = "Aliquot (Specimen Label)") %>%
    mutate(tumor_stage_pathological = case_when(
      Type == "Normal" ~ "Normal",
      tumor_stage_pathological %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
      tumor_stage_pathological %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
      tumor_stage_pathological %in% c("Stage III", "Stage IIIA", "Stage IIIB") ~ "Stage III",
      tumor_stage_pathological %in% c("Stage IV") ~ "Stage IV"
    ))
}

.Scatter.plot <- function(df, x, y, z) {
  p1 <- ggscatter(df,
    x = x,
    y = y, # color = "black",
    add = "reg.line", # conf.int = TRUE,
    add.params = list(color = "black", fill = "lightgray"),
    cor.coef = TRUE,
    cor.method = "pearson",
    color = z,
    xlab = paste(x, "protein expression)"),
    ylab = paste(y, "protein expression)")
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
    path = file.path(output.directory, "LUAD"),
    filename = paste(x, y, "cor.pdf"),
    plot = p1,
    width = 3,
    height = 3,
    useDingbats = FALSE
  )
}

.get.EIF.CPTAC.LUAD.Proteomics <- function(x) {
  CPTAC_LUAD_Proteomics %>%
  #  dplyr::select_if(names(.) %in% c(x, "Sample"))
    select(any_of(x), "Sample")
}

.get.EIF.CPTAC.LUAD.Phos <- function(x) {
  CPTAC_LUAD_Phos %>%
    dplyr::filter(.data$...1 %in% c(x, "Sample")) %>%
    # as.data.frame(.) %>%
    mutate(phosname = paste(.data$...1, .data$...2)) %>%
    tibble::column_to_rownames(var = "phosname") %>%
    dplyr::select(-c(.data$...1, .data$...2)) %>%
    t() %>%
    as_tibble() %>%
    mutate_at(vars(-.data$`Sample na`), funs(as.numeric)) %>%
    rename("Sample" = "Sample na")
}


# Boxplot for phosphor-proteomics across clinic data----------------------------
.protein.boxplot <- function(df, x) {
  hline <- summarise(group_by(df, .data$Gene, .data$Type), MD = 2**median(.data$normalize)) %>%
    ungroup() %>%
    dplyr::filter(.data$Gene == x & .data$Type == "Tumor") %>%
    dplyr::select(.data$MD) %>%
    as.numeric() %>%
    # hline <- dataMedian$MD[dataMedian$Gene == x & dataMedian$Type == "Tumor"] %>%
    round(digits = 3)
  p2 <- ggplot(
    data = df[df$Gene == x, ],
    aes_(
      x = ~ tumor_stage_pathological,
      y = ~ 2**normalize
    )
  ) +
    geom_boxplot(
      data = df[df$Gene == x, ],
      aes_(fill = ~ Gene),
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
    stat_compare_means(
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
    path = file.path(output.directory, "CPTAC"),
    filename = paste0(str_remove(x, ":"), "pro.pdf"),
    plot = p2,
    width = 3,
    height = 3,
    useDingbats = FALSE
  )
}

# protein-protein correlation LUAD ---------------------------------------------
EIF.pro.correlation <- function() {
  LUAD.Pro <- CPTAC_LUAD_Proteomics[CPTAC_LUAD_Proteomics$Type %in% "Tumor", ]

  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "EIF4G1", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "EIF4A1", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "EIF4E", z = "dark blue")

  ## cell division
  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "CKAP2", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "CKAP2", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "CKAP2", z = "dark blue")

  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "CCNA2", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "CCNA2", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "CCNA2", z = "dark blue")

  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "ERCC6L", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "ERCC6L", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "ERCC6L", z = "dark blue")

  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "MCM7", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "MCM7", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "MCM7", z = "dark blue")

  ## translation
  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "RPS2", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "RPS2", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "RPS2", z = "dark blue")

  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "EIF3B", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "EIF3B", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "EIF3B", z = "dark blue")

  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "EIF3G", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "EIF3G", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "EIF3G", z = "dark blue")

  .Scatter.plot(df = LUAD.Pro, x = "EIF4G1", y = "EIF2S3", z = "dark green")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4E", y = "EIF2S3", z = "dark red")
  .Scatter.plot(df = LUAD.Pro, x = "EIF4A1", y = "EIF2S3", z = "dark blue")
}

plot.EIF4.CPTAC.pro.LUAD <- function(EIF_list) {
  EIF.CPTAC.LUAD.Proteomics <- .get.EIF.CPTAC.LUAD.Proteomics(EIF_list)
  EIF.CPTAC.LUAD.Phos <- .get.EIF.CPTAC.LUAD.Phos(EIF_list)
  EIF.LUAD.Phos.Proteomics.Sampletype <- list(
    EIF.CPTAC.LUAD.Phos,
    EIF.CPTAC.LUAD.Proteomics,
    CPTAC_LUAD_Clinic_Sampletype
  ) %>%
    reduce(full_join, by = "Sample") %>%
    pivot_longer(
      cols = -c("Sample", "tumor_stage_pathological", "Type"),
      names_to = "Gene",
      values_to = "value",
      values_drop_na = TRUE
    ) %>%
    mutate_if(is.character, as.factor) %>%
    na.omit(.)

  EIF.LUAD.Phos.Proteomics.Sampletype.Normalization <- EIF.LUAD.Phos.Proteomics.Sampletype %>%
    dplyr::filter(.data$Type == "Normal") %>%
    group_by(.data$Gene) %>%
    # group_by(Gene) %>%
    summarise(NAT.mean = median(.data$value)) %>%
    left_join(EIF.LUAD.Phos.Proteomics.Sampletype, ., by = "Gene") %>% # right_join is possible with the dev dplyr
    # group_by(Gene) %>%
    mutate(normalize = .data$value - .data$NAT.mean) # %>%

  lapply(sort(unique(EIF.LUAD.Phos.Proteomics.Sampletype.Normalization$Gene)),
    .protein.boxplot,
    df = EIF.LUAD.Phos.Proteomics.Sampletype.Normalization
  )
}
