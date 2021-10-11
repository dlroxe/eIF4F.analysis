if (!exists("LUAD.Proteomics")) {
  LUAD.Proteomics <- read_excel(
    file.path(data.file.directory, "Protein.xlsx"),
    col_names = FALSE
  ) %>%
    # as.data.frame(.) %>%
    mutate(...1 = make.unique(...1)) %>% # relabel the duplicates
    column_to_rownames(var = "...1") %>%
    t(.) %>%
    as_tibble(.) %>%
    mutate_at(vars(-Type, -Sample), funs(as.numeric)) # exclude two columns convert character to number
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

#### protin-protein correlation LUAD -------------------------------------------
EIF.pro.correlation <- function() {
  LUAD.Pro <- LUAD.Proteomics[LUAD.Proteomics$Type %in% "Tumor", ]

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

## Boxplots for phosphor-proteomics --------------------------------------------
CPTAC.LUAD.Sampletype <- read_excel(
  file.path(
    data.file.directory,
    "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Samples_r1_May2019.xlsx"
  )
)

CPTAC.LUAD.Clinic <- read_excel(
  file.path(
    data.file.directory,
    "S046_BI_CPTAC3_LUAD_Discovery_Cohort_Clinical_Data_r1_May2019.xlsx"
  ),
  sheet = 2
)

CPTAC.LUAD.Clinic.Sampletype <- merge(CPTAC.LUAD.Clinic,
  CPTAC.LUAD.Sampletype,
  by.x = "case_id",
  by.y = "Participant ID (case_id)"
) %>%
  select("tumor_stage_pathological", "Aliquot (Specimen Label)", "Type") %>%
  rename("Sample" = "Aliquot (Specimen Label)") %>%
  mutate(tumor_stage_pathological = case_when(
    Type == "Normal" ~ "Normal",
    tumor_stage_pathological %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
    tumor_stage_pathological %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
    tumor_stage_pathological %in% c("Stage III", "Stage IIIA", "Stage IIIB") ~ "Stage III",
    tumor_stage_pathological %in% c("Stage IV") ~ "Stage IV"
  ))

.get.EIF.CPTAC.LUAD.Proteomics <- function(x) {
  LUAD.Proteomics %>%
    select_if(names(.) %in% c(x, "Sample"))
  # select(all_of(x), "Sample")
}

.get.EIF.CPTAC.LUAD.Phos <- function(x) {
  read_excel(file.path(data.file.directory, "Phos.xlsx"),
    col_names = FALSE
  ) %>%
    filter(...1 %in% c(x, "Sample")) %>%
    # as.data.frame(.) %>%
    mutate(phosname = paste(...1, ...2)) %>%
    column_to_rownames(var = "phosname") %>%
    select(-c(...1, ...2)) %>%
    t(.) %>%
    as_tibble(.) %>%
    mutate_at(vars(-`Sample na`), funs(as.numeric)) %>%
    rename("Sample" = "Sample na")
}

.protein.boxplot <- function(df, x) {
  hline <- summarise(group_by(df, Gene, Type), MD = 2**median(normalize)) %>%
    ungroup(.) %>%
    filter(Gene == x & Type == "Tumor") %>%
    select(., MD) %>%
    as.numeric(.) %>%
    # hline <- dataMedian$MD[dataMedian$Gene == x & dataMedian$Type == "Tumor"] %>%
    round(., digits = 3)
  p2 <- ggplot(
    data = df[df$Gene == x, ],
    aes(
      x = tumor_stage_pathological,
      y = 2**normalize
    )
  ) +
    geom_boxplot(
      data = df[df$Gene == x, ],
      aes(fill = Gene),
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

plot.EIF4.CPTAC.pro.LUAD <- function(EIF_list) {
  EIF.CPTAC.LUAD.Proteomics <- .get.EIF.CPTAC.LUAD.Proteomics(EIF_list)
  EIF.CPTAC.LUAD.Phos <- .get.EIF.CPTAC.LUAD.Phos(EIF_list)
  EIF.LUAD.Phos.Proteomics.Sampletype <- list(
    EIF.CPTAC.LUAD.Phos,
    EIF.CPTAC.LUAD.Proteomics,
    CPTAC.LUAD.Clinic.Sampletype
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
    filter(Type == "Normal") %>%
    group_by(Gene) %>%
    # group_by(Gene) %>%
    summarise(NAT.mean = median(value)) %>%
    left_join(EIF.LUAD.Phos.Proteomics.Sampletype, ., by = "Gene") %>% # right_join is possible with the dev dplyr
    # group_by(Gene) %>%
    mutate(normalize = value - NAT.mean) # %>%

  lapply(sort(unique(EIF.LUAD.Phos.Proteomics.Sampletype.Normalization$Gene)),
    .protein.boxplot,
    df = EIF.LUAD.Phos.Proteomics.Sampletype.Normalization
  )
}


# Run master functions ---------------------------------------------------------
# EIF.pro.correlation()
# plot.EIF4.CPTAC.pro.LUAD(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "AKT1", "MTOR", "EIF4B", "EIF4H"))
