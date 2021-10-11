## Library Preparation ---------------------------------------------------------
requiredPackages <- c(
  "AnnotationDbi", "base", "BiocGenerics", "circlize",
  "clusterProfiler", "ComplexHeatmap", "corrplot",
  "data.table", "data.table", "dplyr", "EnvStats", "eulerr",
  "factoextra", "FactoMineR", "forcats", "forestplot",
  "ggfortify", "ggplot2", "ggpubr", "graphics", "grDevices",
  "grid", "lattice", "limma", "missMDA", "org.Hs.eg.db",
  "purrr", "RColorBrewer", "ReactomePA", "readr", "readxl",
  "reshape2", "scales", "stats", "stats4", "stringr", "survival",
  "survivalAnalysis", "tibble", "tidyr"
)


ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

ipak(requiredPackages)

## Directory Preparation -------------------------------------------------------
data.file.directory <- "~/Downloads/Test"
output.directory <- "~/Documents/EIF_output"


## Format Preparation ----------------------------------------------------------
black_bold_tahoma_7 <- element_text(
  color = "black",
  face = "bold",
  size = 7
)

black_bold_12 <- element_text(
  color = "black",
  face = "bold",
  size = 12
)

black_bold_12_45 <- element_text(
  color = "black",
  face = "bold",
  size = 12,
  angle = 45,
  hjust = 1
)


black_bold_16 <- element_text(
  color = "black",
  face = "bold",
  size = 16
)

black_bold_16_right <- element_text(
  color = "black",
  face = "bold",
  size = 16,
  angle = 90
)

black_bold_16_45 <- element_text(
  color = "black",
  face = "bold",
  size = 16,
  angle = 45,
  hjust = 1
)


black_bold_16_90 <- element_text(
  color = "black",
  face = "bold",
  size = 16,
  angle = 90,
  hjust = 1,
  vjust = 0.5
)

black_bold_18 <- element_text(
  color = "black",
  face = "bold",
  size = 18
)


color <- function() {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
  col_vector <- unlist(mapply(
    brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
}
col_vector <- color()
