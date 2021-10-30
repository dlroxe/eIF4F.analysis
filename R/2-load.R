## Directory Preparation -------------------------------------------------------
data.file.directory <- "~/Downloads/Test"
output.directory <- "~/Documents/EIF_output"


## Format Preparation ----------------------------------------------------------
initialize.format <- function() {
  black_bold_tahoma_7 <<- element_text(
    color = "black",
    face = "bold",
    size = 7
  )

  black_bold_12 <<- element_text(
    color = "black",
    face = "bold",
    size = 12
  )

  black_bold_12_45 <<- element_text(
    color = "black",
    face = "bold",
    size = 12,
    angle = 45,
    hjust = 1
  )


  black_bold_16 <<- element_text(
    color = "black",
    face = "bold",
    size = 16
  )

  black_bold_16_right <<- element_text(
    color = "black",
    face = "bold",
    size = 16,
    angle = 90
  )

  black_bold_16_45 <<- element_text(
    color = "black",
    face = "bold",
    size = 16,
    angle = 45,
    hjust = 1
  )


  black_bold_16_90 <<- element_text(
    color = "black",
    face = "bold",
    size = 16,
    angle = 90,
    hjust = 1,
    vjust = 0.5
  )

  black_bold_18 <<- element_text(
    color = "black",
    face = "bold",
    size = 18
  )

  col_vector <<- color()
}

color <- function() {
  qual_col_pals <- brewer.pal.info[brewer.pal.info$category == "qual", ]
  col_vector <- unlist(mapply(
    brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
}
