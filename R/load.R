# Directory Preparation -------------------------------------------------------
data_file_directory <- "~/Downloads/Test"
output_directory <- "~/Documents/EIF_output"


#' Set output directories
#'
#' @description This function set up the output directories
#'
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_dir()
#' }
#'
initialize_dir <- function() {
  dir.create(file.path(output_directory, "CNV"), recursive = TRUE)
  dir.create(file.path(output_directory, "DEG"), recursive = TRUE)
  dir.create(file.path(output_directory, "Survival", "KM"), recursive = TRUE)
  dir.create(file.path(output_directory, "Survival", "CoxPH"), recursive = TRUE)
  dir.create(file.path(output_directory, "PCA"))
  dir.create(file.path(output_directory, "PCA", "All"), recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "TCGA"), recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "GTEX"), recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "Lung"), recursive = TRUE)
  dir.create(file.path(output_directory, "CORs"), recursive = TRUE)
  dir.create(file.path(output_directory, "RNApro"), recursive = TRUE)
  dir.create(file.path(output_directory, "Proteomics"), recursive = TRUE)
}


# Format Preparation ----------------------------------------------------------
# due to NSE notes in R CMD check
black_bold_tahoma_7 <- black_bold_12 <- black_bold_12_45 <- black_bold_16 <- NULL
black_bold_16_right <- black_bold_16_45 <- black_bold_16_90 <- black_bold_18 <- NULL
col_vector <- NULL

#' Set format for font and color of plots
#'
#' @description This function set up the font type, size and color for ggplots.
#'
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_format()
#' }
#'
initialize_format <- function() {
  # due to NSE notes in R CMD check
  # black_bold_tahoma_7 <- black_bold_12 <- black_bold_12_45 <- black_bold_16 <- NULL
  # black_bold_16_right <- black_bold_16_45 <- black_bold_16_90 <- black_bold_18 <- NULL

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


#' Generates a number of most distinctive colors for PCA.
#'
#' @description This function should not be used directly,
#' only inside [initialize_format()] function.
#'
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#'
#' @keywords internal
#'
#' @examples \dontrun{
#' color()
#' }
#'
color <- function() {
  qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == "qual", ]
  col_vector <- unlist(mapply(
    RColorBrewer::brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
}
