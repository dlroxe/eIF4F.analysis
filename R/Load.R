#' ### Directory Preparation
## Directory Preparation =======================================================
#' @noRd
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
  dir.create(file.path(output_directory, "CNV"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "DEG"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Survival", "KM"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Survival", "CoxPH"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "All"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "TCGA"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "GTEX"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "PCA", "Lung"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "CORs"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "RNApro"),
             showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_directory, "Proteomics"),
             showWarnings = FALSE, recursive = TRUE)
}


#' ### Data preparation
## Data preparation ============================================================
#' Read datasets from the load data files
#'
#' @description This function read datasets from the load data files.
#'
#'
#' @export
#'
#' @examples \dontrun{
#' initialize_data()
#' }
#'
initialize_data <- function() {
  initialize_cnv_data()
  initialize_RNAseq_data()
  initialize_survival_data()
  initialize_proteomics_data()
  initialize_phosphoproteomics_data()
}


#' ### Format Preparation
## Format Preparation ==========================================================
#' @noRd
## due to NSE notes in R CMD check
black_bold_tahoma_7 <- black_bold_12 <- black_bold_12_45 <- NULL
black_bold_16 <- black_bold_16_right <- NULL
black_bold_16_45 <- black_bold_16_90 <- black_bold_18 <- NULL
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
  qual_col_pals <- RColorBrewer::brewer.pal.info[
    RColorBrewer::brewer.pal.info$category == "qual", ]
  col_vector <- unlist(mapply(
    RColorBrewer::brewer.pal,
    qual_col_pals$maxcolors,
    rownames(qual_col_pals)
  ))
}


