library(eIF4F.analysis)

initialize_data <- function() {
  initialize_cnv_data()
  initialize_RNAseq_data()
  initialize_survival_data()
  initialize_proteomics_data()
  initialize_phosphoproteomics_data()
}
initialize_data()

initialize_format()

# run master functions in 3-CNV.R ----------------------------------------------
EIF4F_CNV_analysis()

# Run master functions in 4-DEG.R ----------------------------------------------
EIF4F_DEG_analysis()

# Run master functions in 5-Survival.R -----------------------------------------
EIF4F_Survival_analysis()

# Run master functions in 6-PCA.R ----------------------------------------------
EIF4F_PCA()

# Run master functions in 7-Correlation.R --------------------------------------
EIF4F_Corrgene_analysis()

# Run master functions in 9-ProRNACorrelation.R --------------------------------------
EIF4F_RNA_pro_correlation()

# Run master functions in 8-Proteincorr.R --------------------------------------
EIF4F_Proteomics_analysis()
