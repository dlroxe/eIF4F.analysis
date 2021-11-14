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
plot_bargraph_CNV_TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))

plot_matrix_CNVcorr_TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))

plot_boxgraph_CNVratio_TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))

# Run master functions in 4-DEG.R ----------------------------------------------
plot_boxgraph_RNAseq_TCGA(c(
  "EIF4G1", "EIF4G2", "EIF4G3", "PABPC1",
  "EIF4A1", "EIF4A2", "EIF4B", "EIF4H",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4EBP1", "EIF3D"
))

plot_boxgraph_RNAratio_TCGA("EIF4G1", "EIF4A1","EIF4E")

# Run master functions in 5-Survival.R -----------------------------------------
plot_KM_RNAseq_TCGA(
  EIF = "EIF4E",
  cutoff = 0.3,
  tumor = "lung adenocarcinoma"
)

plot_KM_RNAseq_TCGA(
  EIF = "EIF4E",
  cutoff = 0.2,
  tumor = "All"
)

plot_CoxPH_RNAseq_TCGA(c(
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", # "RPS6KB1",
  "MYC"
), "All")

plot_CoxPH_RNAseq_TCGA(c(
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", # "RPS6KB1",
  "MYC"
), "lung adenocarcinoma")


# Run master functions in 6-PCA.R ----------------------------------------------
plot_PCA_TCGA_GTEX(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot_PCA_TCGA_GTEX_tumor(
  c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "PABPC1", "MKNK1", "MKNK2"),
  "Lung"
)

plot_PCA_CPTAC_LUAD(c(
  "EIF4E", "EIF4G1", "EIF4A1", "PABPC1",
  "MKNK1", "MKNK2", "EIF4EBP1"
))

# Run master functions in 7-Correlation.R --------------------------------------
lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "PABPC1"),
       plot_scatter_RNApro_CCLE)
lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "PABPC1"),
       plot_scatter_RNApro_LUAD)
plot.Venn.all(x = "All")
plot.Venn.all(x = "Lung")


# Run master functions in 8-Proteincorr.R --------------------------------------
EIF.pro.correlation()
plot_boxgraph_protein_CPTAC(c("EIF4G1", "EIF4A1", "EIF4E"))
plot_boxgraph_protein_CPTAC(c(
  "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
  "AKT1", "MTOR", "EIF4B", "EIF4H",
  "MKNK1", "MKNK2"
))
