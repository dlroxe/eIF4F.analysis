library(eIF4F.analysis)

initialize.data <- function() {
  initialize.cnv.data()
  initialize.RNAseq.data()
  initialize.survival.data()
  initialize.proteomics.data()
  initialize.RNApro.data()
  initialize.phosphoproteomics.data()
}
initialize.data()

initialize.format()

# run master functions in 3-CNV.R ----------------------------------------------
plot.bargraph.CNV.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))

plot.matrix.CNVcorr.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))

plot.boxgraph.CNVratio.TCGA(c("EIF4A1", "EIF4E", "EIF4EBP1", "EIF4G1"))

# Run master functions in 4-DEG.R ----------------------------------------------
plot.boxgraph.RNAseq.TCGA(c(
  "EIF4G1", "EIF4G2", "EIF4G3", "PABPC1",
  "EIF4A1", "EIF4A2", "EIF4B", "EIF4H",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4EBP1", "EIF3D"
))

plot.boxgraph.RNAratio.TCGA("EIF4G1", "EIF4A1","EIF4E")

# Run master functions in 5-Survival.R -----------------------------------------
plot.km.RNAseq.TCGA(
  EIF = "EIF4E",
  cutoff = 0.3,
  tumor = "lung adenocarcinoma"
)

plot.km.RNAseq.TCGA(
  EIF = "EIF4E",
  cutoff = 0.2,
  tumor = "All"
)

plot.coxph.RNAseq.TCGA(c(
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", # "RPS6KB1",
  "MYC"
), "All")

plot.coxph.RNAseq.TCGA(c(
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", # "RPS6KB1",
  "MYC"
), "lung adenocarcinoma")


# Run master functions in 6-PCA.R ----------------------------------------------
plot.PCA.TCGA.GTEX(c(
  "EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.PCA.TCGA.GTEX.tumor(
  c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "PABPC1", "MKNK1", "MKNK2"),
  "Lung"
)

plot.PCA.CPTAC.LUAD(c(
  "EIF4E", "EIF4G1", "EIF4A1", "PABPC1",
  "MKNK1", "MKNK2", "EIF4EBP1"
))

# Run master functions in 7-Correlation.R --------------------------------------
lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "PABPC1"), plot.EIF.cor.CCLE)
lapply(c("EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1", "PABPC1"), plot.EIF.cor.LUAD)
plot.Venn.all(x = "All")
plot.Venn.all(x = "Lung")


# Run master functions in 8-Proteincorr.R --------------------------------------
EIF.pro.correlation()
plot.EIF4.CPTAC.pro.LUAD(c(
  "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
  "AKT1", "MTOR", "EIF4B", "EIF4H",
  "MKNK1", "MKNK2"
))
