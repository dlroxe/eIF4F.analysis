library(eIF4F.analysis)

initialize.data <- initialize.cnv.data()

# run master functions in 3-CNV.R ----------------------------------------------
plot.bargraph.CNV.TCGA(c(
  "TP53", "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "MYC", "EIF3D", "EIF4EBP1",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "PABPC1", "MKNK1", "MKNK2"
))

plot.matrix.CNVcorr.TCGA(c(
  "TP53", "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "MYC", "EIF3D", "EIF4EBP1",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "PABPC1", "MKNK1", "MKNK2"
))


plot.boxgraph.CNVratio.TCGA(c(
  "TP53", "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "MYC", "EIF3D", "EIF4EBP1",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "PABPC1", "MKNK1", "MKNK2"
))

# Run master functions in 4-DEG.R ----------------------------------------------
plot.boxgraph.RNAseq.TCGA(c(
  "EIF4G1", "EIF4G2", "EIF4G3", "PABPC1",
  "EIF4A1", "EIF4A2", "EIF4B", "EIF4H",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4EBP1", "EIF3D"
))

plot.RNAratio.TCGA(c(
  "EIF4E", "EIF4E2", "EIF4E3", "EIF4EBP1",
  "EIF4G1", "EIF4G2", "EIF4G3", "EIF3D", #
  "EIF4A1", "EIF4A2"
))

plot.cormatrix.RNAseq(c(
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4EBP1", "EIF4EBP2", "MTOR",
  "EIF3C", "EIF3D", "EIF3E", "PABPC1",
  "MKNK1", "MKNK2",
  "TP53", "MYC"
))

# Run master functions in 5-Survival.R -----------------------------------------
lapply(c(
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF3D", "EIF3E",
  "EIF4EBP1", "EIF4EBP2",
  "EIF4H", "EIF4B", "MYC",
  "PABPC1", "MKNK1", "MKNK2"
),
plot.km.EIF.tumor,
cutoff = 0.2, tumor = "All"
)

lapply(c(
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF3D", "EIF3E", "EIF4EBP1", "EIF4EBP2",
  "EIF4H", "EIF4B", "MYC",
  "PABPC1", "MKNK1", "MKNK2"
),
plot.km.EIF.tumor,
cutoff = 0.2,
tumor = "lung adenocarcinoma"
)

plot.km.EIF.tumor(
  EIF = "EIF4E",
  cutoff = 0.3,
  tumor = "lung adenocarcinoma"
)

plot.km.EIF.tumor(
  EIF = "EIF4E",
  cutoff = 0.2,
  tumor = "skin cutaneous melanoma"
)

plot.km.EIF.tumor(
  EIF = "EIF4E",
  cutoff = 0.3,
  tumor = "All"
)

plot.coxph.EIF.tumor(c(
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2", "EIF3D",
  "EIF3E", "EIF4EBP1", "EIF4EBP2", # "PABPC1",
  "MKNK1", "MKNK2", "EIF4B", "EIF4H",
  "MTOR", # "RPS6KB1",
  "MYC"
), "All")

plot.coxph.EIF.tumor(c(
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

plot.PCA.TCGA.GTEX(c(
  "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
  "PABPC1", "MKNK1", "MKNK2",
  "EIF4B", "EIF4H",
  "MYC", "JUN"
))

plot.PCA.TCGA.GTEX(c(
  "EIF4G1", "EIF4G2", "EIF4G3",
  "EIF4A1", "EIF4A2",
  "EIF4E", "EIF4E2", "EIF4E3",
  "EIF4EBP1", "EIF4EBP2", "MTOR",
  "EIF3C", "EIF3D", "EIF3E", "PABPC1",
  "MKNK1", "MKNK2",
  "TP53", "MYC"
))

plot.PCA.TCGA.GTEX.tumor(
  c(
    "EIF4G1", "EIF4A1", "EIF4E", "EIF4EBP1",
    "PABPC1", "MKNK1", "MKNK2"
  ),
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
