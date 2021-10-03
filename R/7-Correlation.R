CCLE.RNAseq <- fread(
  file.path(data.file.directory,
            "CCLE_expression_full.csv"),
  data.table = FALSE)

CCLE.Anno <- fread(
  file.path(data.file.directory,
            "sample_info.csv"),
  data.table = FALSE) %>% select(1, 2)

CCLE.EIF.Proteomics <- fread(
  file.path(data.file.directory,
            "protein_quant_current_normalized.csv"),
  data.table = FALSE)


get.CCLE.EIF.RNAseq <- function (EIF) {
  CCLE.RNAseq <- CCLE.RNAseq %>%
    rename("DepMap_ID" = "V1") %>%
    select("DepMap_ID", starts_with((!!paste(EIF,"(ENSG")) ))
  }

get.CCLE.EIF.Proteomics <- function (EIF) {
  CCLE.EIF.Proteomics <- CCLE.EIF.Proteomics %>%
    filter(Gene_Symbol == EIF)
  if (nrow(CCLE.EIF.Proteomics) > 1) {
    df <- CCLE.EIF.Proteomics %>%
      select(contains("_Peptides")) %>%
      mutate(sum = rowSums(., na.rm = T))
    name <- rownames(df[df$sum == max(df$sum),]) # for the maximum value
    CCLE.EIF.Proteomics <- CCLE.EIF.Proteomics  %>%
      filter(row.names(CCLE.EIF.Proteomics) == name)
  }
  else TRUE
  CCLE.EIF.Proteomics <- CCLE.EIF.Proteomics  %>%
    column_to_rownames(var = 'Gene_Symbol') %>%
    select(contains("TenPx"), -contains("_Peptides")) %>%
    t(.) %>%
    as_tibble(.) %>%
    mutate(Celline = sub("\\_.*", "", row.names(.)),
           Type = sub(".*_ *(.*?) *_.*", "\\1", row.names(.)))%>%
    mutate_if(is.character, as.factor) %>%
    #select(!!paste0(EIF,".pro") := EIF)
    rename(!!paste0(EIF,".pro") := EIF) # rename with dynamic variables
  return(CCLE.EIF.Proteomics)
  }

protein.RNA.scatter.plot <- function(df, x, y){
  p1 <- ggscatter(df,
                  x = paste0(x, ".pro"),
                  y = paste0(x, ".rna"), #color = "Type",
                  add = "reg.line", #conf.int = TRUE,
                  cor.coef = TRUE,
                  cor.method = "pearson",
                  title = x,
                  xlab = "Protein expresion",
                  ylab = "RNA expression")+
    theme_bw() +
    theme(
      plot.title = black_bold_12,
      axis.title.x = black_bold_12,
      axis.title.y = black_bold_12,
      axis.text.x = black_bold_12,
      axis.text.y = black_bold_12,
      panel.grid = element_blank(),
      legend.position = "none",
      strip.text = black_bold_12,
      strip.background = element_rect(fill = "white")
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output.directory, "LUAD"),
    filename = paste(y, x ,"cor",".pdf"),
    plot = p1,
    width = 3,
    height = 3,
    useDingbats = FALSE
  )
}

## Figure 5 ##
######################################
#### RNA protein correlation CCLE ####
######################################
plot.EIF.cor.CCLE <- function (EIF){
  get.CCLE.EIF.RNAseq(EIF) %>%
    full_join(CCLE.Anno,
              by = "DepMap_ID") %>%
    na.omit(.) %>%
    select(-DepMap_ID) %>%
    rename("Celline" = "stripped_cell_line_name",
           !!paste0(EIF,".rna") := contains(EIF)) %>%
    full_join(get.CCLE.EIF.Proteomics(EIF),
              by = "Celline") %>%
    na.omit(.) %>%
    protein.RNA.scatter.plot(df = ., x = EIF, y = "CCLE")
}
#plot.EIF.cor.CCLE("EIF4A1")

#lapply(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1","PABPC1"), plot.EIF.cor.CCLE)


######################################
#### RNA protein correlation LUAD ####
######################################
LUAD.Proteomics <- read_excel(
  file.path(data.file.directory, "Protein.xlsx"),
  col_names = FALSE) %>%
  #as.data.frame(.) %>%
  mutate(...1 = make.unique(...1)) %>% # relabel the duplicates
  column_to_rownames(var = '...1') %>%
  t(.) %>%
  as_tibble(.) %>%
  mutate_at(vars(-Type, -Sample),funs(as.numeric)) # exclude two columns convert character to number


LUAD.RNAseq <- read_excel(
  file.path(data.file.directory, "RNA.xlsx"),
  col_names = FALSE) %>%
  #as_tibble(.) %>%
  mutate(...1 = make.unique(...1)) %>% # relabel the duplicates
  column_to_rownames(var = '...1') %>%
  t(.) %>%
  as_tibble(.) %>%
  mutate_at(vars(-Type, -Sample),funs(as.numeric)) # exclude two columns convert character to number


plot.EIF.cor.LUAD <- function(EIF) {
  EIF.LUAD.Proteomics <- LUAD.Proteomics %>%
    select(all_of(EIF),"Type","Sample") %>%
    filter(Type == "Tumor")

  EIF.LUAD.RNAseq <- LUAD.RNAseq %>%
    select(all_of(EIF),"Type","Sample")%>%
    filter(Type == "Tumor")

  df <-  merge(EIF.LUAD.Proteomics,
               EIF.LUAD.RNAseq,
               by = c("Sample","Type"),
               suffixes = c(".pro",".rna"))

  protein.RNA.scatter.plot(df = df, x = EIF, y = "LUAD")
  }


## Heatmap of correlation analysis ##

# prepare RNA-seq related dataset from TCGA and GTEx----------------------------
# RNA-seq data were import in 4-DEG.R
if (!exists("TCGA.GTEX.RNAseq")){
  get.TCGA.GTEX.RNAseq <- function() {
  TCGA.pancancer <- data.table::fread(
    file.path(data.file.directory,
              "TcgaTargetGtex_RSEM_Hugo_norm_count"),
    data.table = FALSE
  )  %>%
    as_tibble(.) %>%
    distinct(., sample, .keep_all = TRUE) %>%
    na.omit(.) %>%
    remove_rownames(.) %>%
    column_to_rownames(var = 'sample')

  # transpose function from the data.table library keeps numeric values as numeric.
  TCGA.pancancer_transpose <- data.table::transpose(TCGA.pancancer)
  # get row and colnames in order
  rownames(TCGA.pancancer_transpose) <- colnames(TCGA.pancancer)
  colnames(TCGA.pancancer_transpose) <- rownames(TCGA.pancancer)
  return (TCGA.pancancer_transpose)
}
  TCGA.GTEX.RNAseq <- get.TCGA.GTEX.RNAseq()}

if (!exists("TCGA.GTEX.sampletype")){
TCGA.GTEX.sampletype <- readr::read_tsv(
  file.path(data.file.directory,
            "TcgaTargetGTEX_phenotype.txt")) %>% {
              as_tibble(.) %>%
                distinct(., sample, .keep_all = TRUE) %>%
                #na.omit(.) %>%
                remove_rownames() %>%
                column_to_rownames(var = 'sample') %>%
                select("_sample_type",
                       "primary disease or tissue",
                       "_primary_site",
                       "_study") %>%
                rename("sample.type" = "_sample_type",
                       "primary.disease" = "primary disease or tissue",
                       "primary.site" = "_primary_site",
                       "study" = "_study")}
}

if (!exists("TCGA.GTEX.RNAseq.sampletype")){
TCGA.GTEX.RNAseq.sampletype <- merge(TCGA.GTEX.RNAseq,
                                     TCGA.GTEX.sampletype,
                                     by    = "row.names",
                                     all.x = TRUE) %>% {
                                       remove_rownames(.) %>%
                                         column_to_rownames(var = 'Row.names')}
}

EIF.correlation <- function(df, y) {
  TCGA.GTEX.tumor <- df[
    df$sample.type %in% y,
  ] %>% na.omit(.) # select tumor or healthy samples
  Gene.ID <- names(df) %>%
    .[!.%in% c("Row.names",
               "sample.type",
               "primary.disease",
               "primary.site",
               "study")]

  correlation.coefficient <- function(x, y) {
    result <- cor.test(TCGA.GTEX.tumor[[x]],
                       TCGA.GTEX.tumor[[y]],
                       method = "pearson"
    )
    res <- data.frame(x,
                      y,
                      result[c("estimate",
                               "p.value",
                               "statistic",
                               "method")],
                      stringsAsFactors = FALSE
    )
  }
  # find all genes positively correlate with EIF4F expression
  # lapply function gives a large list, need to convert it to a dataframe
  EIF.cor.list <- function(EIF) {
    cor.data <- do.call(rbind.data.frame,
                        lapply(Gene.ID,
                               correlation.coefficient,
                               y = EIF))
    rownames(cor.data) <- cor.data[, 1]
    # cor.data1 <- cor.data[cor.data[, "p.value"] <= 0.05,]
    return(cor.data)
  }


  EIF4E.cor <- EIF.cor.list("EIF4E")
  EIF4G1.cor <- EIF.cor.list("EIF4G1")
  EIF4A1.cor <- EIF.cor.list("EIF4A1")
  EIF4EBP1.cor <- EIF.cor.list("EIF4EBP1")


  cor.data <- cbind(
    setNames(data.frame(EIF4E.cor[, c(3, 4)]), c("EIF4E", "EIF4E.p")),
    setNames(data.frame(EIF4G1.cor[, c(3, 4)]), c("EIF4G1", "EIF4G1.p")),
    setNames(data.frame(EIF4A1.cor[, c(3, 4)]), c("EIF4A1", "EIF4A1.p")),
    setNames(data.frame(EIF4EBP1.cor[, c(3, 4)]), c("EIF4EBP1", "EIF4EBP1.p"))
  )


  c4.posCOR <- cbind(EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
                     EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
                     EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
                     EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05)


  c4.negCOR <- cbind(EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
                     EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
                     EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
                     EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05)


  count.CORs <- function (){
    #c4 <- cbind(
    #  EIF4E.cor$estimate > 0.3 & EIF4E.cor$p.value <= 0.05,
    #  EIF4G1.cor$estimate > 0.3 & EIF4G1.cor$p.value <= 0.05,
    #  EIF4A1.cor$estimate > 0.3 & EIF4A1.cor$p.value <= 0.05,
    #  EIF4EBP1.cor$estimate > 0.3 & EIF4EBP1.cor$p.value <= 0.05
    #)
    c4 <- c4.posCOR
    colnames(c4) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    df <- as.data.frame(summary(c4))
    df1 <- df[df$Freq %like% "TRUE", ]
    df1$Var1 <- NULL
    df1$Var2 <- gsub(" ", "", df1$Var2)
    row.names(df1) <- df1$Var2
    df1$Var2 <- NULL
    df1$Freq <- gsub("TRUE :", "", df1$Freq)
    df1$Freq <- as.numeric(df1$Freq)
    colnames(df1) <- "posCORs"

    #c5 <- cbind(
    #  EIF4E.cor$estimate < -0.3 & EIF4E.cor$p.value <= 0.05,
    #  EIF4G1.cor$estimate < -0.3 & EIF4G1.cor$p.value <= 0.05,
    #  EIF4A1.cor$estimate < -0.3 & EIF4A1.cor$p.value <= 0.05,
    #  EIF4EBP1.cor$estimate < -0.3 & EIF4EBP1.cor$p.value <= 0.05
    #)
    c5 <- c4.negCOR
    colnames(c5) <- c("EIF4E", "EIF4G1", "EIF4A1", "EIF4EBP1")
    dt <- as.data.frame(summary(c5))
    dt1 <- dt[dt$Freq %like% "TRUE", ]
    dt1$Var1 <- NULL
    dt1$Var2 <- gsub(" ", "", dt1$Var2)
    row.names(dt1) <- dt1$Var2
    dt1$Var2 <- NULL
    dt1$Freq <- gsub("TRUE :", "", dt1$Freq)
    dt1$Freq <- as.numeric(dt1$Freq)
    colnames(dt1) <- "negCORs"
    df2 <- cbind(df1, dt1)
    return(df2)}
  CORs.counts <- count.CORs()
  ## output four variables: cor.data for the heatmap function and get.CORs for
  ## bargraph function, c4.posCOR, c4.negCOR for Venn plots
  output <- list(cor.data, CORs.counts, c4.posCOR, c4.negCOR)
  return(output)
}
Venn.plot <- function(df, x, z, CORs) {
  b <- limma::vennCounts(df)
  colnames(b) <- c("EIF4E","EIF4G1","EIF4A1","EIF4EBP1","Counts")
  vennDiagram(b)
  ## eulerr generates area-proportional Euler diagrams that display set
  ## relationships (intersections, unions, and disjoints) with circles or ellipses.
  pos.Venn2 <- eulerr::euler(c(
    EIF4E       = b[9, "Counts"], # EIF4E
    EIF4G1      = b[5, "Counts"], # EIF4G1
    EIF4A1      = b[3, "Counts"], # EIF4A1
    EIF4EBP1    = b[2, "Counts"], # EIF4EBP1
    "EIF4E&EIF4G1"   = b[13, "Counts"],
    "EIF4E&EIF4A1"   = b[11, "Counts"],
    "EIF4E&EIF4EBP1" = b[10, "Counts"],
    "EIF4G1&EIF4A1"   = b[7, "Counts"],
    "EIF4G1&EIF4EBP1"   = b[6, "Counts"],
    "EIF4A1&EIF4EBP1"   = b[4, "Counts"],
    "EIF4E&EIF4G1&EIF4A1" = b[15, "Counts"],
    "EIF4E&EIF4G1&EIF4EBP1" = b[14, "Counts"],
    "EIF4E&EIF4A1&EIF4EBP1" = b[12, "Counts"],
    "EIF4G1&EIF4A1&EIF4EBP1" = b[8, "Counts"],
    "EIF4E&EIF4G1&EIF4A1&EIF4EBP1" = b[16, "Counts"]),
    #shape = "ellipse"
  )
  p2 <- plot(pos.Venn2,
             # key = TRUE,
             main = paste(x, z, CORs),
             #main = paste(z, Cor),
             lwd = 0,
             fill = c("#999999", "#009E73",
                      "#56B4E9", "#E69F00"),
             quantities = list(cex = 1.25),
             labels = list(labels = c("EIF4E","EIF4G1",
                                      "EIF4A1","EIF4EBP1"),
                           cex = 1.25)
  )
  print(p2)
  ggplot2::ggsave(
    path = file.path(output.directory, "Heatmap"),
    filename = paste("all", x ,z, CORs, "Venn.pdf"),
    plot = p2,
    width = 8,
    height = 8,
    useDingbats = FALSE)
}


count.CORs.tumor.normal <- function (df1, df2) {
  EIF.cor.counts.tumor <- df1 %>%
    add_column(label = "tumor") %>%
    tibble::rownames_to_column(., "gene")

  EIF.cor.counts.normal <- df2%>%
    add_column(label = "normal") %>%
    tibble::rownames_to_column(., "gene")

  EIF.cor <- rbind(EIF.cor.counts.tumor,
                   EIF.cor.counts.normal,
                   make.row.names = F) %>%
    mutate(label = factor(label,
                          levels = c("tumor", "normal"),
                          labels = c("Tumors", "Healthy tissues"))) %>%
    mutate(gene = factor(gene,
                         levels = c("EIF4EBP1", "EIF4A1", "EIF4G1", "EIF4E")))
}
CORs.bargraph <- function(df, x, CORs, coord_flip.ylim) {
  p1 <- ggplot(data = df,
               aes(x = gene,
                   y = !!sym(CORs),# quote the passed variable CORs
                   fill = label), color = label) +
    geom_bar(stat = "identity",
             position = position_dodge()) +
    geom_text(aes(label = !!sym(CORs)),
              position = position_dodge(width = 0.9),
              size = 3.5
    ) +
    scale_fill_manual(values = c("#CC79A7", "#0072B2", "#E69F00",
                                 "#009E73", "#D55E00")) + # for color-blind palettes
    labs(y = paste("number of ", x, CORs)) +
    coord_flip(ylim = c(0, coord_flip.ylim)) +
    guides(fill = guide_legend(reverse = TRUE)) + # Flip ordering of legend without altering ordering in plot
    theme_bw() +
    theme(
      plot.title = black_bold_18,
      axis.title.x = black_bold_18,
      axis.title.y = element_blank(),
      axis.text.x = black_bold_18,
      axis.text.y = black_bold_18,
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      panel.grid = element_blank(),
      legend.title = element_blank(),
      legend.text = black_bold_18,
      legend.position = "top",
      legend.justification = "left",
      legend.box = "horizontal",
      strip.text = black_bold_18
    )
  print(p1)
  ggplot2::ggsave(
    path = file.path(output.directory, "Heatmap"),
    filename = paste0("all ", x ,CORs, ".pdf"),
    plot = p1,
    width = 8,
    height = 8,
    useDingbats = FALSE
  )
}


EIF.cor.normal.tumor <- function (df1, df2){
  cor.data <- cbind(
    setNames(
      data.frame(df1[1:8]),
      c("EIF4E.tumor", "EIF4E.p.tumor",
        "EIF4G1.tumor", "EIF4G1.p.tumor",
        "EIF4A1.tumor", "EIF4A1.p.tumor",
        "EIF4EBP1.tumor", "EIF4EBP1.p.tumor")),
    setNames(
      data.frame(df2[1:8]),
      c("EIF4E.normal", "EIF4E.p.normal",
        "EIF4G1.normal", "EIF4G1.p.normal",
        "EIF4A1.normal", "EIF4A1.p.normal",
        "EIF4EBP1.normal", "EIF4EBP1.p.normal"))
  )
  DF <- as.matrix(na.omit(cor.data[
    cor.data$EIF4E.tumor > 0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4E.tumor < -0.3 & cor.data$EIF4E.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor > 0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4G1.tumor < -0.3 & cor.data$EIF4G1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor > 0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4A1.tumor < -0.3 & cor.data$EIF4A1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor > 0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4EBP1.tumor < -0.3 & cor.data$EIF4EBP1.p.tumor <= 0.05 |
      cor.data$EIF4E.normal > 0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4E.normal < -0.3 & cor.data$EIF4E.p.normal <= 0.05 |
      cor.data$EIF4G1.normal > 0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4G1.normal < -0.3 & cor.data$EIF4G1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal > 0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4A1.normal < -0.3 & cor.data$EIF4A1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal > 0.3 & cor.data$EIF4EBP1.p.normal <= 0.05 |
      cor.data$EIF4EBP1.normal < -0.3 & cor.data$EIF4EBP1.p.normal <= 0.05,
  ]))
  DF <- DF[, c(1, 3, 5, 7, 9, 11, 13, 15)]
  return(DF)
}
Heatmap.pathway <- function (df, x){
  # DF_scaled = t(scale(t(DF)))
  ## Creating heatmap with three clusters (See the ComplexHeatmap documentation for more options
  ht1 <- Heatmap(df,
                 name = paste("Correlation Coefficient Heatmap", x),
                 heatmap_legend_param = list(
                   labels_gp = gpar(font = 15),
                   legend_width = unit(6, "cm"),
                   direction = "horizontal"
                 ),
                 show_row_names = FALSE,
                 show_column_names = FALSE,
                 bottom_annotation = HeatmapAnnotation(
                   annotation_legend_param = list(direction = "horizontal"),
                   type = c("tumor", "tumor", "tumor", "tumor",
                            "normal", "normal", "normal", "normal"),
                   col = list(type = c(
                     "normal" = "royalblue",
                     "tumor" = "pink"
                   )),
                   cn = anno_text(gsub("\\..*", "", colnames(df)),
                                  location = 0,
                                  rot = 0,
                                  just = "center",
                                  gp = gpar(fontsize = 15,
                                            fontface = "bold")
                   )
                 ),
                 # cluster_rows  = as.dendrogram(hclust(dist(DF))),
                 # row_split     = 3,
                 row_km = 3,
                 # row_km_repeats = 100,
                 row_title = "cluster_%s",
                 row_title_gp = gpar(
                   fontsize = 15,
                   fontface = "bold"
                 ),
                 border = TRUE,
                 col = circlize::colorRamp2(
                   c(-1, 0, 1),
                   c("blue", "#EEEEEE", "red")
                 )
  )
  ht <- draw(ht1,
             merge_legends = TRUE,
             heatmap_legend_side = "top",
             annotation_legend_side = "top"
  )

  pdf(file.path(
    path = file.path(output.directory, "Heatmap"),
    filename = paste(x ," tumors heatmap.pdf")
  ),
  width = 8,
  height = 8,
  useDingbats = FALSE
  )
  ht <- draw(ht1,
             merge_legends = TRUE,
             heatmap_legend_side = "top",
             annotation_legend_side = "top"
  )
  dev.off()
  return(ht1)
  }


get.cluster.pathway.data <- function(df1, df2){
  cluster.geneID.list <- function(x) {
    c1 <- t(t(row.names(df1[row_order(df2)[[x]], ])))
    c1 <- as.data.frame(c1)
    c1$V1 <- as.character(c1$V1)
    c1$entrez <- mapIds(org.Hs.eg.db,
                        keys = c1$V1,
                        column = "ENTREZID",
                        keytype = "SYMBOL",
                        multiVals = "first"
    )
    # c1 <- c1[!is.na(c1)]
    c1 <- na.omit(c1)
    return(c1$entrez)
  }
  cluster.num <- as.character(c(1:3))
  names(cluster.num) <- paste("cluster", 1:3)
  cluster.data <- lapply(cluster.num, cluster.geneID.list)
  return(cluster.data)
}
pathway.dotplot <- function(df, p, x) {
  p1 <- dotplot(df,
                title = paste("The Most Enriched", p, "Pathways"),
                showCategory = 8,
                font.size = 18,
                includeAll = FALSE
    ) +
      theme_bw() +
      theme(
        plot.title = black_bold_16,
        axis.title = black_bold_16,
        axis.text.x = black_bold_16,
        axis.text.y = black_bold_16,
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"),
        panel.grid = element_blank(),
        legend.title = black_bold_16,
        legend.text = black_bold_16,
        strip.text = black_bold_16
      )
    print(p1)
    ggplot2::ggsave(
      path = file.path(output.directory, "Heatmap"),
      filename = paste(x, " tumors", p, ".pdf"),
      plot = p1,
      width = 10,
      height = 8,
      useDingbats = FALSE
    )}


### find posCOR and negCOR in the overlapping CORs from all cancer cases
plot.Venn.all <- function(x) {
  TCGA.GTEX.RNAseq.sampletype <- TCGA.GTEX.RNAseq.sampletype %>%
    filter(if (x != "All") primary.site == x else TRUE) %>%
    na.omit(.) %>%
    #mutate_if(is.character, as.factor) %>%
    mutate_at(c("sample.type",
                "primary.disease",
                "primary.site",
                "study"), factor)

  all.tumor.type <- TCGA.GTEX.RNAseq.sampletype %>%
    select(sample.type) %>%
    mutate_if(is.character, as.factor) %>%
    {levels(.$sample.type)}%>%
    .[!.%in% c("Cell Line","Normal Tissue","Solid Tissue Normal")]


  EIF.cor.tumor <- EIF.correlation(df = TCGA.GTEX.RNAseq.sampletype,
                                    y = all.tumor.type)
  Venn.plot(df = EIF.cor.tumor[[3]], x = x, z = "tumor", CORs = "posCOR")
  Venn.plot(df = EIF.cor.tumor[[4]], x = x, z = "tumor", CORs = "negCOR")


  EIF.cor.normal <- EIF.correlation(df = TCGA.GTEX.RNAseq.sampletype,
                                     y = c("Normal Tissue"))
  Venn.plot(df = EIF.cor.normal[[3]], x = x, z = "normal", CORs = "posCOR")
  Venn.plot(df = EIF.cor.normal[[4]], x = x, z = "normal", CORs = "negCOR")


  EIF.cor <- count.CORs.tumor.normal(df1 = EIF.cor.tumor[[2]],
                                     df2 = EIF.cor.normal[[2]])
  CORs.bargraph(df = EIF.cor, x = x,
                CORs = "posCORs",
                coord_flip.ylim = 14000)
  CORs.bargraph(df = EIF.cor, x = x,
                CORs = "negCORs",
                coord_flip.ylim = 14000)


  DF <- EIF.cor.normal.tumor (df1 = EIF.cor.tumor[[1]],
                              df2 = EIF.cor.normal[[1]])
  ht1 <- Heatmap.pathway(df = DF, x = x)


  cluster.data <- get.cluster.pathway.data(df1 = DF, df2 = ht1)
  ck.GO <- compareCluster(
    geneCluster = cluster.data,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db"
  )
  ck.KEGG <- compareCluster(
    geneCluster = cluster.data,
    fun = "enrichKEGG"
  )
  ck.REACTOME <- compareCluster(
    geneCluster = cluster.data,
    fun = "enrichPathway"
  )
  pathway.dotplot(df = ck.GO, p ="GO", x = x)
  pathway.dotplot(df = ck.KEGG, p ="KEGG", x = x)
  pathway.dotplot(df = ck.REACTOME, p ="REACTOME", x = x)
}


# Run master functions ---------------------------------------------------------
#lapply(c("EIF4G1","EIF4A1","EIF4E","EIF4EBP1","PABPC1"), plot.EIF.cor.LUAD)

#plot.Venn.all(x = "All")
#plot.Venn.all(x = "Lung")

