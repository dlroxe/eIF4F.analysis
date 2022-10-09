library(eIF4F.analysis)
library(log4r)
library(R.utils)

require(ReactomePA)
require(reactome.db)

R.utils::setOption("clusterProfiler.download.method", "wget")

my_logfile = "~/eIF4F.analysis.log.txt"

my_logger <- log4r::logger(threshold = "INFO",
                           appenders = list(
                             console_appender(layout = default_log_layout()),
                             file_appender(my_logfile, append = TRUE,
                                           layout = default_log_layout())
                           ))

log4r::info(my_logger, "initializing output dir")
initialize_dir()

log4r::info(my_logger, "initializing formats")
initialize_format()

log4r::info(my_logger, "initializing CNV data")
initialize_cnv_data()
log4r::info(my_logger, "initializing RNAseq data")
initialize_RNAseq_data()
log4r::info(my_logger, "initializing survival data")
initialize_survival_data()
log4r::info(my_logger, "initializing proteomics data")
initialize_proteomics_data()
log4r::info(my_logger, "initializing phosphoproteomics data")
initialize_phosphoproteomics_data()

log4r::info(my_logger, "analyzing CNV data")
EIF4F_CNV_analysis()

log4r::info(my_logger, "analyzing DEG data")
EIF4F_DEG_analysis()

log4r::info(my_logger, "analyzing Survival data")
EIF4F_Survival_analysis()

log4r::info(my_logger, "performing PCA")
EIF4F_PCA()

log4r::info(my_logger, "RNA pro_correlation")
EIF4F_RNA_pro_correlation()

log4r::info(my_logger, "proteomics analysis")
EIF4F_Proteomics_analysis()

log4r::info(my_logger, "correlation analysis")
EIF4F_Corrgene_analysis()

log4r::info(my_logger, "all analyses complete")
