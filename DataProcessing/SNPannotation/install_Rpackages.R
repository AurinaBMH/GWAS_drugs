if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install(c("biomaRt"))
