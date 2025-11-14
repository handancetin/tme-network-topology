setwd("~/Documents/CSBL/github-projects/tme-cfm-gems-mesb/r")

# Load packages
library("TCGAbiolinks")
library("survival")

install.packages(survival)



TCGAbiolinks:::getProjectSummary("TCGA-COAD")

query_TCGA = GDCquery(
  project = "TCGA-LIHC",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")