  rm(list=ls())
  setwd("/Users/hecate/研一/RNA-seq/DESeq2/venn")
  library(VennDiagram)
### 比较一下
  source("~/test_function.R")
### MS4A4A-vs-WT
  vennplot(Gene1=MS4A4A_WT,Gene2=res,taskname = "MS4A4A-vs-WT_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_WT,Gene2=res,taskname = "MS4A4A-vs-WT_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_WT,Gene2=res,taskname = "MS4A4A-vs-WT_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### WT_F-vs-WT_M
  vennplot(Gene1=WT_F_vs_WT_M,Gene2=res,taskname = "WT_F-vs-WT_M_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=WT_F_vs_WT_M,Gene2=res,taskname = "WT_F-vs-WT_M_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=WT_F_vs_WT_M,Gene2=res,taskname = "WT_F-vs-WT_M_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### MS4A4A_F-vs-WT_F
  vennplot(Gene1=MS4A4A_F_vs_WT_F,Gene2=res,taskname = "MS4A4A_F_vs_WT_F_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_F_vs_WT_F,Gene2=res,taskname = "MS4A4A_F_vs_WT_F_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_F_vs_WT_F,Gene2=res,taskname = "MS4A4A_F_vs_WT_F_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### MS4A4A_F-vs-MS4A4A_M
  vennplot(Gene1=MS4A4A_F_vs_MS4A4A_M,Gene2=res,taskname = "MS4A4A_F_vs_MS4A4A_M_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_F_vs_MS4A4A_M,Gene2=res,taskname = "MS4A4A_F_vs_MS4A4A_M_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_F_vs_MS4A4A_M,Gene2=res,taskname = "MS4A4A_F_vs_MS4A4A_M_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### MS4A4A_M -vs-WT_M
  vennplot(Gene1=MS4A4A_M_vs_WT_M,Gene2=res,taskname = "MS4A4A_M_vs_WT_M_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_M_vs_WT_M,Gene2=res,taskname = "MS4A4A_M_vs_WT_M_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  vennplot(Gene1=MS4A4A_M_vs_WT_M,Gene2=res,taskname = "MS4A4A_M_vs_WT_M_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
  


### overlap
### MS4A4A-vs-WT
overlap(Gene1=MS4A4A_WT,Gene2=res,taskname = "MS4A4A-vs-WT_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_WT,Gene2=res,taskname = "MS4A4A-vs-WT_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_WT,Gene2=res,taskname = "MS4A4A-vs-WT_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### WT_F-vs-WT_M
overlap(Gene1=WT_F_vs_WT_M,Gene2=res,taskname = "WT_F-vs-WT_M_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=WT_F_vs_WT_M,Gene2=res,taskname = "WT_F-vs-WT_M_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=WT_F_vs_WT_M,Gene2=res,taskname = "WT_F-vs-WT_M_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### MS4A4A_F-vs-WT_F
overlap(Gene1=MS4A4A_F_vs_WT_F,Gene2=res,taskname = "MS4A4A_F_vs_WT_F_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_F_vs_WT_F,Gene2=res,taskname = "MS4A4A_F_vs_WT_F_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_F_vs_WT_F,Gene2=res,taskname = "MS4A4A_F_vs_WT_F_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### MS4A4A_F-vs-MS4A4A_M
overlap(Gene1=MS4A4A_F_vs_MS4A4A_M,Gene2=res,taskname = "MS4A4A_F_vs_MS4A4A_M_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_F_vs_MS4A4A_M,Gene2=res,taskname = "MS4A4A_F_vs_MS4A4A_M_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_F_vs_MS4A4A_M,Gene2=res,taskname = "MS4A4A_F_vs_MS4A4A_M_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
### MS4A4A_M -vs-WT_M
overlap(Gene1=MS4A4A_M_vs_WT_M,Gene2=res,taskname = "MS4A4A_M_vs_WT_M_0.05",Q_value = 0.05,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_M_vs_WT_M,Gene2=res,taskname = "MS4A4A_M_vs_WT_M_0.1",Q_value = 0.1,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")
overlap(Gene1=MS4A4A_M_vs_WT_M,Gene2=res,taskname = "MS4A4A_M_vs_WT_M_0.2",Q_value = 0.2,dir_path = "/Users/hecate/研一/RNA-seq/DESeq2/venn")







