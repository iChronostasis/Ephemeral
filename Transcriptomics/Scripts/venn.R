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


#---------------------------------------------------------------------------------------------------------------
#从CRAN安装ggVennDiagram包；
  install.packages("ggVennDiagram")
  install.packages("ggsci")

#载入所需的R包；
  library(ggplot2)
  library(ggsci)
  library(sf)
  library(ggVennDiagram)

#自定义颜色；
  color1 <- alpha("#f8766d",0.9)
  color2 <- alpha("#FF99CC",0.7)
  color3 <- alpha("#c77cff",0.5)
  color4 <- alpha("#99CC00",0.5)

#绘制常见的4组数据venn图；
  ggVennDiagram(x[1:4], label_alpha=0) +
    scale_fill_gradient(low="white",high =color4 ,guide="none")
p <- ggVennDiagram(venn_list,size=1,,lty="longdash") +
  scale_fill_gradient(low="white",high = color4 ,guide="none") +
  labs(title = "Overlap of highly variable features",
       subtitle = "After the standard pre-processing workflow") +
  theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
                                                      size = 20, margin = margin(t = 1, b = 12)),
                            plot.subtitle = element_text(hjust = 0,vjust = .5,size=15),
                            plot.caption = element_text(face = 'bold',size = 12),
                            axis.text.x = element_blank(),
                            axis.text.y = element_blank(),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank())
  ggsave("venn.pdf", plot = p, width = 15, height = 10)  

#--------------------------------------------------------------
### Visualize the venn in another type[4 dataset]
venn_list <- readRDS(paste0(out_data_dir,"venn_list_QC.rds"))
seurat_list <- readRDS(paste0(out_data_dir,"seurat_list_QC.rds"))
install.packages("UpSetR")
library(UpSetR)
pdf(paste0(fig_dir,"Visualization/UpSet_QC_4.pdf"), width=10, height=10)
upset(fromList(venn_list), order.by = "freq")
# upset(fromExpression(expressionInput), order.by = "freq")
dev.off()

