  rm(list=ls())
  setwd("/Users/hecate/研一/RNA-seq")
  gene_diff<-read.table("gene_diff.txt",header = T)
### MS4A4A-vs-WT
  MS4A4A_WT<-gene_diff[,c(1:3,7:9)]
### WT_F-vs-WT_M
  WT_F_vs_WT_M<-gene_diff[,c(1:3,16:18)]
### MS4A4A_F-vs-WT_F
  MS4A4A_F_vs_WT_F<-gene_diff[,c(1:3,10:12)]
### MS4A4A_F-vs-MS4A4A_M
  MS4A4A_F_vs_MS4A4A_M<-gene_diff[,c(1:3,4:6)]
### MS4A4A_M -vs-WT_M
  MS4A4A_M_vs_WT_M<-gene_diff[,c(1:3,13:15)]

  a<-na.omit(MS4A4A_F_vs_WT_F[MS4A4A_F_vs_WT_F$diffexp_deseq2_qvalue_WT_F.vs.MS4A4A_F<0.2&MS4A4A_F_vs_WT_F$diffexp_log2fc_WT_F.vs.MS4A4A_F>=0,])
  b<-na.omit(MS4A4A_M_vs_WT_M[MS4A4A_M_vs_WT_M$diffexp_deseq2_qvalue_WT_M.vs.MS4A4A_M<0.2&MS4A4A_M_vs_WT_M$diffexp_log2fc_WT_M.vs.MS4A4A_M>=0,])
  intersect(a$gene_symbol,b$gene_symbol)





  source("test_function.R")
# MS4A4A-vs-WT
  selectQ(data=MS4A4A_WT,taskname = "MS4A4A_WT",Q_value = 0.05)
  
  selectQ(data=MS4A4A_WT,taskname = "MS4A4A_WT",Q_value = 0.1)
  
  selectQ(data=MS4A4A_WT,taskname = "MS4A4A_WT",Q_value = 0.2)
# WT_F-vs-WT_M
  selectQ(data=WT_F_vs_WT_M,taskname = "WT_F_vs_WT_M",Q_value = 0.05)
  
  selectQ(data=WT_F_vs_WT_M,taskname = "WT_F_vs_WT_M",Q_value = 0.1)
  
  selectQ(data=WT_F_vs_WT_M,taskname = "WT_F_vs_WT_M",Q_value = 0.2)

# MS4A4A_F-vs-MS4A4A_M
  selectQ(data=MS4A4A_F_vs_MS4A4A_M,taskname = "MS4A4A_F_vs_MS4A4A_M",Q_value = 0.05)
  
  selectQ(data=MS4A4A_F_vs_MS4A4A_M,taskname = "MS4A4A_F_vs_MS4A4A_M",Q_value = 0.1)
  
  selectQ(data=MS4A4A_F_vs_MS4A4A_M,taskname = "MS4A4A_F_vs_MS4A4A_M",Q_value = 0.2)

# MS4A4A_F-vs-WT_F
  selectQ(data=MS4A4A_F_vs_WT_F,taskname = "MS4A4A_F_vs_WT_F",Q_value = 0.05)
  
  selectQ(data=MS4A4A_F_vs_WT_F,taskname = "MS4A4A_F_vs_WT_F",Q_value = 0.1)
  
  selectQ(data=MS4A4A_F_vs_WT_F,taskname = "MS4A4A_F_vs_WT_F",Q_value = 0.2)

# MS4A4A_M -vs-WT_M
  selectQ(data=MS4A4A_M_vs_WT_M,taskname = "MS4A4A_M_vs_WT_M",Q_value = 0.05)
  
  selectQ(data=MS4A4A_M_vs_WT_M,taskname = "MS4A4A_M_vs_WT_M",Q_value = 0.1)
  
  selectQ(data=MS4A4A_M_vs_WT_M,taskname = "MS4A4A_M_vs_WT_M",Q_value = 0.2)
  
  library(DOSE)
  library(org.Hs.eg.db)
  library(topGO)
  library(clusterProfiler)
  library(pathview)
  library(org.Mm.eg.db)
### hallmark
  toupper()

#### GSEA
##### hallmark
  gsea_hallmark <- read.gmt("h.all.v7.4.entrez.gmt") #读gmt文件

##### reactome
  gsea_reactome <- read.gmt("c2.cp.reactome.v7.4.entrez.gmt") #读gmt文件

  source("test_function.R")

# 0.05
# MS4A4A_WT
  GSEA_plot(taskname="upMS4A4A_WT0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_WT/",reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt"
            ,p_value=0.2)
  GSEA_plot(taskname="downMS4A4A_WT0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_WT0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)
  GSEA_plot(taskname="downMS4A4A_WT0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_WT0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_WT/")
  GO_BP(taskname="downMS4A4A_WT0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_WT/")

# WT_F-vs-WT_M
  GSEA_plot(taskname="upWT_F_vs_WT_M0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downWT_F_vs_WT_M0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upWT_F_vs_WT_M0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downWT_F_vs_WT_M0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upWT_F_vs_WT_M0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/WT_F_vs_WT_M/")
  GO_BP(taskname="downWT_F_vs_WT_M0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/WT_F_vs_WT_M/")

# MS4A4A_F-vs-WT_F
  GSEA_plot(taskname="upMS4A4A_F_vs_WT_F0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_WT_F0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_F_vs_WT_F0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_WT_F0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_F_vs_WT_F0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_WT_F/")
  GO_BP(taskname="downMS4A4A_F_vs_WT_F0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_WT_F/")

# MS4A4A_F-vs-MS4A4A_M
  GSEA_plot(taskname="upMS4A4A_F_vs_MS4A4A_M0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_MS4A4A_M0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_F_vs_MS4A4A_M0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_MS4A4A_M0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_F_vs_MS4A4A_M0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_MS4A4A_M/")
  GO_BP(taskname="downMS4A4A_F_vs_MS4A4A_M0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_F_vs_MS4A4A_M/")

# MS4A4A_M -vs-WT_M
  GSEA_plot(taskname="upMS4A4A_M_vs_WT_M0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_M_vs_WT_M0.05",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_M_vs_WT_M0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_M_vs_WT_M0.05",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_M_vs_WT_M0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_M_vs_WT_M/")
  GO_BP(taskname="downMS4A4A_M_vs_WT_M0.05",dir_path="/Users/hecate/研一/RNA-seq/result/0.05/MS4A4A_M_vs_WT_M/")#无结果,只有CC

# 0.1
# MS4A4A_WT
  GSEA_plot(taskname="upMS4A4A_WT0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_WT/",reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt"
            ,p_value=0.2)
  GSEA_plot(taskname="downMS4A4A_WT0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_WT0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)
  GSEA_plot(taskname="downMS4A4A_WT0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_WT0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_WT/")
  GO_BP(taskname="downMS4A4A_WT0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_WT/")

# WT_F-vs-WT_M
  GSEA_plot(taskname="upWT_F_vs_WT_M0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downWT_F_vs_WT_M0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upWT_F_vs_WT_M0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downWT_F_vs_WT_M0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upWT_F_vs_WT_M0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/WT_F_vs_WT_M/")
  GO_BP(taskname="downWT_F_vs_WT_M0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/WT_F_vs_WT_M/")

# MS4A4A_F-vs-WT_F
  GSEA_plot(taskname="upMS4A4A_F_vs_WT_F0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_WT_F0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_F_vs_WT_F0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_WT_F0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_F_vs_WT_F0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_WT_F/")
  GO_BP(taskname="downMS4A4A_F_vs_WT_F0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_WT_F/")

# MS4A4A_F-vs-MS4A4A_M
  GSEA_plot(taskname="upMS4A4A_F_vs_MS4A4A_M0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_MS4A4A_M0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_F_vs_MS4A4A_M0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_MS4A4A_M0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_F_vs_MS4A4A_M0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_MS4A4A_M/")
  GO_BP(taskname="downMS4A4A_F_vs_MS4A4A_M0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_F_vs_MS4A4A_M/")

# MS4A4A_M -vs-WT_M
  GSEA_plot(taskname="upMS4A4A_M_vs_WT_M0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_M_vs_WT_M0.1",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_M_vs_WT_M0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_M_vs_WT_M0.1",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_M_vs_WT_M0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_M_vs_WT_M/")
  GO_BP(taskname="downMS4A4A_M_vs_WT_M0.1",dir_path="/Users/hecate/研一/RNA-seq/result/0.1/MS4A4A_M_vs_WT_M/")

# 0.2
# MS4A4A_WT
  GSEA_plot(taskname="upMS4A4A_WT0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_WT/",reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt"
            ,p_value=0.2)
  GSEA_plot(taskname="downMS4A4A_WT0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_WT0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研p_value=0.2一/RNA-seq/result/0.2/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)
  GSEA_plot(taskname="downMS4A4A_WT0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_WT/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_WT0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_WT/")
  GO_BP(taskname="downMS4A4A_WT0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_WT/")

# WT_F-vs-WT_M
  GSEA_plot(taskname="upWT_F_vs_WT_M0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downWT_F_vs_WT_M0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upWT_F_vs_WT_M0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downWT_F_vs_WT_M0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/WT_F_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upWT_F_vs_WT_M0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/WT_F_vs_WT_M/")
  GO_BP(taskname="downWT_F_vs_WT_M0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/WT_F_vs_WT_M/")

# MS4A4A_F-vs-WT_F
  GSEA_plot(taskname="upMS4A4A_F_vs_WT_F0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_WT_F0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_F_vs_WT_F0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_WT_F0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_WT_F/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_F_vs_WT_F0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_WT_F/")
  GO_BP(taskname="downMS4A4A_F_vs_WT_F0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_WT_F/")

# MS4A4A_F-vs-MS4A4A_M
  GSEA_plot(taskname="upMS4A4A_F_vs_MS4A4A_M0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_MS4A4A_M0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_F_vs_MS4A4A_M0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_F_vs_MS4A4A_M0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_MS4A4A_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_F_vs_MS4A4A_M0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_MS4A4A_M/")
  GO_BP(taskname="downMS4A4A_F_vs_MS4A4A_M0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_F_vs_MS4A4A_M/")

# MS4A4A_M -vs-WT_M
  GSEA_plot(taskname="upMS4A4A_M_vs_WT_M0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_M_vs_WT_M0.2",type="_GSEA_hallmark",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/mh.all.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GSEA_plot(taskname="upMS4A4A_M_vs_WT_M0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  GSEA_plot(taskname="downMS4A4A_M_vs_WT_M0.2",type="_GSEA_reactome",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_M_vs_WT_M/",
            reference_path="/Users/hecate/研一/RNA-seq/m2.cp.reactome.v0.2.symbols.gmt",p_value=0.2)#无结果
  
  GO_BP(taskname="upMS4A4A_M_vs_WT_M0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_M_vs_WT_M/")
  GO_BP(taskname="downMS4A4A_M_vs_WT_M0.2",dir_path="/Users/hecate/研一/RNA-seq/result/0.2/MS4A4A_M_vs_WT_M/")#无结果,只有CC
