  rm(list=ls())  
  setwd("")
## 导入所有reads的count值，header = T保证后续计算不会将列名算入而造成错误。如果file和脚本不在同一个路径，记得添加file的绝对路径。 
  count <- read.table(file="/Users/hecate/研一/Term_2/mouse placental/2022.03.31/out_feature_again.txt",header = T,sep="\t")[,-12:-16] #从第二列开始，将每个样本的count循环转化成FPKM 
  length <- read.table(file="/Users/hecate/研一/Term_2/mouse placental/2022.03.31/out_feature_again.txt",header = T,sep="\t")$Length 
## 导入所有reads的count值，header = T保证后续计算不会将列名算入而造成错误。如果file和脚本不在同一个路径，记得添加file的绝对路径。
  n <- as.data.frame(length)
# 1. 转换为FPKM值
## 从第二列开始，将每个样本的count循环转化成FPKM
  merge <- cbind(count,n)
## 从第二列开始，将每个样本的count循环转化成FPKM
  i <- 2
  gene_num <- nrow(count)
  sample_num <- 10
  mapped_reads <- sum(merge[1:gene_num,2],na.rm = T)
  
  repeat{
    mapped_reads <- sum(merge[1:gene_num,i],na.rm = T)#计算每个样本的mapped reads数
    FPKM <- merge[1:gene_num,i]/(10^-9*mapped_reads*merge[1:gene_num,dim(merge)[2]])#计算FPKM值
    FPKM <- pmax(FPKM,0)#去掉矫正带来的负值
    count = data.frame(count[1:gene_num,],FPKM)#添加到count表格后面i
    i <- i + 1
    
    if(i > sample_num+1){
      break
    }
  }
  

## 生成表格列名称
  count_colname <- read.table("/Users/hecate/研一/Term_2/mouse placental/2022.03.31/out_feature_again.txt",header = F,nrow = 1,as.is=TRUE)[,2:11]
  FPKM_colname <- paste(count_colname[1,],"_FPKM",sep="")
  colname <- c("Genename",count_colname,FPKM_colname)
  names(count) <- colname

## 生成表格
  write.csv(count,"ALL_FPKM.csv",row.names = FALSE)#row.names是为了去除第一列自动生成的行名，同理col.names=FALSE可以去除第一行自动生成的列名
  head(read.csv("ALL_FPKM.csv"))

  dev.new()
## 样品间相关性分析
  library(pheatmap)
  #data=read.csv("generead_counts.csv",header = T,row.names = 1)
  data=read.table("/Users/hecate/研一/Term_2/mouse placental/2022.03.31/count.txt",header = T,sep="\t")[,-1]
  matrix=cor(data, use = "complete.obs")
  write.csv(matrix,"gene_coefficient_matrix.csv")
  pdf(file="gene_cor.pdf")
  pheatmap(matrix,cluster_rows=F,cluster_cols = F ,display_numbers = T,fontsize = 5)
  dev.off()


## 2. 基因表达分布(箱线图)
  library(RSkittleBrewer)
  library(genefilter)
  library(dplyr)
  library(devtools)
  pheno_data = read.table("sample.csv",header=T,sep=",",colClasses=c("character","factor"))
  #fpkm = read.table("FPKM_table.csv",header=T,sep=",")
  label <- read.csv("~/label.csv",header = T,sep=",")
  fpkm <- count[,c(1,12:21)]
  fpkm<-fpkm[which(rowSums(fpkm[,2:11]) > 0),]
  gene <- fpkm[,1]
  fpkm<-fpkm[,-1]
  fpkm = log10(fpkm+1)
  #fpkm = log2(fpkm)
  pdf("box1.pdf")
  boxplot(fpkm,las=2,ylab='log10(FPKM+1)', col=rainbow(12),cex.axis=0.5)
  dev.off()
  
  
  library(ggplot2)
  library(reshape2)
  fpkm <- cbind(gene,fpkm)
  data_m <- melt(fpkm)
  head(data_m)

## 可以利用strsplit分割，取出其前面的字符串
## R中复杂的输出结果多数以列表的形式体现，在之前的矩阵操作教程中
## 提到过用str函数来查看复杂结果的结构，并从中获取信息
  group = unlist(lapply(substr(data_m$variable,6,8), function(x) x[1]))
  data_m$group = group
  pdf("boxplot.pdf")
  ggplot(data_m, aes(x=variable, y=value),color=variable) + 
    geom_boxplot(aes(fill=factor(group))) + labs(x="Sample",y="log10(FPKM+1)") + 
    theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) 
  dev.off()

# 3.DEGs做pca
  library(DESeq2)
  diff<-rbind(D10,D12)
  differ<-count[match(diff[,1],count[,1]),]
  differ<-unique(differ)
  differ<-differ[,1:11]
  rownames(differ) <- differ[,1]
  differ <- differ[,-1]
  condition <- factor(c(rep("D10", 6), rep("D12", 4)), levels = c("D10","D12"))
  colData <- data.frame(row.names = colnames(differ), condition)
  dds <- DESeqDataSetFromMatrix(differ, colData, design = ~condition)
#归一化，因为样本量超过了30，因此用vst
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
#返回样本名和分组
  pcaData <- plotPCA(vsd, intgroup=c("condition"),returnData = T)
#这里按照condition排序了，原因见下。
  pcaData <- pcaData[order(pcaData$condition,decreasing=F),]
#知道每一个组有多少样本
  table(pcaData$condition)
# plot
  library(ggrepel)
  library(ggplot2)
  pca<-ggplot(pcaData, aes(PC1, PC2, color=group, shape=group)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",pcaData[,1],"% variance")) +
    ylab(paste0("PC2: ",pcaData[,2],"% variance")) +
    coord_fixed()+
    geom_text_repel(label = pcaData$name)
  
  ggsave("PCA.png",pca)
  dev.off()
  
####################################################################################################
# RNA-seq的counts值，RPM, RPKM, FPKM, TPM 的异同 https://cloud.tencent.com/developer/article/1484078
# Counts FPKM RPKM TPM CPM 的转化 https://cloud.tencent.com/developer/article/2018630
# 获取基因有效长度的N种方法 https://cloud.tencent.com/developer/article/2031969
