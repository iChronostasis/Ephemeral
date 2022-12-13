# RNA-seq pipeline

## 分析流程概览及所需软件

> 1. 上游分析

创建单独的分析环境：miniconda3

数据下载和格式转换：sra-tools(sratoolskit)

质控清洗：fastqc，multiqc，trim-galore，trimmomatic，cutadapt，fastp

序列比对：hisat2，subread，star，bwa，bowtie2，tophat

计数：featurecounts，htseq，bedtools，deeptools

> 2. 下游分析

DEG差异分析：DESeq2，edgeR，limma

差异基因富集注释：

（1）KEGG 基因通路分析 与 GO 基因功能分析

（2）GSEA 基因集富集分析，GSVA 基因集变异分析

（3）PPI 蛋白质互作网络

（4）WGCNA 加权网络共表达网络分析

## 1. 参考基因索引

### 1.1 下载基因组

### 1.1.1 小鼠参考基因组

```sh
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz
```

### 1.1.2 人类参考基因组

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz

wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
```

## 1.2 建立索引

**Reads=150bp**

```sh
STAR \
--runMode genomeGenerate \
--genomeDir /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/reads_150bp \
--genomeFastaFiles /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz \
--sjdbGTFfile /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz \
--runThreadN 8 \
--sjdbOverhang 149   ###150bp
```

## 2. 检测数据质量(fastqc)

```sh
module load conda

File1=/share/data3/xmu_chenxiaofen/4atask/CleanData
#File2=/home/chengxin/RNA-seq_slurm/List.txt
File3=/home/chengxin/RNA-seq_slurm/result

for i in *_1.fastq.gz; 
do
i=${i%_1.fastq.gz*}; 
cd $File3/01.fastqc
mkdir $i
fastqc --t 10 -o $File3/01.fastqc/$b $File1/$i/$i\_1.fq.gz $File2/$i/$i\_2.fq.gz
done
```

## 3. 比对fastq

```sh
### 加载软件
module load conda
File1=/share/data3/xmu_chenxiaofen/4atask/CleanData
File2=/home/chengxin/RNA-seq_slurm/List.txt
File3=/home/chengxin/RNA-seq_slurm/result
cd /share/Projects/chengxin/02.align

for i in *_1.fastq.gz; 
do
i=${i%_1.fastq.gz*}; 
STAR \
--genomeDir /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/reads_150bp \
--readFilesIn /share/data3/xmu_chenxiaofen/4atask/CleanData/$i/$i\_1.fq.gz /share/data3/xmu_chenxiaofen/4atask/CleanData/$i/$i\_2.fq.gz \
--readFilesCommand zcat \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix $i
done
```

## 4. 定量featureCounts

```sh
module load conda

File1=/share/Projects/chengxin/02.align

featureCounts \
	-T 16 \
	-p \
	-t exon \
	-g gene_id \
	-a /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/reads_150bp/GCF_000001635.26_GRCm38.p6_genomic.gtf \
	-o /share/Projects/chengxin/out_feature_.txt \
	$File1/*Aligned.sortedByCoord.out.bam 
```

## 5. DESeq2

该步骤在R中运行，参照该包的Bioconductor的官网流程。

https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

```r
rm(list=ls())
library(DESeq2)
setwd("/Users/hecate/研一/RNA-seq/DESeq2")
data<-read.table("/Users/hecate/研一/RNA-seq/DESeq2/all_feature.txt", header=TRUE, quote="\t")[,c(1,7:20)]
### 质控
data<-data[rowSums(data[,c(2:15)]) > 0, ]
mycounts <- data[,c(2:15)]
#显示mycounts信息
head(mycounts)
rownames(mycounts)<-data$Geneid
#设置样品组别、重复数
# MS4A4A-vs-WT
condition <- factor(c(rep("MS4A4A", 7), rep("WT", 7)), levels = c("MS4A4A","WT"))
#显示condition设置
condition
#设置colData值
colData <- data.frame(row.names = colnames(mycounts), condition)
#显示colData值
colData
#在R里面用于构建公式对象，~左边为因变量，右边为自变量。
#标准流程：dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design= ~ batch + condition) 
#countData为表达矩阵即countdata
#colData为样品信息矩阵即coldata
#design为差异表达矩阵即批次和条件（对照、处理）等
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)
#对原始dds进行normalize
dds <- DESeq(dds)
#显示dds信息
dds

#使用DESeq2包中的results()函数，提取差异分析的结果
#Usage:results(object, contrast, name, .....）
#将提取的差异分析结果定义为变量"res" 
#contrast: 定义谁和谁比较
res = results(dds, contrast=c("condition","MS4A4A","WT"))

#对结果res利用order()函数按pvalue值进行排序
#创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)
# lfcShrink 不改变p值q值，但改变了fc，使 foldchange范围变小，
#所以选择DEG时会有不同结果，一般会偏少！所以，根据数据情况，本次分析DEG还是不做shrink
#resultsNames(dds)
#resLFC <- lfcShrink(dds, coef="condition_WT_vs_MS4A4A", type="apeglm")
#将分析的所有结果进行输出保存
#write.csv(res, file="All_results.csv")
#显示显著差异的数目
table(res$padj<0.05)
table(res$log2FoldChange>=0&res$padj<0.05)
table(res$log2FoldChange<0&res$padj<0.05)

table(res$padj<0.1)
table(res$log2FoldChange>=0&res$padj<0.1)
table(res$log2FoldChange<0&res$padj<0.1)

table(res$padj<0.2)
table(res$log2FoldChange>=0&res$padj<0.2)
table(res$log2FoldChange<0&res$padj<0.2)

res<-as.matrix(res[,c(2,5,6)])
colnames(res)<-c("MS4A4A-vs-WT_log2FC","MS4A4A-vs-WT_p_value","MS4A4A-vs-WT_FDR")
all_result<-res


# WT_F-vs-WT_M
condition <- factor(c(rep("WT_F", 3), rep("WT_M", 4)), levels = c("WT_F","WT_M"))
mycounts <- data[,c(9:15)]
#显示mycounts信息
head(mycounts)
rownames(mycounts)<-data$Geneid
colData <- data.frame(row.names = colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)

#对原始dds进行normalize
dds <- DESeq(dds)

res = results(dds, contrast=c("condition","WT_F","WT_M"))
#对结果res利用order()函数按pvalue值进行排序
#创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)

volcano(DEG_data=res,taskname="WT_F-vs-WT_M",dir_path="/Users/hecate/研一/RNA-seq/DESeq2/")
#将分析的所有结果进行输出保存
#write.csv(res, file="All_results.csv")
#显示显著差异的数目
table(res$padj<0.05)
table(res$log2FoldChange>=0&res$padj<0.05)
table(res$log2FoldChange<0&res$padj<0.05)

table(res$padj<0.1)
table(res$log2FoldChange>=0&res$padj<0.1)
table(res$log2FoldChange<0&res$padj<0.1)

table(res$padj<0.2)
table(res$log2FoldChange>=0&res$padj<0.2)
table(res$log2FoldChange<0&res$padj<0.2)

res<-as.matrix(res[,c(2,5,6)])
colnames(res)<-c("WT_F-vs-WT_M_log2FC","WT_F-vs-WT_M_p_value","WT_F-vs-WT_M_FDR")
all_result<-cbind(all_result,res)

# MS4A4A_F-vs-WT_F
condition <- factor(c(rep("MS4A4A_F", 4), rep("WT_F", 3)), levels = c("MS4A4A_F","WT_F"))
mycounts <- data[,c(2:5,9:11)]
#显示mycounts信息
head(mycounts)
rownames(mycounts)<-data$Geneid
colData <- data.frame(row.names = colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)

#对原始dds进行normalize
dds <- DESeq(dds)
res = results(dds, contrast=c("condition","MS4A4A_F","WT_F"))
#对结果res利用order()函数按pvalue值进行排序
#创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)

#显示显著差异的数目
table(res$padj<0.05)
table(res$log2FoldChange>=0&res$padj<0.05)
table(res$log2FoldChange<0&res$padj<0.05)

table(res$padj<0.1)
table(res$log2FoldChange>=0&res$padj<0.1)
table(res$log2FoldChange<0&res$padj<0.1)

table(res$padj<0.2)
table(res$log2FoldChange>=0&res$padj<0.2)
table(res$log2FoldChange<0&res$padj<0.2)

res<-as.matrix(res[,c(2,5,6)])
colnames(res)<-c("MS4A4A_F-vs-WT_F_log2FC","MS4A4A_F-vs-WT_F_p_value","MS4A4A_F-vs-WT_F_FDR")
all_result<-cbind(all_result,res)

# MS4A4A_F-vs-MS4A4A_M
condition <- factor(c(rep("MS4A4A_F", 4), rep("MS4A4A_M", 3)), levels = c("MS4A4A_F","MS4A4A_M"))
mycounts <- data[,c(2:8)]
#显示mycounts信息
head(mycounts)
rownames(mycounts)<-data$Geneid
colData <- data.frame(row.names = colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)

#对原始dds进行normalize
dds <- DESeq(dds)
res = results(dds, contrast=c("condition","MS4A4A_F","MS4A4A_M"))
#对结果res利用order()函数按pvalue值进行排序
#创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)
volcano(DEG_data=res,taskname="MS4A4A_F-vs-MS4A4A_M",dir_path="/Users/hecate/研一/RNA-seq/DESeq2/")

#显示显著差异的数目
table(res$padj<0.05)
table(res$log2FoldChange>=0&res$padj<0.05)
table(res$log2FoldChange<0&res$padj<0.05)

table(res$padj<0.1)
table(res$log2FoldChange>=0&res$padj<0.1)
table(res$log2FoldChange<0&res$padj<0.1)

table(res$padj<0.2)
table(res$log2FoldChange>=0&res$padj<0.2)
table(res$log2FoldChange<0&res$padj<0.2)

res<-as.matrix(res[,c(2,5,6)])
colnames(res)<-c("MS4A4A_F-vs-MS4A4A_M_log2FC","MS4A4A_F-vs-MS4A4A_M_p_value","MS4A4A_F-vs-MS4A4A_M_FDR")
all_result<-cbind(all_result,res)

# MS4A4A_M -vs-WT_M
condition <- factor(c(rep("MS4A4A_M", 3), rep("WT_M", 4)), levels = c("MS4A4A_M","WT_M"))
mycounts <-data[,c(6:8,12:15)]
#显示mycounts信息
head(mycounts)
rownames(mycounts)<-data$Geneid
colData <- data.frame(row.names = colnames(mycounts), condition)
dds <- DESeqDataSetFromMatrix(mycounts, colData, design = ~condition)

#对原始dds进行normalize
dds <- DESeq(dds)
res = results(dds, contrast=c("condition","MS4A4A_M","WT_M"))
#对结果res利用order()函数按pvalue值进行排序
#创建矩阵时，X[i,]指矩阵X中的第i行，X[,j]指矩阵X中的第j列
#order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
res = res[order(res$pvalue),]
#显示res结果首信息
head(res)
#对res矩阵进行总结，利用summary命令统计显示一共多少个genes上调和下调
summary(res)

#显示显著差异的数目
table(res$padj<0.05)
table(res$log2FoldChange>=0&res$padj<0.05)
table(res$log2FoldChange<0&res$padj<0.05)

table(res$padj<0.1)
table(res$log2FoldChange>=0&res$padj<0.1)
table(res$log2FoldChange<0&res$padj<0.1)

table(res$padj<0.2)
table(res$log2FoldChange>=0&res$padj<0.2)
table(res$log2FoldChange<0&res$padj<0.2)

res<-as.matrix(res[,c(2,5,6)])
colnames(res)<-c("MS4A4A_M -vs-WT_M_log2FC","MS4A4A_M -vs-WT_M_p_value","MS4A4A_M -vs-WT_M_FDR")
all_result<-cbind(all_result,res)

#将分析的所有结果进行输出保存
write.csv(all_result,"all_result.csv")
```



## 6. 后续的流程分析

```R
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

#————————————————————————————————————————
# if error,then do following command to delet 0 row
  which(apply(count,1,sum)==0)
  count1<-count[-which(apply(count,1,sum)==0),]
  tcount<-t(count1)
  which(apply(tcount,2,sum)==0)
  pca<- prcomp(tcount, scale. = TRUE)
# if still error, then scale tcount to find out NA and delet them
  x <- tcount
  x.scale <- scale(x)
  which(is.na(apply(x.scale,2,var))==TRUE)
  rname <- names(which(is.na(apply(x.scale,2,var))==TRUE))
  x.save <- which(colnames(x) %in% rname == FALSE)
  x.filter <- x[,x.save]
  prcomp(x.filter,scale. = TRUE)

  xlab <- paste("PC1","(",round((summary(pca))$importance[2,1]*100,1),"%)",sep="")
  ylab <- paste("PC2","(",round((summary(pca))$importance[2,2]*100,1),"%)",sep="")
  x<-"PC1"
  y<-"PC2"
  df <- data.frame(pc1 = pca$x[,1],
                  pc2 = pca$x[,2])
  df <- cbind(df,col)

      ggplot(data=df,aes(x =pc1,y=pc2 ,color=condition))+
        geom_point(size=3)+
        theme_bw()+theme(panel.grid=element_blank())+
        theme(legend.position = c(0.80,0.9))+stat_ellipse(lwd=0.5)

#————————————————————————————————————————
### 查看某一个感兴趣的gene在组间的差别(Before)
  plotCounts(dds, gene="ENSMUSG00000024045", intgroup="condition", returnData=TRUE) %>% ggplot(aes(condition, count)) + geom_boxplot(aes(fill=condition)) + scale_y_log10() + ggtitle("selected gene")



```



