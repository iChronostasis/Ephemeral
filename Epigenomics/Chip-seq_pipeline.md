# Chip-seq 分析流程
<mark>ATAC-seq和Chip-seq分析流程基本一致，具体个性化分析可以参考CUT_Tag_pipeline的部分内容。</mark>

<!-- TOC -->

- [Chip-seq 分析流程](#chip-seq-分析流程)
  - [（一）原理](#一原理)
  - [（二）数据下载](#二数据下载)
  - [（三）分析流程](#三分析流程)
    - [1. 质控fastqc](#1-质控fastqc)
    - [2. 过滤](#2-过滤)
    - [3.  比对](#3--比对)
      - [3.1 构建索引](#31-构建索引)
      - [3.2 比对](#32-比对)
    - [4. 去除PCR重复](#4-去除pcr重复)
    - [5. call peak](#5-call-peak)
    - [6. peak注释（ChIPseeker）](#6-peak注释chipseeker)
    - [7. motif分析](#7-motif分析)
  - [参考文章：](#参考文章)

<!-- /TOC -->

<mark>**进行分析前需要先安装相关的软件：trimmomatic，bowtie2，GATK，MACS2**</mark>

## （一）原理
ChIP被用来研究细胞内DNA与蛋白质相互作用，具体来说就是确定特定蛋白（如转录因子）是否结合特定基因组区域（如启动子或其它DNA结合位点）——可能定义顺反组。ChIP还被用来确定基因组上与组蛋白修饰相关的特定位点（即组蛋白修饰酶类的靶标）。
**流程：**
* 第一步： 将蛋白交联到DNA上。 也就是保证蛋白和DNA能够结合，找到互作位点。
* 第二步： 通过超声波剪切DNA链。
* 第三步： 加上附上抗体的磁珠用于免疫沉淀靶蛋白。抗体很重要
* 第四步： 接触蛋白交联；纯化DNA
![](Epigenomics/images/Chip_seq.png)
## （二）数据下载
一个例子...
实际应用时进行更改...
```bash
module load sratoolkit/v2.11.3
# download
mkdir -p /home/User/ChIP-seq && cd /home/User/ChIP-seq  
mkdir rawData && cd rawData  
prefetch SRR8073294

# 批量下载
prefetch --option-file SRR_Acc_List.txt

# .sra to .fastq
mkdir cleanData && cd cleanData

# a sample 
fastq-dump --gzip -O ~/ChIP-seq/cleanData SRR8073294.sra

# many samples
for i in SRR*
do
echo $i
fastq-dump --gzip --split-files -O /home/User/ChIP-seq/cleanData $i/$i.sra
done
```

## （三）分析流程
### 1. 质控fastqc

```bash
module load conda/miniconda3
module load Java

# 原始数据存放路径：/home/User/ChIP-seq/rawData
# -t 线程数 
# -o 输出路径
Work_directory=/home/User/ChIP-seq
mkdir $Work_directory/result
mkdir $Work_directory/result/1.fastqc
#find /home/User/ChIP-seq/cleanData/ -name *fq.gz| xargs fastqc -t 20 -o 
find /home/User/ChIP-seq/cleanData/ -name *fastq.gz| xargs fastqc -t 20 -o $Work_directory/result/1.fastqc/
## 用multiqc整合报告
pip install git+https://github.com/ewels/MultiQC.git
pip install --user multiqc
cd $Work_directory/result/1.fastqc
multiqc *zip -o multiqc/
```

![](Epigenomics/images/Chip_seq1.png)

### 2. 过滤

```bash
Work_directory=/home/User/ChIP-seq
cd $Work_directory/cleanData
mkdir $Work_directory/result/2.Trim
# Trimmomatic
for i in *_1.fastq.gz; 
do
i=${i%_1.fastq.gz*}; 
trimmomatic PE -baseout $Work_directory/result/2.Trim/${i}.fq.gz ${i}_1.fastq.gz ${i}_2.fastq.gz \
ILLUMINACLIP:/share/old_apps/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:8
done

#或者用trim_galore过滤
# for i in *_1.fastq.gz; 
# do
# i=${i%_1.fastq.gz*}; 
# trim_galore --phred33 -q 5 --length 36 --stringency 3 --fastqc -o ./ $Work_directory/result/2.trimmed_fq/${i}.fastq.gz
# done
```

![](Epigenomics/images/Chip_seq2.png)

### 3.  比对
#### 3.1 构建索引
* UCSC官网
```bash
mkdir -p /home/User/refdata/bowtie2/ && cd /home/User/refdata/bowtie2/
# hg19
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz  

gunzip hg19.fa.gz 
# mm10
wget -c http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz 
 
##构建索引  
nohup bowtie2-build --threads 2 -f hg19.fa hg19 &

```

* [bowtie2官网](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
![](Epigenomics/images/bowtie2_install_ref.jpg)

```bash
# install
conda install -c bioconda bowtie2


mkdir refdata
cd refdata
unzip -o /home/User/refdata/mm10 mm10.zip
rm mm10.zip make_mm10.sh

# index directory
#/home/User/refdata/mm10
```

#### 3.2 比对

```bash
# Run
module load conda/miniconda3
module load Java
module load bowtie/v2.4.4

# 去bowtie2的官网下载相关物种版本的bowtie2_index, 此次mm10 

bowtie2_index='/share/data3/Cuttag_data_wangxinLab/mm10/bowtie/mm10'
Work_directory=/home/User/ChIP-seq
cd $Work_directory/result/2.Trim
mkdir $Work_directory/result/3.Align
mkdir $Work_directory/result/3.Align/mapped_fragments_bam

for i in *_1P.fq.gz; 
  do
  i=${i%_1P.fq.gz*};  
  (bowtie2  --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -p 20 -x $bowtie2_index -1 $Work_directory/result/2.Trim/${i}_1P.fq.gz -2 $Work_directory/result/2.Trim/${i}_2P.fq.gz) 2> $Work_directory/result/3.Align/${i}.bowtie2 | samtools view -bS - >$Work_directory/result/3.Align/${i}_aligned_reads.bam  
  # 保留匹配的reads 
  samtools view -bh -f 3 -F 4 -F 8 $Work_directory/result/3.Align/${i}_aligned_reads.bam > $Work_directory/result/3.Align/mapped_fragments_bam/${i}.step1.bam
done
```

### 4. 去除PCR重复

```bash
module load conda/miniconda3
module load Java
module load GATK/v4.2.2.0

Output_directory=/home/User/ChIP-seq/result/3.Align
File=/home/User/ChIP-seq/SRR_Acc_List.txt
mkdir $Output_directory/dup_marked
mkdir $Output_directory/rmdup

# Run
cat $File | while read LINE || [[ -n ${LINE} ]]
do
  a=$(echo $LINE |sed 's/\n//g')
  sample=$(echo $a |sed 's/\r//g')
  echo $sample 
  
  gatk SortSam -I $Output_directory/mapped_fragments_bam/${sample}.step1.bam -O $Output_directory/mapped_fragments_bam/${sample}.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY SILENT
  
  ## Marking duplicates
  gatk MarkDuplicates -INPUT $Output_directory/mapped_fragments_bam/${sample}.bam -OUTPUT $Output_directory/dup_marked/${sample}.bam -VALIDATION_STRINGENCY SILENT -METRICS_FILE $Output_directory/dup_marked/metrics.${sample}.txt
  
  ## Removing duplicates 2  #REMOVE_DUPLICATES=true
  gatk MarkDuplicates I=$Output_directory/mapped_fragments_bam/${sample}.bam  O=$Output_directory/rmdup/${sample}.bam REMOVE_DUPLICATES=true METRICS_FILE=$Output_directory/rmdup/metrics.${sample}_picard.rmDup.txt
  
  ## Creating bam index files
  ##samtools index $Output_directory/Align/mapped_fragments_bam/${sample}.bam
  samtools index $Output_directory/rmdup/${sample}.bam
  #### Convert into bigwig file format
  bamCoverage -b $Output_directory/rmdup/${sample}.bam -o $Output_directory/rmdup/${sample}.bw
  ##bamCoverage -b $Output_directory/dup_marked/${sample}.bam -o $Output_directory/dup_marked/bigwig/${sample}.bw
  ##bamCoverage -b $Output_directory/dup_marked/${sample}.bam --normalizeUsing RPKM -o $Output_directory/dup_marked/bigwigNorm/${sample}.bw 
done

```

### 5. call peak

```bash
module load conda/miniconda3
source activate macs2

Input_directory=/home/User/ChIP-seq/result/3.Align/rmdup
Output_directory=/home/User/ChIP-seq/result

mkdir $Output_directory/4.peak
mkdir $Output_directory/4.peak/filter_q0.05
mkdir $Output_directory/4.peak/filter_p0.05

# Run
cd $Input_directory
for i in *.bam; 
do
  i=${i%.bam*}; 
  # q_value=0.05
  mkdir $Output_directory/4.peak/filter_q0.05/${i}
  macs2 callpeak -t $Input_directory/${i}.bam -g mm -f BAMPE -n ${i} --outdir $Output_directory/4.peak/filter_q0.05/${i} -q 0.05 -B --SPMR --keep-dup all 2> $Output_directory/4.peak/filter_q0.05/${i}/${i}.macs2
  ### add control "-c"
  # p_value=0.05
  mkdir $Output_directory/4.peak/filter_p0.05/${i}
  macs2 callpeak -t $Input_directory/${i}.bam -g mm -f BAMPE -n ${i} --outdir $Output_directory/4.peak/filter_p0.05/${i} -p 0.05 -B --SPMR --keep-dup all 2> $Output_directory/4.peak/filter_p0.05/${i}/${i}.macs2
done
# Finished
conda deactivate
```

**macs2的输出文件解读**

得到的文件有：

**NAME_peaks.xls**

虽然后缀名是xls，但实际上就是一个普通的文本文件。包含peak信息的tab分割的文件，前几行会显示callpeak时的命令。输出信息包含：

- 染色体号
- peak起始位点
- peak结束位点
- peak区域长度
- peak的峰值位点（summit position）
- peak 峰值的属性（包括pileup峰高和可信度）（**pileup height** at peak summit, -log10(pvalue) **for the peak summit**）都是值越高越好
- peak的富集倍数（相对于random Poisson distribution with local lambda）

### 6. peak注释（ChIPseeker）
[ChIPseeker](https://www.plob.org/article/24683.html)虽然最初是为了ChIP-seq注释而写的一个R包，但它不只局限于ChIP-seq，也可用于ATAC-Seq等其他富集peaks注释，也可用于lincRNA注释、DNA breakpoints的断点位置注释等所有genomic coordination的注释，另外提供了丰富的可视化功能。
**使用方法：**
使用ChIPseeker需要准备两个文件：
-   要注释的peaks的文件，需满足BED格式。
-   注释参考文件，即需要一个包含注释信息的TxDb对象
-   **从UCSC下载参考基因信息做注释，从UCSC生成TxDb**

```R
## install the packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPseeker")
## loading packages
library(ChIPseeker)
library(GenomicFeatures)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(clusterProfiler)

WT="/home/～_peaks.narrowPeak"
KO="/home/～_peaks.narrowPeak"
peaks=c(WT,KO)
outname=c("WT","KO")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
head(seqlevels(txdb))

## read the bed files
peak=list()
for(i in 1:length(peaks)){
  peak[[i]]=readPeakFile(peakfile=peaks[i], header=F, as = 'GRanges')
}
head(seqlevels(peak[[1]]))

## annotatePeak
peakAnno=list()
for(i in 1:length(peak)){
  peakAnno[[i]]=annotatePeak(peak = peak[[i]], tssRegion=c(-3000,3000), TxDb = txdb, annoDb='org.Mm.eg.db',assignGenomicAnnotation=T,flankDistance=5000)
}####tssRegion=c(-3000,3000):想看peak上下游某个范围内（比如说-3k到3k的距离）都有什么基因

## save the result
peakAnno_df=list()
for(i in 1:length(peakAnno)){
  peakAnno_df[[i]]=as.data.frame(peakAnno[[i]])
}

for(i in 1:length(peakAnno_df)){
  peakAnno_df[[i]]$strand=NULL
  colnames(peakAnno_df[[i]])[1:11] = c("chr","start","end","width","peak_name","score","strand","fold_change","-log10(pvalue)","-log10(qvalue)","peak_value")
}

for(i in 1:length(peakAnno_df)){
  write.table(peakAnno_df[[i]],paste(outname[i], sep="_","peaks_anotate.txt"),r=F,c=T,sep="\t",quote=F)
}

for(i in 1:length(peakAnno)){
	pdf(file=paste0(outname[i],".peakAnnotation.pdf"))
	plotAnnoPie(peakAnno[[i]], main=paste0(outname[i],"\nDistribution of Peaks"), line=-8)
	dev.off()
}
```


**输出文件解读:**
1.  genomic annotation注释：`annotation列`，注释的是peak的位置，它落在什么地方了，可以是UTR，可以是内含子或者外显子    
2.  nearest gene annotation（最近基因）注释：`多出来的其它列`，注释的是peak相对于转录起始位点（TSS）的距离，不管这个peak是落在内含子或者别的什么位置上，即使它落在基因间区上，都能够注释到一个离它最近的基因（即使它可能非常远）
	-   针对于不同的转录本，一个基因可能有多个转录起始位点，所以注释是在转录本的水平上进行的，可以看到输出有一列是**transcriptId**
3.  两种注释策略面对不同的研究对象，**第一种策略，peak所在的位置可能就是调控本身**（比如你要做可变剪切的，最近基因的注释显然不是你关注的点）；而做基因表达调控的，当然**promoter区域是重点，离结合位点最近的基因更有可能被调控**。

![](Epigenomics/images/Chip_seq3.png)

### 7. motif分析
Homer软件集成了许多的功能，包括peak calling, peak注释，motif分析等等，通过这一个软件，就可以完成chip_seq的绝大部分分析内容。
```bash

# 默认安装 homer 最新版

conda install -c bioconda homer -y # -y 表示直接安装，不询问确认与否

# 下载配置文件

mkdir Homer

cd Homer

wget http://homer.ucsd.edu/homer/configureHomer.pl

perl configureHomer.pl -install homer

```

在homer中通过`annotatePeaks.pl`这个脚本进行peak的注释，分为以下两步：

* **准备参考基因组的注释信息**

```bash

perl configureHomer.pl -list

  

perl configureHomer.pl -install mm10

```

* 寻找motif
```bash
# 利用homer进行寻找motif：
# 因为homer对输入文件有格式要求，所以需要对call peaks的输出文件进行格式转换：
# awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' Xu_MUT_rep2_rmdup_summits.bed >homer_peaks.tmp 

Work_directory=/home/User/ChIP-seq
File=/home/User/ChIP-seq/SRR_Acc_List.txt
mkdir $Work_directory/result/bedfile
cd $Work_directory/result/4.peak/filter_q0.05

cat $File | while read LINE || [[ -n ${LINE} ]]
do
	a=$(echo $LINE |sed 's/\n//g')
	sample=$(echo $a |sed 's/\r//g')
	echo $sample 
	cp $Work_directory/result/4.peak/filter_q0.05/$sample/$sample\_summits.bed $Work_directory/result/bedfile
done

mkdir $Work_directory/result/motif
cd $Work_directory/result/motif
cat $File | while read LINE || [[ -n ${LINE} ]]
do
	a=$(echo $LINE |sed 's/\n//g')
	sample=$(echo $a |sed 's/\r//g')
	echo $sample 
	awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' $Work_directory/result/bedfile/$sample\_summits.bed>homer_peaks.tmp  
	findMotifsGenome.pl homer_peaks.tmp mm10 ${sample}_motifDir -p 10 -len 8,10,12 2>${sample}.motif.log 
	### 注释
	#annotatePeaks.pl $Work_directory/result/4.peak/filter_p0.05/${i}/${i}_peaks.narrowPeak mm10 >${sample}_peak_annotation.txt
	annotatePeaks.pl homer_peaks.tmp mm10 1>${sample}.peakAnn.xls 2>${sample}.ann_err.log
done

## -p：线程
## -len：motif 长度, 默认=8,10,12

```

homer输出结果：

homer找motif结果文件：

![](Epigenomics/images/Chip_seq4.png)


## 参考文章：
1. https://mp.weixin.qq.com/s/VNYNGj8LbJnzcELsMZtOjg
2. https://mp.weixin.qq.com/s/Ole5PPTb4nfUFYGKS9bmgg
3. http://homer.ucsd.edu/homer/ngs/annotation.html
_____

