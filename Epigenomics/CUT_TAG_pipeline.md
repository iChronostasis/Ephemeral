<!-- TOC -->

- [CUT_Tag_pipeline](#cut_tag_pipeline)
  - [（一）原理](#一原理)
  - [（二）分析概要](#二分析概要)
  - [（三）数据处理](#三数据处理)
    - [1. 数据预处理  Data Pre-processing](#1-数据预处理--data-pre-processing)
      - [1.1 Quality Control using FastQC [Optional]](#11-quality-control-usingfastqcoptional)
      - [1.2 合并多个泳道 Merge technical replicates/lanes if needed [Optional]](#12-合并多个泳道-merge-technical-replicateslanes-if-needed-optional)
    - [2.  过滤数据(filter data)](#2--过滤数据filter-data)
    - [3. 比对 Alignment](#3-比对-alignment)
      - [3.1. 构建索引](#31-构建索引)
      - [3.2 比对](#32-比对)
        - [3.2.1 与参考基因组对齐](#321-与参考基因组对齐)
        - [3.2.2. 与spike-in基因组对齐以进行spike-in校准[可选/推荐]](#322-与spike-in基因组对齐以进行spike-in校准可选推荐)
      - [3.3 报告排序映射汇总【必填】](#33-报告排序映射汇总必填)
        - [1. 测序深度](#1-测序深度)
        - [2. Spike-in alignment](#2-spike-in-alignment)
      - [3.4 Filter BAM(Samtools)](#34-filter-bamsamtools)
      - [3.4.1 去除PCR重复](#341-去除pcr重复)
      - [3.4.2 评估映射片段大小分布[required]](#342-评估映射片段大小分布required)
      - [3.4.3 评估实验重复性](#343-评估实验重复性)
        - [3.4.3.1 文件格式转换](#3431-文件格式转换)
        - [3.4.3.2. 评估重复实验的重现性 Assess replicate reproducibility](#3432-评估重复实验的重现性-assess-replicate-reproducibility)
          - [1.  使用R语言进行评估重复重现性](#1--使用r语言进行评估重复重现性)
          - [2. deeptools 相关性分析](#2-deeptools-相关性分析)
      - [3.5  Spike-in calibration[optional]](#35--spike-in-calibrationoptional)
        - [3.5.1  Scaling factor](#351--scaling-factor)
      - [3.6 Peak Calling](#36-peak-calling)
        - [3.6.1  MACS2](#361--macs2)
        - [3.6.2 SEACR](#362-seacr)
          - [3.6.2.1 使用SEACR进行peak calling](#3621-使用seacr进行peak-calling)
          - [3.6.2.2  Number of peaks called](#3622--number-of-peaks-called)
          - [3.6.2.3 生物重复峰的重现性 Reproducibility of the peak across biological replicates](#3623-生物重复峰的重现性-reproducibility-of-the-peak-across-biological-replicates)
          - [3.6.2.4 峰区的片段比例（FRiPs）](#3624-峰区的片段比例frips)
          - [3.6.2.5 峰数、峰宽、峰重现性和 FRiP 的可视化](#3625-峰数峰宽峰重现性和-frip-的可视化)
      - [3.7. 可视化](#37-可视化)
        - [3.7.1 IGV](#371-igv)
        - [3.7.2 特定区域的热图可视化](#372-特定区域的热图可视化)
          - [3.7.2.1 转录单位的热图](#3721-转录单位的热图)
          - [3.7.2.2 CUT&Tag 峰的热图](#3722-cuttag-峰的热图)
      - [3.8 Peak Annotation](#38-peak-annotation)
        - [3.8.1  ChIPseeker 进行峰注释](#381--chipseeker-进行峰注释)
        - [3.8.2  超几何检验](#382--超几何检验)
      - [3.9 差异peaks分析](#39-差异peaks分析)
        - [3.9.1 DESeq2](#391-deseq2)
        - [3.9.2 差异peaks分析--DiffBind](#392-差异peaks分析--diffbind)
    - [参考文章：](#参考文章)

<!-- /TOC -->

# CUT_Tag_pipeline

## （一）原理
真核细胞核中 DNA 上发生的所有动态过程都发生在染色质景观的背景下，该景观包括核小体及其修饰、转录因子和染色质相关复合物。多种染色质特征标记激活和沉默转录调控元件和染色质结构域的位点，这些位点在细胞类型之间有所不同，并在发育过程中发生变化。
全基因组染色质特征图谱传统上是使用**染色质免疫沉淀 (ChIP) 进行的，其中染色质被交联和溶解，蛋白质的抗体或感兴趣的修饰用于免疫沉淀结合的 DNA（图 1a ）**。另一种染色质分析策略是原位酶系，其中染色质蛋白或感兴趣的修饰被抗体或融合蛋白靶向。然后，潜在的 DNA 被标记或切割，并且在过去的二十年中引入了一系列酶系方法。
<mark>目标和标记下的切割 (CUT&Tag) 是一种使用蛋白-A-Tn5 (pA-Tn5) 转座融合蛋白的束缚方法（图 1b）。</mark>在 CUT&Tag 中，透化的细胞或细胞核与特定染色质蛋白的抗体一起孵育，然后装载有马赛克末端接头的 pA-Tn5 连续连接到抗体结合位点。通过添加镁离子激活转座体导致接头整合到附近的 DNA 中。然后将它们放大以生成测序文库。由于 pA-Tn5 束缚后样品的严格洗涤和适配器整合的高效率，基于抗体束缚 Tn5 的方法实现了高灵敏度。**相对于 ChIP-seq 改进的信噪比意味着映射染色质特征所需的测序量减少了一个数量级。**

![](Epigenomics/images/ChIPseqCUTTag1.png)
**图 1. 免疫沉淀和抗体靶向染色质分析策略之间的差异。** **A.** ChIP-seq 实验程序。**B.** CUT&Tag 实验程序。细胞和细胞核以灰色表示，染色质以红色核小体表示，特定染色质蛋白以绿色表示。

## （二）分析概要

![](Epigenomics/images/CUT_TAG_pipeline.png)


**官网数据：**
```bash
# 直接下载 or 类似ChIP-seq和ATAC-seq的数据下载使用SRA Toolkit
Work_directory=/home/User/CUT_Tag
cd $Work_directory/cleanData

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_1.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_2.fastq.gz



# 下载时根据实验设计更改fastq文件名字
projPath=/home/User/CUT_Tag
mkdir $projPath/data && mkdir $projPath/data/IgG_rep2
wget -O $projPath/data/IgG_rep2/IgG_rep2_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_1.fastq.gz

wget -O $projPath/data/IgG_rep2/IgG_rep2_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_2.fastq.gz

wget -O $projPath/data/IgG_rep2/IgG_rep2_R1_002.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_1.fastq.gz

wget -O $projPath/data/IgG_rep2/IgG_rep2_R2_002.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_2.fastq.gz

```
## （三）数据处理
### 1. 数据预处理  Data Pre-processing
#### 1.1 Quality Control using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) [Optional]
```bash
##== linux command ==##
## if you don't have the FastQC,download the software
mkdir -p /home/User/tools
wget -P /home/User/tools https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
cd /home/User/tools
unzip fastqc_v0.11.9.zip
#############################


### Quality Control
module load conda/miniconda3
module load Java
Work_directory=/home/User/CUT_Tag
mkdir $Work_directory/result
mkdir $Work_directory/result/1.QC

find $Work_directory/cleanData -name *fastq.gz| xargs fastqc -t 20 -o $Work_directory/result/1.QC/
```

![](Epigenomics/images/R1_sequenceContent.png)

结果解读：
**reads开头的序列内容不一致是CUT&Tag reads的普遍现象。未能通过 <mark>Per base seuqnence </mark>内容并不意味着您的数据失败。**
-   这可能是由于 Tn5 偏好。
-   可能会检测到 10 bp 的周期性，它在长度分布中显示为锯齿模式。如果是这样，这是正常的，不会影响对齐或峰值调用。在任何情况下，我们都不建议修剪，因为我们列出的 bowtie2 参数将在不修剪的情况下提供准确的映射信息。

#### 1.2 合并多个泳道 Merge technical replicates/lanes if needed [Optional] 
有时，为了提高效率，样本通常会跨多个泳道进行测序，并且可以在比对之前合并。如果要检查同一样本不同lane的序列之间的重现性，可以跳过这一步，分别对齐每个测序文件（fastq文件）。
```bash
### for example
histName="K27me3_rep1"

mkdir -p ${projPath}/fastq
cat ${projPath}/data/${histName}/*_R1_*.fastq.gz >${projPath}/fastq/${histName}_R1.fastq.gz
cat ${projPath}/data/${histName}/*_R2_*.fastq.gz >${projPath}/fastq/${histName}_R2.fastq.gz
```

### 2.  过滤数据(filter data)
```bash
Work_directory=/home/User/CUT_Tag
cd $Work_directory/cleanData
mkdir $Work_directory/result/2.trimmed_fq
# Trimmomatic
for i in *_1.fastq.gz; 
do
i=${i%_1.fastq.gz*}; 
trimmomatic PE  $Work_directory/result/2.trimmed_fq/${i}.fq.gz ${i}_1.fastq.gz ${i}_2.fastq.gz \
ILLUMINACLIP:/share/old_apps/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true \
LEADING:3 \
TRAILING:3 \
SLIDINGWINDOW:4:15 \
MINLEN:8
done

## ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#或者用trim_galore过滤
# for i in *_1.fastq.gz; 
# do
# i=${i%_1.fastq.gz*}; 
# trim_galore --phred33 -q 5 --length 36 --stringency 3 --fastqc -o ./ $Work_directory/result/2.trimmed_fq/${i}.fastq.gz
# done
```


### 3. 比对 Alignment
带有 Tn5 接头和条形码 PCR 引物的 CUT&Tag 插入文库的结构如下所示：

![](Epigenomics/images/BarcodedCnTLibrary.png)

标准流程是在单个 HiSeq 2500 流动槽上对多达 90 个混合样本执行单索引 25x25 PE Illumina 测序，其中每个样本都有一个独特的 PCR 引物条形码。每个文库的数量都经过调整，以提供约 500 万条双末端读数，从而为丰富的染色质特征提供高质量的分析，并具有特异性和高产量的抗体。不太丰富的特征通常需要更少的读取，而质量较低的抗体可能会增加生成稳健染色质谱所需的读取数量。关于 CUT&Tag 的特征召回和测序深度的全面讨论已经发表（Kaya-Okur 等人 2020）。
#### 3.1. 构建索引
* [bowtie2官网](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

![](Epigenomics/images/bowtie2_install_ref.jpg)

```bash
# e.g.mm10
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
其中以下注释的代码在 **4.3.1 文件格式转换** 会运行，该部分代码是将sam文件转换为后期可以进行评估实验的重复性。

##### 3.2.1 与参考基因组对齐
```bash
# Run
module load conda/miniconda3
module load Java
module load bowtie/v2.4.4
# 去bowtie2的官网下载相关物种版本的bowtie2_index, 此次mm10 
bowtie2_index=/home/User/refdata/mm10
Work_directory=/home/User/CUT_Tag
cd $Work_directory/result/2.trimmed_fq
mkdir $Work_directory/result/3.Align
mkdir $Work_directory/result/3.Align/mapped_fragments_bam
mkdir $Work_directory/result/3.Align/mapped_fragments_bam/bed

for i in *_1P.fq.gz; 
  do
	i=${i%_1P.fq.gz*};  
	(bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 20 -x $bowtie2_index -1 $Work_directory/result/2.trimmed_fq/${i}_1P.fq.gz -2 $Work_directory/result/2.trimmed_fq/${i}_2P.fq.gz) 2> $Work_directory/result/3.Align/${i}.bowtie2 | samtools view -bS - >$Work_directory/result/3.Align/${i}_aligned_reads.bam

#for i in *_1P.fq.gz; 
  #do
	#i=${i%_1P.fq.gz*};  
	# Filter and keep the mapped read pairs 保留匹配的reads 
	#samtools view -bh -f 3 -F 4 -F 8 $Work_directory/result/3.Align/${i}_aligned_reads.bam > $Work_directory/result/3.Align/mapped_fragments_bam/${i}.step1.bam
	
	## Convert into bed file format
	#bedtools bamtobed -i $Work_directory/result/3.Align/mapped_fragments_bam/${i}.step1.bam -bedpe >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.bed
	
	## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
	#awk '$1==$4 && $6-$2 < 1000 {print $0}' $Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.bed >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.clean.bed
	
	## Only extract the fragment related columns
	#cut -f 1,2,6 $Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.fragments.bed

done
```

`--end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700`Bowtie2 使用用于映射长度为 10-700 bp 的插入物的参数对齐配对末端读数。

##### 3.2.2. 与spike-in基因组对齐以进行spike-in校准[可选/推荐]
大肠杆菌 DNA 与细菌产生的 pA-Tn5 蛋白一起携带，并在反应过程中被非特异性标记。映射到大肠杆菌基因组的总读数的比例取决于表位靶向 CUT&Tag 的产量，因此取决于使用的细胞数量和染色质中该表位的丰度。由于将恒定量的 pATn5 添加到 CUT&Tag 反应中并带来固定量的大肠杆菌 DNA，因此大肠杆菌读数可用于在一组实验中标准化表位丰度。
```bash
##== linux command ==##
spikeInRef="/shared/ngs/illumina/henikoff/Bowtie2/Ecoli"
chromSize="/fh/fast/gottardo_r/yezheng_working/SupplementaryData/hg38/chromSize/hg38.chrom.size"

## bowtie2-build path/to/Ecoli/fasta/Ecoli.fa /path/to/bowtie2Index/Ecoli
bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${cores} -x ${spikeInRef} -1 ${projPath}/fastq/${histName}_R1.fastq.gz -2 ${projPath}/fastq/${histName}_R2.fastq.gz -S $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam &> $projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.txt

seqDepthDouble=`samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2_spikeIn.sam | wc -l`
seqDepth=$((seqDepthDouble/2))
echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDepth
```
-   对于spike-in标准化，读数与大肠杆菌基因组U00096.3比对，另外两个参数`--no-overlap`和`--no-dovetail`（`--end-to-end --very-sensitive --no-overlap --no-dovetail --no-mixed --no-discordant --phred33 -I 10 -X 700`校准。

#### 3.3 报告排序映射汇总【必填】
总结原始读取和唯一映射读取以报告对齐效率。对于高质量数据，对齐频率预计将 > 80%。CUT&Tag 数据通常具有非常低的背景，因此只需 100 万个映射片段就可以为人类基因组中的组蛋白修饰提供可靠的配置文件。对不太丰富的转录因子和染色质蛋白的分析可能需要 10 倍于下游分析的映射片段。

我们可以评估以下指标：
-   测序深度
-   对齐率
-   可映射片段数
-   复制率
-   独特的库大小
-   片段大小分布

##### 1. 测序深度
```R
##=== R command ===## 
## Path to the project and histone list
projPath = "
/home/User/CUT_Tag"
sampleList = c("K27me3_rep1", "K27me3_rep2", "K4me3_rep1", "K4me3_rep2", "IgG_rep1", "IgG_rep2")
histList = c("K27me3", "K4me3", "IgG")

## Collect the alignment results from the bowtie2 alignment summary files
alignResult = c()
for(hist in sampleList){
  alignRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, ".bowtie2"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))




##=== R command ===## 
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(viridis)
## Generate sequencing depth boxplot
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (hg38)")

fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (hg38)")

fig3D = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Spike-in Alignment Rate") +
    xlab("") +
    ggtitle("D. Alignment Rate (E.coli)")

ggarrange(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
```

![](Epigenomics/images/Sequencing_depth.jpg)
![](Epigenomics/images/unnamed-chunk-12-1.png)


##### 2. Spike-in alignment
```R
##=== R command ===## 
spikeAlign = c()
for(hist in sampleList){
  spikeRes = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(spikeRes$V1[6], 1, nchar(as.character(spikeRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  spikeAlign = data.frame(Histone = histInfo[1], Replicate = histInfo[2], 
                          SequencingDepth = spikeRes$V1[1] %>% as.character %>% as.numeric, 
                          MappedFragNum_spikeIn = spikeRes$V1[4] %>% as.character %>% as.numeric + spikeRes$V1[5] %>% as.character %>% as.numeric, 
                          AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
}
spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
```

#### 3.4 Filter BAM(Samtools)
#### 3.4.1 去除PCR重复
CUT&Tag 将接头整合到与抗体连接的 pA-Tn5 附近的 DNA 中，整合的确切位点受到周围 DNA 可及性的影响。由于这个原因，共享确切起始和结束位置的片段应该是常见的，并且这种“重复”可能不是由于 PCR 期间的重复。在实践中，我们发现高质量 CUT&Tag 数据集的表观重复率较低，即使是明显的“重复”片段也可能是真正的片段。因此，我们不建议删除重复项。在使用非常少量材料或怀疑 PCR 重复的实验中，可以去除重复。以下命令显示如何使用[Picard](https://broadinstitute.github.io/picard/)检查重复率。

PCR扩增和一些重复序列（如微卫星、着丝粒）会产生重复，干扰真实的富集信号，所以在call peaks前需要先去除重复，这里先用picard去除PCR重复。picard去除PCR重复时要加上参数`REMOVE_DUPLICATES=true`,否则只是标记了duplicates，并没有去除。

<mark>**.bw文件最后放在IGV对结果进行可视化，因为bam文件太大。**</mark>

```bash
module load conda/miniconda3
module load Java
module load GATK/v4.2.2.0

# Set the directory
Output_directory=/home/User/CUT_Tag/result/3.Align
File=/home/User/CUT_Tag/SRR_Acc_List.txt
mkdir $Output_directory/dup_marked
mkdir $Output_directory/rmdup
# Run
cat $File | while read LINE || [[ -n ${LINE} ]]
do
  a=$(echo $LINE |sed 's/\n//g')
  sample=$(echo $a |sed 's/\r//g')
  echo $sample 
  
  #gatk SortSam -I $Output_directory/mapped_fragments_bam/${sample}.step1.bam -O $Output_directory/mapped_fragments_bam/${sample}.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY SILENT
  
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

> **总结了明显的重复率并计算了没有重复的唯一库大小。**

```R
##=== R command ===## 
## Summarize the duplication information from the picard summary outputs.
dupResult = c()
for(hist in sampleList){
  dupRes = read.table(paste0(projPath, "/alignment/removeDuplicate/picard_summary/", hist, "_picard.rmDup.txt"), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}
dupResult$Histone = factor(dupResult$Histone, levels = histList)
alignDupSummary = left_join(alignSummary, dupResult, by = c("Histone", "Replicate", "MappedFragNum_hg38")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))
alignDupSummary
```
![](Epigenomics/images/removeDuplicate.jpg)
-   在这些示例数据集中，IgG 对照样本具有相对较高的重复率，因为该样本中的读数源自 CUT&Tag 反应中的非特异性标记。因此，在下游分析之前从 IgG 数据集中删除重复项是合适的。
    
-   估计的文库大小是基于 Picard 计算的 PE 重复估计的文库中独特分子的数量。
    
-   估计的文库大小与目标表位的丰度和所用抗体的质量成正比，而 IgG 样品的估计文库大小预计会非常低。
    
-   唯一片段数由 MappedFragNum_hg38 * (1-DuplicationRate/100) 计算。

```R
##=== R command ===## 
## generate boxplot figure for the  duplication rate
fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Duplication Rate (*100%)") +
    xlab("") 

fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Estimated Library Size") +
    xlab("") 

fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("# of Unique Fragments") +
    xlab("")

ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
```
![](Epigenomics/images/unnamed-chunk-15-1.png)

>**实际应用例子（仅供参考）**
```R
#---------------------
###Remove duplicates
#---------------------
setwd("/home/User/CUT_Tag/result/3.Align/rmdup")
##=== R command ===## 
## Summarize the duplication information from the picard summary outputs.
all_sampleList <- read.table("../samplelist.txt") 
#--------------------------------------------
########## H3K9ac analysis  #################
#--------------------------------------------
sampleList1 <- all_sampleList[grep("H3K9ac",all_sampleList$V1),]

histList1 = rep(c("H3K9ac-Dp16","H3K9ac-Dp16-VPA","H3K9ac-WT"),each=4)

dupResult = c()
for(hist in sampleList1){
  dupRes = read.table(paste0("/home/User/CUT_Tag/result/3.Align/rmdup.", hist, "_picard.rmDup.txt"), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Group = histInfo[1], Replicate = histInfo[1], MappedFragNum_mm10 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric,
  						DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, 
  						EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_mm10 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}

#dupResult$Histone = factor(dupResult$Histone, levels = histList)
dupResult$Group = histList1
dupResult$name <- dupResult$Replicate
dupResult$Replicate <- rep(c("rep1","rep2","rep3","rep4"),3)
dupResult <- dupResult[,c(1,7,2:6)]
write.table(dupResult,"H3K9ac_sample_duplicated_rate.txt",quote=F,sep="\t",col.names=T,row.names=F)

##=== R command ===## 
## generate boxplot figure for the  duplication rate
dupResult$Group <- gsub("H3K9ac-","",dupResult$Group)

fig4A = dupResult %>% ggplot(aes(x = Group, y = DuplicationRate, fill = Group)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Duplication Rate (*100%)") +
    xlab("") 

fig4B = dupResult %>% ggplot(aes(x = Group, y = EstimatedLibrarySize, fill = Group)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Estimated Library Size") +
    xlab("") 

fig4C = dupResult %>% ggplot(aes(x = Group, y = UniqueFragNum, fill = Group)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("# of Unique Fragments") +
    xlab("")

p1 <- ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
ggsave("H3K9ac_duplicated_rate_fig4ABC.png", p1, width =16 , height = 8)

```

#### 3.4.2 评估映射片段大小分布[required]
CUT&Tag 在连接酶附近的染色质颗粒两侧插入接头，尽管染色质颗粒内的标记也可能发生。因此，针对组蛋白修饰的 CUT&Tag 反应主要产生核小体长度 (~180 bp) 或该长度的倍数的片段。靶向转录因子的 CUT&Tag 主要分别从相邻的核小体和因子结合位点产生核小体大小的片段和可变数量的较短片段。核小体表面的 DNA 标记也会发生，用单碱基对分辨率绘制片段长度显示 10 bp 的锯齿周期，这是成功的 CUT&Tag 实验的典型特征。

```bash
##== linux command ==##

mkdir -p $projPath/alignment/sam/fragmentLen

## Extract the 9th column from the alignment sam file which is the fragment length
samtools view -F 0x04 $projPath/alignment/sam/${histName}_bowtie2.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >$projPath/alignment/sam/fragmentLen/${histName}_fragmentLen.txt
```

```R
##=== R command ===## 
## Collect the fragment size information
fragLen = c()
for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, "/alignment/sam/fragmentLen/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[1], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)
## Generate the fragment size density plot (violin plot)
fig5A = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(bw = 5) +
    scale_y_continuous(breaks = seq(0, 800, 50)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ggpubr::rotate_x_text(angle = 20) +
    ylab("Fragment Length") +
    xlab("")

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Histone, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

ggarrange(fig5A, fig5B, ncol = 2)
```
![](Epigenomics/images/unnamed-chunk-17-1.png)


#### 3.4.3 评估实验重复性
通过对整个基因组的映射读取计数进行相关性分析来评估重复之间的数据再现性。

##### 3.4.3.1 文件格式转换
```bash
# Run
module load conda/miniconda3
module load Java
module load bowtie/v2.4.4
# 去bowtie2的官网下载相关物种版本的bowtie2_index, 此次mm10 
Work_directory=/home/User/CUT_Tag
cd $Work_directory/result/2.trimmed_fq
mkdir $Work_directory/result/3.Align
mkdir $Work_directory/result/3.Align/mapped_fragments_bam
mkdir $Work_directory/result/3.Align/mapped_fragments_bam/bed

for i in *_1P.fq.gz; 
  do
	i=${i%_1P.fq.gz*};  
	# Filter and keep the mapped read pairs 保留匹配的reads 
	samtools view -bh -f 3 -F 4 -F 8 $Work_directory/result/3.Align/${i}_aligned_reads.bam > $Work_directory/result/3.Align/mapped_fragments_bam/${i}.step1.bam
	
	## Convert into bed file format
	bedtools bamtobed -i $Work_directory/result/3.Align/mapped_fragments_bam/${i}.step1.bam -bedpe >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.bed
	
	## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.bed >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.clean.bed
	
	## Only extract the fragment related columns
	cut -f 1,2,6 $Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.fragments.bed

done
```

##### 3.4.3.2. 评估重复实验的重现性 Assess replicate reproducibility

为了研究重复之间和跨条件的可重复性，将基因组分成 500 bp 的 bin，并在重复数据集之间计算<mark>每个 bin 中读取计数的 log2 转换值的 Pearson 相关性。</mark>多个重复和 IgG 控制数据集显示在分层聚类相关矩阵中。
>我们使用每个片段的中点来推断该片段属于哪个 500bp 的 bin。
```bash
##== linux command ==##
## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.

cd $Work_directory/result/3.Align/mapped_fragments_bam/bed

for i in *_bowtie2.fragments.bed; 
do
i=${i%_bowtie2.fragments.bed*}; 
echo ${i}

binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >$Work_directory/result/3.Align/mapped_fragments_bam/bed/${i}_bowtie2.fragmentsCount.bin$binLen.bed

done
```
###### 1.  使用R语言进行评估重复重现性
```R
# directory
setwd("/home/User/CUT_Tag/result/3.Align/mapped_fragments_bam/bed")
bam_dir <- "/home/User/CUT_Tag/result/3.Align/mapped_fragments_bam"
sampleList <- unique(strsplit(list.files(bam_dir,pattern = "*.bam"),".bam")[[1]])

##== R command ==##
library(corrplot)

reprod = c()
fragCount = NULL
for(sample in sampleList){
  print(sample)
  sampleInfo = strsplit(sample, "_")[[1]]
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0("",sample, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", sampleInfo[1])
  
  }else{
    
    fragCountTmp = read.table(paste0("",sample, "_bowtie2.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", sampleInfo[1])
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}
head(fragCount)
write.table(fragCount,"sample_bowtie2_fragmentsCount_bin500.txt",sep="\t",col.names=T,row.names=F,quote=F)

# pearson
M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs",method= "pearson") 
png("sample_duplicated_corrplot_pearson.png")
pic <- corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()

# spearman
M1= cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs",method= "spearman") 
png("sample_duplicated_corrplot_spearman.png")
pic <- corrplot(M1, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))
dev.off()
```
###### 2. deeptools 相关性分析
deeptools的相关介绍见**ATAC-seq pipeline.md**。
```bash
# install
conda install -c bioconda deeptools
# Command line installation without pip
git clone https://github.com/deeptools/deepTools.git
cd deepTools
## --prefix 根据真实路径更改
python setup.py install --prefix /User/Tools/deepTools2.0
```
____
```bash
#---------------------------------------------------------------
## Deeptools to calculate the correlation of duplicated samples
#---------------------------------------------------------------
cd deepTools2.0
module load homor/v4.11

# Set the directory
Work_directory=/home/User/CUT_Tag
mkdir $Work_directory/deepTools
mkdir $Work_directory/deepTools/correlation

#1 calculat the coverage of bam file
# according to the data
multiBamSummary bins \
--bamfiles *****.bam *****.bam \
--binSize 10000 \
--numberOfProcessors 10 \
--labels WT1 WT2 WT3 PS19_1 PS19_2 PS19_3 PS19_Usp25_1 PS19_Usp25_2 PS19_Usp25_3 \
--outRawCounts $Work_directory/deepTools/correlation/results1.txt \
-o $Work_directory/deepTools/correlation/results1.npz 

#2 visualization
cd $Work_directory/deepTools/correlation

## spearman 
plotCorrelation \
 -in reads.npz \
 --corMethod spearman \
 --skipZeros \
 --plotTitle "Sperman Correlation of Read Counts" \
 --whatToPlot heatmap \
 --colorMap RdYlBu \
 --plotNumbers \
 -o heatmap_SpearmanCorr.png \
 --outFileCorMatrix SpearmanCorr_readCounts.tab

## pearson
plotCorrelation \
 -in results.npz \
 --corMethod pearson \
 --skipZeros \
 --plotTitle "Pearson Correlation of Average Scores Per Transcript" \
 --whatToPlot scatterplot \
 -o scatterplot_PearsonCorr.png \
 --outFileCorMatrix PearsonCorr_bigwigScores.tab
```
![](Epigenomics/images/unnamed-chunk-21-1.png)

______
#### 3.5  Spike-in calibration[optional]
此部分是**可选**的，但根据您的实验方案**推荐。**在第 3.2.1 节中展示了与spike-in 基因组的比对，在第3.3节中展示了spike-in 比对摘要。
基本假设是定位到初级基因组的片段与大肠杆菌基因组的比率对于一系列样本是相同的，每个样本使用相同数量的细胞。由于这个假设，我们不会在实验之间或纯化的 pATn5 批次之间进行标准化，这可能具有非常不同数量的携带大肠杆菌 DNA。使用常数 C 来避免标准化数据中的小部分，我们将比例因子 S 定义为
`S = C / (fragments mapped to E. coli genome)`
然后将归一化的覆盖率计算为：
`Normalized coverage = (primary_genome_coverage) * S`
常数是一个任意乘数，通常为 10,000。作为基因组覆盖图文件，生成的文件将相对较小。
```bash
##== linux command ==##
if [[ "$seqDepth" -gt "1" ]]; then
    
    mkdir -p $projPath/alignment/bedgraph

    scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph
    
fi
```
##### 3.5.1  Scaling factor
```R
##=== R command ===## 
scaleFactor = c()
multiplier = 10000
for(hist in sampleList){
  spikeDepth = read.table(paste0(projPath, "/alignment/sam/bowtie2_summary/", hist, "_bowtie2_spikeIn.seqDepth"), header = FALSE, fill = TRUE)$V1[1]
  
  histInfo = strsplit(hist, "_")[[1]]
  scaleFactor = data.frame(scaleFactor = multiplier/spikeDepth, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(scaleFactor, .)
}
scaleFactor$Histone = factor(scaleFactor$Histone, levels = histList)
left_join(alignDupSummary, scaleFactor, by = c("Histone", "Replicate"))
```
![](Epigenomics/images/Scalingfactor.jpg)
```R
##=== R command ===##
## Generate sequencing depth boxplot
fig6A = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ylab("Spike-in Scalling Factor") +
    xlab("")

normDepth = inner_join(scaleFactor, alignResult, by = c("Histone", "Replicate")) %>% mutate(normDepth = MappedFragNum_hg38 * scaleFactor)

fig6B = normDepth %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ylab("Normalization Fragment Count") +
    xlab("") + 
    coord_cartesian(ylim = c(1000000, 130000000))
ggarrange(fig6A, fig6B, ncol = 2, common.legend = TRUE, legend="bottom")
```
![](Epigenomics/images/unnamed-chunk-24-1.png)
______

#### 3.6 Peak Calling 
##### 3.6.1  MACS2
peaks calling 有不同的方法，MACS2是最常用的call peaks工具。 [MACS全称Model-based Analysis of ChIP-Seq](https://github.com/taoliu/MACS)，最初的设计是用来鉴定转录因子的结合位点，但是它也可以用于其他类型的富集方式测序。  

```bash
module load conda/miniconda3
source activate macs2
# Set the directory
Input_directory=/home/User/CUT_Tag/result/3.Align
Output_directory=/home/User/CUT_Tag/result
mkdir $Output_directory/peak
mkdir $Output_directory/peak/filter_q0.05
mkdir $Output_directory/peak/filter_p0.05
# Run
cd $Input_directory
for i in *.bam; 
do
  i=${i%.bam*}; 
  # q_value=0.05
  mkdir $Output_directory/peak/filter_q0.05/${i}
  macs2 callpeak -t $Input_directory/${i}.bam -g mm -f BAMPE -n ${i} --outdir $Output_directory/peak/filter_q0.05/${i} -q 0.05 -B --SPMR --keep-dup all 2> $Output_directory/peak/filter_q0.05/${i}/${i}.macs2
  ### add control "-c"
  # p_value=0.05
  mkdir $Output_directory/peak/filter_p0.05/${i}
  macs2 callpeak -t $Input_directory/${i}.bam -g mm -f BAMPE -n ${i} --outdir $Output_directory/peak/filter_p0.05/${i} -p 0.05 -B --SPMR --keep-dup all 2> $Output_directory/peak/filter_p0.05/${i}/${i}.macs2
done

# Finished
conda deactivate
```

##### 3.6.2 SEACR
###### 3.6.2.1 使用SEACR进行peak calling
```bash
##== linux command ==##
seacr="/fh/fast/gottardo_r/yezheng_working/Software/SEACR/SEACR_1.3.sh"
histControl=$2
mkdir -p $projPath/peakCalling/SEACR

bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph \
     $projPath/alignment/bedgraph/${histControl}_bowtie2.fragments.normalized.bedgraph \
     non stringent $projPath/peakCalling/SEACR/${histName}_seacr_control.peaks

bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks
```

###### 3.6.2.2  Number of peaks called
```R
##=== R command ===## 
peakN = c()
peakWidth = c()
peakType = c("control", "top0.01")
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[1], Replicate = histInfo[2]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(peakWidth, .)
    }
  }
}
peakN %>% select(Histone, Replicate, peakType, peakN)
```
![](Epigenomics/images/Number_of_peaks_called.jpg)

###### 3.6.2.3 生物重复峰的重现性 Reproducibility of the peak across biological replicates
比较重复数据集上的峰值调用以定义可重复的峰值。前 1% 的峰（按每个块中的总信号排序）被选为高置信度位点。

```R
##=== R command ===## 
histL = c("K27me3", "K4me3")
repL = paste0("rep", 1:2)
peakType = c("control", "top0.01")
peakOverlap = c()
for(type in peakType){
  for(hist in histL){
    overlap.gr = GRanges()
    for(rep in repL){
      peakInfo = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_", type, ".peaks.stringent.bed"), header = FALSE, fill = TRUE)
      peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      }else{
        overlap.gr = peakInfo.gr
        
      }
    }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Histone = hist, peakType = type) %>% rbind(peakOverlap, .)
  }
}

peakReprod = left_join(peakN, peakOverlap, by = c("Histone", "peakType")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)
peakReprod %>% select(Histone, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)
```

![](Epigenomics/images/Reproducibilityofthepeakacrossbiologicalreplicates.jpg)

重现性由下式计算
``# peaks overlapping rep1 and rep2/# peaks of rep1 or rep2 * 100`
因此，它对每次重复调用的峰总数很敏感。

###### 3.6.2.4 峰区的片段比例（FRiPs）
计算峰中读取的分数 (FRiP) 作为信噪比的度量，并将其与 IgG 控制数据集中的 FRiP 进行对比以进行说明。尽管 CUT&Tag 的测序深度通常只有 1-5 百万次读取，但该方法的低背景会导致 FRiP 得分较高。
```R
##=== R command ===## 
library(chromVAR)

bamDir = paste0(projPath, "/alignment/bam")
inPeakData = c()
## overlap with bam file to get count
for(hist in histL){
  for(rep in repL){
    peakRes = read.table(paste0(projPath, "/peakCalling/SEACR/", hist, "_", rep, "_seacr_control.peaks.stringent.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0(bamDir, "/", hist, "_", rep, "_bowtie2.mapped.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Histone = hist, Replicate = rep))
  }
}

frip = left_join(inPeakData, alignResult, by = c("Histone", "Replicate")) %>% mutate(frip = inPeakN/MappedFragNum_hg38 * 100)
frip %>% select(Histone, Replicate, SequencingDepth, MappedFragNum_hg38, AlignmentRate_hg38, FragInPeakNum = inPeakN, FRiPs = frip)
```

![](Epigenomics/images/FRiPs.jpg)

###### 3.6.2.5 峰数、峰宽、峰重现性和 FRiP 的可视化
```R
fig7A = peakN %>% ggplot(aes(x = Histone, y = peakN, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    facet_grid(~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Number of Peaks") +
    xlab("")

fig7B = peakWidth %>% ggplot(aes(x = Histone, y = width, fill = Histone)) +
    geom_violin() +
    facet_grid(Replicate~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
    theme_bw(base_size = 18) +
    ylab("Width of Peaks") +
    xlab("")

fig7C = peakReprod %>% ggplot(aes(x = Histone, y = peakReprodRate, fill = Histone, label = round(peakReprodRate, 2))) +
    geom_bar(stat = "identity") +
    geom_text(vjust = 0.1) +
    facet_grid(Replicate~peakType) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Peaks Reproduced") +
    xlab("")

fig7D = frip %>% ggplot(aes(x = Histone, y = frip, fill = Histone, label = round(frip, 2))) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Fragments in Peaks") +
    xlab("")

ggarrange(fig7A, fig7B, fig7C, fig7D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
```

![](Epigenomics/images/unnamed-chunk-29-1.png)

#### 3.7. 可视化
##### 3.7.1 IGV
通常，我们有兴趣使用基因组浏览器可视化区域中的染色质景观。[Integrative Genomic Viewer](http://software.broadinstitute.org/software/igv/home)提供易于使用的 Web 应用程序版本和本地桌面版本。[UCSC 基因组浏览器](https://genome.ucsc.edu/)提供最全面的补充基因组信息。
![](Epigenomics/images/chr7.png)

##### 3.7.2 特定区域的热图可视化
我们还对查看注释位点列表中的染色质特征感兴趣，例如基因启动子处的组蛋白修饰信号。我们将使用[deepTools](https://deeptools.readthedocs.io/en/develop/)`computeMatrix`中的和`plotHeatmap`函数来生成热图。

>**参考基因组的GTF文件的下载**
```bash
# 1. 小鼠的参考基因组
mkdir RefGenome
cd /home/User/RefGenome

wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gtf.gz

gzip -dk GCF_000001635.26_GRCm38.p6_genomic.gtf.gz

# 2. 人类的参考基因组
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz

gzip -dk GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
```

###### 3.7.2.1 转录单位的热图
```bash
##== linux command ==##
cores=8
computeMatrix scale-regions -S $projPath/alignment/bigwig/K27me3_rep1_raw.bw \
                               $projPath/alignment/bigwig/K27me3_rep2_raw.bw \
                               $projPath/alignment/bigwig/K4me3_rep1_raw.bw \
                               $projPath/alignment/bigwig/K4me3_rep2_raw.bw \
                              -R $projPath/data/hg38_gene/hg38_gene.tsv \
                              --beforeRegionStartLength 3000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 3000 \
                              --skipZeros -o $projPath/data/hg38_gene/matrix_gene.mat.gz -p $cores

plotHeatmap -m $projPath/data/hg38_gene/matrix_gene.mat.gz -out $projPath/data/hg38_gene/Histone_gene.png --sortUsing sum
```

![](Epigenomics/images/Histone_gene.png)

>实际例子参考
```bash
# for many samples
module load homor/v4.11
# Set the directory
Work_directory=/home/User/CUT_Tag
mkdir $Work_directory/deepTools/heatmap
cd $Work_directory/3.Align/rmdup

# e.g. H3K9ac
mv H3K9ac*.bw plot/H3K9ac/

computeMatrix reference-point --referencePoint TSS \
				  -S $Work_directory/Align/rmdup/*.bw \
				  -R /home/User/RefGenome/Mus_musculus.GRCm38.94.gtf \
				  -a 3000 \
				  -b 3000 \
				  --skipZeros \
				  -p 10 \
				  -o $Work_directory/deepTools/H3K9ac_sample_TSS.mat_4duplication.gz
  
  plotHeatmap -m $Work_directory/deepTools/H3K9ac_sample_TSS.mat_4duplication.gz \
  -out $Work_directory/deepTools/heatmap/HEATMAP_H3K9ac_sample_TSS.mat_4duplication.pdf \
  --colorMap 'RdBu' \
  --whatToShow 'plot, heatmap and colorbar' \
  --samplesLabel WT1 WT2 WT3 WT4 Dp16_rep1 Dp16_rep2 Dp16_rep3 Dp16_rep4 Dp16-VPA_rep1 Dp16-VPA_rep2 Dp16-VPA_rep3 Dp16-VPA_rep4 \
  --heatmapWidth 6 \
  --dpi 720
done

#--colorMap 'RdBu' 'RdBu' 'RdBu' 'jet' 'jet' 'jet' 'jet' 'jet' 'jet' \
```
###### 3.7.2.2 CUT&Tag 峰的热图
我们使用从 SEACR 返回的信号块的中点来对齐热图中的信号。SEACR 输出的第六列是 chr:start-end 形式的条目，表示该区域具有最大信号的区域的第一个和结束碱基。我们首先在第 6 列中生成一个包含此中点信息的新床文件，并使用 deeptools 进行热图可视化。

```bash
##== linux command ==##
awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $projPath/peakCalling/SEACR/${histName}_${repName}_seacr_control.pe\
aks.stringent.bed >$projPath/peakCalling/SEACR/${histName}_${repName}_seacr_control.peaks.summitRegion.bed

computeMatrix reference-point -S $projPath/alignment/bigwig/${histName}_${repName}_raw.bw \
              -R $projPath/peakCalling/SEACR/${histName}_${repName}_seacr_control.peaks.summitRegion.bed \
              --skipZeros -o $projPath/peakCalling/SEACR/${histName}_${repName}_SEACR.mat.gz -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $projPath/peakCalling/SEACR/${histName}_SEACR.mat.gz -out $projPath/peakCalling/SEACR/${histName}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${histName} ${repName}"
```

![](Epigenomics/images/Histone_peak.png)



___
>可视化两类样本中reads的分布特征
```bash
#### Genome-wide comparison of KO/WT
#png
plotProfile -m ./H3K9ac_sample_TSS.mat_4duplication.gz \
            -out ./H3K9ac_sample_TSS.mat_4duplication_TSS_binding4.pdf \
            --samplesLabel WT1 WT2 WT3 WT4 Dp16_rep1 Dp16_rep2 Dp16_rep3 Dp16_rep4 Dp16-VPA_rep1 Dp16-VPA_rep2 Dp16-VPA_rep3 Dp16-VPA_rep4 \
            --legendLocation upper-left \
            --plotHeight 10 \
            --plotWidth 10 \
            --yMin 2 \
            --perGroup \
            --yAxisLabel "CUT&Tag signal" 
 
plotProfile -m ./H3K9ac_sample_TSS.mat_4duplication.gz \
            -out ./H3K9ac_sample_TSS.mat_4duplication_TSS_binding4.png \
            --samplesLabel WT1 WT2 WT3 WT4 Dp16_rep1 Dp16_rep2 Dp16_rep3 Dp16_rep4 Dp16-VPA_rep1 Dp16-VPA_rep2 Dp16-VPA_rep3 Dp16-VPA_rep4 \
            --legendLocation upper-left \
            --plotHeight 12 \
            --plotWidth 12 \
            --yMin 2 \
            --perGroup \
            --yAxisLabel "CUT&Tag signal" 
```

![](Epigenomics/images/plotprofile.png)


#### 3.8 Peak Annotation
```bash
# linux
# prepare for annotation
Work_directory=/home/User/CUT_Tag/result4.peak/filter_q0.05
mkdir $Work_directory/annotation && mkdir $Work_directory/annotation/H3K9ac
system('ls ../../*/*narrowPeak|grep "H3K9ac" >H3K9ac_samplelist.txt')
```

##### 3.8.1  ChIPseeker 进行峰注释
>**该部分的内容可以见ATAC-seq pipeline.md，以下的代码是实际例子的应用**

```R
# annotation
library(ChIPseeker)
library(GenomicFeatures)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(clusterProfiler)

setwd("/home/User/CUT_Tag/result4.peak/filter_q0.05/annotation/H3K9ac")
gtf='/home/User/RefGenome/mm10.refGene.gtf'
sample <- read.table("H3K9ac_samplelist.txt")
peaks <- sample$V1[c(9:12,1:8)]

outname=c("WT_1","WT_2","WT_3","WT_4","Dp16_1","Dp16_2","Dp16_3","Dp16_4","Dp16_VPA_1","Dp16_VPA_2","Dp16_VPA_3","Dp16_VPA_4")

tx=makeTxDbFromGFF(file = gtf)
head(seqlevels(tx))
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
head(seqlevels(txdb))
peak=list()
for(i in 1:length(peaks)){
  peak[[i]]=readPeakFile(peakfile=peaks[i], header=F, as = 'GRanges')
}
head(seqlevels(peak[[1]]))

peakAnno=list()
for(i in 1:length(peak)){
  peakAnno[[i]]=annotatePeak(peak = peak[[i]], tssRegion=c(-3000,3000), TxDb = txdb, annoDb='org.Mm.eg.db',assignGenomicAnnotation=T,flankDistance=5000)
}
#peakAnno[[1]]=as.data.frame(peakAnno[[1]])
#peakAnno=annotatePeak(peak = peak, tssRegion=c(-3000, 3000),TxDb = tx, annoDb='org.Mm.eg.db', assignGenomicAnnotation = T)
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

# for(i in 1:length(peakAnno)){
# 	pdf(file=paste0(outname[i],".peakAnnotation.pdf"))
# 	plotAnnoPie(peakAnno[[i]], main=paste0(outname[i],"\nDistribution of Peaks"), line=-8)
# 	dev.off()
# 	#write.table(as.data.frame(peakAnno[[i]]), file = paste0(outname[i], ".peakAnno.txt"), sep='\t', row.names=F,quote=F)
# }

peakAnno1 <- peakAnno
# 以下部分可以用loop赋值
names(peakAnno1)[1] <- "WT_1"
names(peakAnno1)[2] <- "WT_2"
names(peakAnno1)[3] <- "WT_3"
names(peakAnno1)[4] <- "WT_4"
names(peakAnno1)[5] <- "Dp16_1"
names(peakAnno1)[6] <- "Dp16_2"
names(peakAnno1)[7] <- "Dp16_3"
names(peakAnno1)[8] <- "Dp16_4"
names(peakAnno1)[9] <- "Dp16_VPA_1"
names(peakAnno1)[10] <- "Dp16_VPA_2"
names(peakAnno1)[11] <- "Dp16_VPA_3"
names(peakAnno1)[12] <- "Dp16_VPA_4"

pdf("allsample_Bar_peakAnnotation.pdf")
plotAnnoBar(peakAnno1, title="Feature Distribution")
dev.off()
save(peakAnno1, file="12allsample_Bar_peakAnnotation.rda")

#select samples
peakAnno <- peakAnno1[c(3,4,7,8,11,12)]
 names(peakAnno)
[1] "WT_3"       "WT_4"       "Dp16_3"     "Dp16_4"     "Dp16_VPA_3"
[6] "Dp16_VPA_4"
pdf("6samples_Bar_peakAnnotation.pdf")
plotAnnoBar(peakAnno, title="Feature Distribution")
dev.off()
```

##### 3.8.2  超几何检验
```bash
# linux
# prepare for annotation
Work_directory=/home/User/CUT_Tag/result4.peak/filter_q0.05/annotation/H3K9ac
mkdir $Work_directory/Fisher_Exact_Test && cd $Work_directory/Fisher_Exact_Test
```

```R
setwd("/home/User/CUT_Tag/result4.peak/filter_q0.05/annotation/H3K9ac/Fisher_Exact_Test")
load("../12allsample_Bar_peakAnnotation.rda")
peakAnno_df=list()
for(i in 1:length(peakAnno1)){
  peakAnno_df[[i]]=as.data.frame(peakAnno1[[i]])
}

getGenomicAnnoStat1 <- function(peakAnno) {
    anno <- peakAnno$annotation
    e1 <- getOption("ChIPseeker.ignore_1st_exon")
    i1 <- getOption("ChIPseeker.ignore_1st_intron")
    ids <- getOption("ChIPseeker.ignore_downstream")

    if (is.null(e1) || !e1) {
        e1lab <- "1st Exon"
        anno[grep("exon 1 of", anno)] <- e1lab
        exonlab <- "Other Exon"
    } else {
        e1lab <- NULL
        exonlab <- "Exon"
    }

    if (is.null(i1) || !i1) {
        i1lab <- "1st Intron"
        anno[grep("intron 1 of", anno)] <- i1lab
        intronlab <- "Other Intron"
    } else {
        i1lab <- NULL
        intronlab <- "Intron"
    }

    anno[grep("Exon \\(", anno)] <- exonlab
    anno[grep("Intron \\(", anno)] <- intronlab
    if (is.null(ids) || !ids) {
        dsd <- getOption("ChIPseeker.downstreamDistance")
        if (is.null(dsd))
            dsd <- 3000 ## downstream 3k by default
        if (dsd > 1000) {
            dsd <- round(dsd/1000, 1)
            dsd <- paste0(dsd, "kb")
        }
        dslab <- paste0("Downstream (<=", dsd, ")")

        anno[grep("Downstream", anno)] <- dslab
        iglab <- "Distal Intergenic"
    } else {
        dslab <- NULL
        iglab <- "Intergenic"
        anno[grep("Downstream", anno)] <- iglab
    }
    anno[grep("^Distal", anno)] <- iglab
    return(anno)
}

peakAnno_df1 <- lapply(peakAnno_df,function(x){
    anno1 <- getGenomicAnnoStat1(x)
    x$annotation1 <- anno1
    return(x)
})


anno_table <- lapply(peakAnno_df1,function(x){
    anno <- data.frame(table(x$annotation1))
    colnames(anno) <- c("name","Freq")
    rownames(anno) <- anno$name
    anno$name <- NULL
    return(anno)
})

anno_sum <- Reduce(cbind,anno_table)
colnames(anno_sum) <- c("WT_1","WT_2","WT_3","WT_4","Dp16_1","Dp16_2","Dp16_3","Dp16_4","Dp16_VPA_1","Dp16_VPA_2","Dp16_VPA_3","Dp16_VPA_4") 

anno_sum$WT <- apply(anno_sum[,1:4],1,mean)
anno_sum$Dp16 <- apply(anno_sum[,5:8],1,mean)
anno_sum$Dp16VPA <- apply(anno_sum[,9:12],1,mean)
anno_sum1 <- apply(anno_sum[,13:15],2,round)
anno_sum1 <- data.frame(anno_sum1)
anno_sum1["sum",] <- apply(anno_sum1,2,sum)

#save.data
anno_sum2 <- anno_sum
anno_sum2 <- apply(anno_sum2,2,round)
anno_sum2 <- data.frame(anno_sum2)
anno_sum2$name <- rownames(anno_sum2)
anno_sum2 <- anno_sum2[,c(ncol(anno_sum2),1:ncol(anno_sum2)-1)]
write.table(anno_sum2,"allsample_peak_annotation.txt",sep="\t",col.names=T,row.names=F,quote=F)

```


```R
### 利用超几何分布检验分析各分组 consensus peak 注释图

library(data.table)

# anno_sum <- fread("../allsample_peak_annotation.txt",header=T)
# anno_sum <- data.frame(anno_sum)
# rownames(anno_sum) <- anno_sum$name
# anno_sum$name <- NULL
# anno_sum1 <- anno_sum[,c(2,3,5:7,9:11)]
# anno_sum1$WT <- apply(anno_sum1[,1:2],1,mean)
# anno_sum1$PS19 <- apply(anno_sum1[,3:5],1,mean)
# anno_sum1$PS19Usp25 <- apply(anno_sum1[,6:8],1,mean)
# anno_sum1 <- apply(anno_sum1,2,round)
# anno_sum1 <- data.frame(anno_sum1)
# anno_sum2 <- anno_sum1[,c("WT","PS19","PS19Usp25")]
# anno_sum2["sum",] <- apply(anno_sum2,2,sum)

#DP16 vs WT
Dp16_wt <- anno_sum1[,c("WT","Dp16")]
hyper_Dp16_wt <- NULL
for(i in 1:(nrow(Dp16_wt)-1)){
    print(i)
    hyper1 <- matrix(c(Dp16_wt[i,2],Dp16_wt[i,1],Dp16_wt[12,2]-Dp16_wt[i,2],Dp16_wt[12,1]-Dp16_wt[i,1]),
    nrow=2,ncol=2,byrow = TRUE,dimnames = list(c(rownames(Dp16_wt)[i],paste0("not ",rownames(Dp16_wt)[i])),c("Dp16", "wt")))
    print(hyper1[1,1])
    print(hyper1)
    n = hyper1[1,2] + hyper1[2,2]
    k = hyper1[1,1] + hyper1[1,2]
    m = hyper1[1,1] + hyper1[2,1]
    q = hyper1[1,1]
    p1 <- data.frame(phyper(q,m,n,k, lower.tail=FALSE))
    colnames(p1) <- rownames(Dp16_wt)[i]
    rownames(p1) <- "pvalue"

    a = hyper1[1,1]
    b = hyper1[1,2] 
    c = hyper1[2,1]
    d = hyper1[2,2]
    or <- data.frame((a*d)/(b*c))

    colnames(or) <- rownames(Dp16_wt)[i]
    rownames(or) <- "OR"

    k1 <- data.frame(t(p1),t(or))
    hyper_Dp16_wt <- rbind(hyper_Dp16_wt,k1)
}

hyper_Dp16_wt
                        pvalue        OR
1st Exon           0.027768516 1.1441040
1st Intron         0.010470543 1.0803023
3' UTR             0.034695059 1.1571822
5' UTR             0.159318468 1.1815359
Distal Intergenic  0.004673889 1.0681540
Downstream (<=300) 0.222368292 1.2517477
Other Exon         0.046338148 1.1121964
Other Intron       0.997998480 0.9256381
Promoter (<=1kb)   0.985356061 0.9603045
Promoter (1-2kb)   0.925344386 0.9274562
Promoter (2-3kb)   0.364088969 1.0195296

hyper_Dp16_wt$name <- rownames(hyper_Dp16_wt)
hyper_Dp16_wt <- hyper_Dp16_wt[,c("name","pvalue","OR")]

#DP16 vs Dp16VPA
Dp16_DP16VPA <- anno_sum1[,c("Dp16VPA","Dp16")]

hyper_Dp16_DP16VPA <- NULL
for(i in 1:(nrow(Dp16_DP16VPA)-1)){
    print(i)
    hyper2 <- matrix(c(Dp16_DP16VPA[i,2],Dp16_DP16VPA[i,1],Dp16_DP16VPA[12,2]-Dp16_DP16VPA[i,2],Dp16_DP16VPA[12,1]-Dp16_DP16VPA[i,1]),
    nrow=2,ncol=2,byrow = TRUE,dimnames = list(c(rownames(Dp16_DP16VPA)[i],paste0("not ",rownames(Dp16_DP16VPA)[i])),c("usp25", "ps19")))
    print(hyper2[1,1])
    print(hyper2)
    n = hyper2[1,2] + hyper2[2,2]
    k = hyper2[1,1] + hyper2[1,2]
    m = hyper2[1,1] + hyper2[2,1]
    q = hyper2[1,1]
    p1 <- data.frame(phyper(q,m,n,k, lower.tail=FALSE))
    colnames(p1) <- rownames(Dp16_DP16VPA)[i]
    rownames(p1) <- "pvalue"

    a = hyper2[1,1]
    b = hyper2[1,2] 
    c = hyper2[2,1]
    d = hyper2[2,2]
    or <- data.frame((a*d)/(b*c))

    colnames(or) <- rownames(Dp16_DP16VPA)[i]
    rownames(or) <- "OR"

    k1 <- data.frame(t(p1),t(or))
    hyper_Dp16_DP16VPA <- rbind(hyper_Dp16_DP16VPA,k1)
}
hyper_Dp16_DP16VPA
                         pvalue        OR
1st Exon           1.871038e-02 1.1501716
1st Intron         9.307318e-01 0.9541985
3' UTR             6.663518e-01 0.9655426
5' UTR             2.145765e-01 1.1296007
Distal Intergenic  9.801398e-01 0.9518718
Downstream (<=300) 5.308163e-01 0.9176808
Other Exon         6.350119e-01 0.9778947
Other Intron       1.000000e+00 0.8051396
Promoter (<=1kb)   7.040535e-22 1.1834381
Promoter (1-2kb)   9.759325e-01 0.9070602
Promoter (2-3kb)   9.726001e-01 0.8960744

hyper_Dp16_DP16VPA$name <- rownames(hyper_Dp16_DP16VPA)
hyper_Dp16_DP16VPA <- hyper_Dp16_DP16VPA[,c("name","pvalue","OR")]

#PS19Usp25_WT
Dp16VPA_wt <- anno_sum1[,c("WT","Dp16VPA")]

hyper_Dp16VPA_wt <- NULL
for(i in 1:(nrow(Dp16VPA_wt)-1)){
    print(i)
    hyper3 <- matrix(c(Dp16VPA_wt[i,2],Dp16VPA_wt[i,1],Dp16VPA_wt[12,2]-Dp16VPA_wt[i,2],Dp16VPA_wt[12,1]-Dp16VPA_wt[i,1]),
    nrow=2,ncol=2,byrow = TRUE,dimnames = list(c(rownames(Dp16VPA_wt)[i],paste0("not ",rownames(Dp16VPA_wt)[i])),c("Dp16VPA", "wt")))
    print(hyper3[1,1])
    print(hyper3)
    n = hyper3[1,2] + hyper3[2,2]
    k = hyper3[1,1] + hyper3[1,2]
    m = hyper3[1,1] + hyper3[2,1]
    q = hyper3[1,1]
    p1 <- data.frame(phyper(q,m,n,k, lower.tail=FALSE))
    colnames(p1) <- rownames(Dp16VPA_wt)[i]
    rownames(p1) <- "pvalue"

    a = hyper3[1,1]
    b = hyper3[1,2] 
    c = hyper3[2,1]
    d = hyper3[2,2]
    or <- data.frame((a*d)/(b*c))

    colnames(or) <- rownames(Dp16VPA_wt)[i]
    rownames(or) <- "OR"

    k1 <- data.frame(t(p1),t(or))
    hyper_Dp16VPA_wt <- rbind(hyper_Dp16VPA_wt,k1)
}

hyper_Dp16VPA_wt
                         pvalue        OR
1st Exon           5.144904e-01 0.9947246
1st Intron         8.422278e-05 1.1321568
3' UTR             1.142939e-02 1.1984787
5' UTR             3.706171e-01 1.0459766
Distal Intergenic  2.072905e-06 1.1221616
Downstream (<=300) 1.585078e-01 1.3640339
Other Exon         2.005518e-02 1.1373376
Other Intron       2.274044e-08 1.1496615
Promoter (<=1kb)   1.000000e+00 0.8114531
Promoter (1-2kb)   3.190569e-01 1.0224858
Promoter (2-3kb)   1.320795e-02 1.1377734

hyper_Dp16VPA_wt$name <- rownames(hyper_Dp16VPA_wt)
hyper_Dp16VPA_wt <- hyper_Dp16VPA_wt[,c("name","pvalue","OR")]

hyper_Dp16_wt$group <- "Dp16/WT"
hyper_Dp16VPA_wt$group <- "Dp16VPA/WT"
hyper_Dp16_DP16VPA$group <- "Dp16/Dp16VPA"

OR_ratio <- rbind(hyper_Dp16_wt,hyper_Dp16VPA_wt,hyper_Dp16_DP16VPA)

write.table(OR_ratio,"H3k9ac_Dp16_wt_Dp16VPA_peak_annotation_hyper_result.txt",sep="\t",quote=F,col.names=T,row.names=F)
save(OR_ratio,file="H3k9ac_Dp16_wt_Dp16VPA_peak_annotation_hyper_result.rda")

library(data.table)
library(ggsci)
library(ggplot2)
png("H3k9ac_Dp16_wt_Dp16VPA_annotation_OR.png",width = 600, height = 480)
OR_ratio$name <- factor(OR_ratio$name,levels=c(
"Promoter (2-3kb)","Promoter (1-2kb)","Promoter (<=1kb)","5' UTR","3' UTR","1st Exon","Other Exon","1st Intron","Other Intron","Downstream (<=300)","Distal Intergenic"))

p <- ggplot(OR_ratio, aes(name,OR,fill=factor(group)),stat = "identity")+
            geom_bar(stat="identity",position="dodge")+coord_flip()+
            geom_hline(aes(yintercept=1),colour="dimgray",linetype="solid")+ #dimgray red
            xlab("")+ 
            ylab("Fold Enrichment")+
            labs(fill = "Group")+
            ylim(0,5)+
            theme_bw()+
            scale_fill_npg()+
            #scale_fill_tron()+
            #scale_fill_lancet()+
            theme( axis.title=element_text(size=8,face="bold"),
                    #axis.ticks.length=unit(0.2,'cm'),
                    axis.text.x=element_text(colour="black",size=10,face="bold"), 
                    axis.text.y=element_text(colour="black",size=10,face="bold"), 
                    axis.title.y=element_text(size = 10,face="bold"), 
                    panel.border = element_blank(),axis.line = element_line(colour = "black",size=0.5), 
                    legend.text=element_text(colour="black",size=10,face="bold"),  
                    legend.title=element_text(colour="black", size=10,face="bold"),
                    #legend.position = c(0.9,0.85),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    plot.title = element_text(hjust = 0.5))

print(p)
dev.off()

```

#### 3.9 差异peaks分析
##### 3.9.1 DESeq2
```R
dir.create("/home/User/CUT_Tag/result/4.peak/filter_q0.05/DEG_narrowpeak")
dir.create("/home/User/CUT_Tag/result/4.peak/filter_q0.05/DEG_narrowpeak/H3K9ac")
setwd("/home/User/CUT_Tag/result/4.peak/filter_q0.05/DEG_narrowpeak/H3K9ac")

library(dplyr)
library(stringr)
library(tidyverse)
library(ggplot2)
library(GenomicRanges)
library(chromVAR) ## For FRiP analysis and differential analysis(Fraction of reads in peaks)
library(DESeq2) ## For differential analysis section
library(ggpubr) ## For customizing figures
library(clusterProfiler)

### Peak结果注释
system('ls ../../*/*narrowPeak|grep "H3K9ac" >H3K9ac_narrowpeak_samplelist.txt')

sample <- read.table("H3K9ac_samplelist.txt")
histL<- sample$V1[c(9:12,1:8)]

mPeak = GRanges()
## overlap with bam file to get count
for(hist in histL){
  #peakRes = read.table(paste0("",peaksDir,"",hist,"/",hist,"_peaks.narrowPeak"), header = FALSE, fill = TRUE)
  peakRes = read.table(hist, header = FALSE, fill = TRUE)
  mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}
masterPeak = reduce(mPeak)
##=== Annotation ===##
library(ChIPseeker)
library(GenomicFeatures)
require(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(rtracklayer)
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
head(seqlevels(txdb))
peakAnno=annotatePeak(peak = masterPeak, tssRegion=c(-3000,3000), TxDb = txdb, annoDb='org.Mm.eg.db',assignGenomicAnnotation=T, flankDistance=5000)  
peakAnno=as.data.frame(peakAnno)
index = as.data.frame(masterPeak)

## From bamfile 
system('ls /home/User/CUT_Tag/result/3.Align/rmdup/*bam|grep "H3K9ac" >H3K9ac_bam_samplelist.txt')
sample1 <- read.table("H3K9ac_bam_samplelist.txt")
histL2<- sample1$V1[c(9:12,1:8)]

countMat = matrix(NA, length(masterPeak), length(histL2))
## overlap with bam file to get count

i = 1
for(hist in histL2){
  	#bamFile = paste0(bamDir, "/", hist)
  	bamFile = hist
    fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
    countMat[, i] = counts(fragment_counts)[,1]
    i = i + 1
}

colnames(countMat) = histL2
colnames(countMat) = gsub('.bam','',colnames(countMat))
colnames(countMat) = gsub('/home/User/CUT_Tag/result/3.Align/rmdup/','',colnames(countMat))
colnames(countMat) = stringr::str_split_fixed(colnames(countMat),"_",2)[,1]

dif = setdiff(index$start,peakAnno$start)
subscript = c()
for (i in 1:length(dif)) {
  subscript[i] = which(index$start == dif[i])
}

countMat=countMat[-subscript,]

rownames(countMat) = paste(peakAnno$seqnames, peakAnno$start, peakAnno$end, peakAnno$SYMBOL,sep = "_")
#rownames(countMat) = paste0('chr', rownames(countMat))
#countMat = countMat[!grepl('^chrMT',rownames(countMat)),]
countMat = countMat[!grepl('^chrMT|chr1_GL456211_random|chr1_GL456221_random|chr4_GL456216_random|chr4_GL456350_random|chr4_JH584293_random|chr4_JH584294_random|chr5_GL456354_random|chrUn_JH584304|chrX_GL456233_random',rownames(countMat)),]

#colnames(countMat) <- c("WT_1","WT_2","WT_3","WT_4","Dp16_1","Dp16_2","Dp16_3","Dp16_4","Dp16_VPA_1","Dp16_VPA_2","Dp16_VPA_3","Dp16_VPA_4")
#countMat <- countMat[,c(1,4:5,11,6:8,12,9,2:3,10)]
##=== R command ===## 
save(countMat,file="H3K9ac_12sample_countMat_narrowpeaks.rda")


### 以下同RNA-seq的DESeq2
#==============
###  Dp16/WT     
#==============
load("H3K9ac_12sample_countMat_narrowpeaks.rda")

countMat_Dp16 <- countMat[,c(1:8)] #four duplicated

selectR = which(rowSums(countMat_Dp16) > 3) ## remove low count genes
dataS = countMat_Dp16[selectR,]
condition = factor(rep(c('WT','KO'), times=c(4,4)), levels = c('WT',"KO"))

dds = DESeqDataSetFromMatrix(countData = dataS,
                              colData = DataFrame(condition),
                              design = ~ condition)
DDS = DESeq(dds)
normDDS = counts(DDS, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS) = paste0(colnames(normDDS), "_norm")
res = results(DDS, independentFiltering = FALSE, altHypothesis = "greaterAbs")
summary(res$pvalue)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.00000 0.03513 0.21474 0.31534 0.55231 1.00000
countMatDiff = data.frame(cbind(dataS, normDDS, res),check.names=F)

countMatDiff <- countMatDiff[!grepl("^chrMT|chr1_GL|chr4_JH|chr4_GL|chr5_GL|chrUn_JH|chrX_GL|_random_",rownames(countMatDiff)),]
countMatDiff$narrow_peak <- rownames(countMatDiff)
countMatDiff <- separate(countMatDiff,narrow_peak,c("chr","peak_start","peak_end","gene"),sep="_")
countMatDiff <- countMatDiff[,c((ncol(countMatDiff)-3):ncol(countMatDiff),1:(ncol(countMatDiff)-4))]
dim(countMatDiff)
#[1] 70809     26

save(countMatDiff,file="H3K9ac_Dp16_WT_countMatDiff_narrowPeaks.rda")
#write.table(countMatDiff, file="PS19_WT_countMatDiff_narrowPeaks.txt", sep = "\t", quote=F,col.names=T,row.names=F)

sig_countMatDiff <- countMatDiff[countMatDiff$pvalue<0.05,]
#20116
write.table(sig_countMatDiff, file="H3K9ac_Dp16_WT_countMatDiff_narrowPeaks_significant_p0.05.txt", sep = "\t", quote=F,col.names=T,row.names=F)

sig_logFC1_countMatDiff <- sig_countMatDiff[sig_countMatDiff$log2FoldChange>=1 | sig_countMatDiff$log2FoldChange<= -1,]
#3060
write.table(sig_logFC1_countMatDiff, file="H3K9ac_Dp16_WT_countMatDiff_narrowPeaks_significant_p0.05_logFC_1.txt", sep = "\t", quote=F,col.names=T,row.names=F)

sig_countMatDiff1 <- countMatDiff[countMatDiff$padj<0.05,]
#10428
write.table(sig_countMatDiff1, file="H3K9ac_Dp16_WT_countMatDiff_narrowPeaks_significant_q0.05.txt", sep = "\t", quote=F,col.names=T,row.names=F)

sig_logFC1_countMatDiff1 <- sig_countMatDiff1[sig_countMatDiff1$log2FoldChange>=1 | sig_countMatDiff1$log2FoldChange<= -1,]
#1733(1327 uniq genes)
write.table(sig_logFC1_countMatDiff1, file="H3K9ac_Dp16_WT_countMatDiff_narrowPeaks_significant_q0.05_logFC_1.txt", sep = "\t", quote=F,col.names=T,row.names=F)

```

##### 3.9.2 差异peaks分析--DiffBind
DiffBind是鉴定两个样本间差异结合位点的一个R包。主要用于peak数据集，包括对peaks的重叠和合并的处理，计算peaks重复间隔的测序reads数，并基于结合亲和力鉴定具有统计显著性的差异结合位点。适用的统计模型有DESeq、DESeq2、edgeR。详细内容可参考DiffBind的文档：[DiffBind](http://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf)

```R
# install the packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DiffBind")
# differential analysis
library(DiffBind)
dbObj <- dba(sampleSheet="SampleSheet.csv")
## 计算每个peaks/regions的count信息
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
dba.plotPCA(dbObj, attributes=DBA_FACTOR, label=DBA_ID)
plot(dbObj)
# Establishing a contrast
dbObj <- dba.contrast(dbObj, categories=DBA_FACTOR,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# summary of results
dba.show(dbObj, bContrasts=T)
# overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
# Save the Result
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
## EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="results/Nanog_vs_Pou5f1_edgeR.txt", sep="t", quote=F, col.names = NA)
## DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="results/Nanog_vs_Pou5f1_deseq2.txt", sep="t", quote=F, col.names = NA)
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05),c("seqnames", "start", "end", "strand", "Fold")]
write.table(edge.bed, file="results/Nanog_vs_Pou5f1_edgeR_sig.bed", sep="t", quote=F, row.names=F, col.names=F)
# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05),
c("seqnames", "start", "end", "strand", "Fold")]
write.table(deseq.bed, file="results/Nanog_vs_Pou5f1_deseq2_sig.bed", sep="t", quote=F, row.names=F, col.names=F)

```


### 参考文章：
[CUTTag_tutorial](https://yezhengstat.github.io/CUTTag_tutorial/index.html)
