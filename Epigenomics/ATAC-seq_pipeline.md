# ATAC-seq pipeline

<mark>ATAC-seq和Chip-seq分析流程基本一致，具体个性化分析可以参考CUT_Tag_pipeline的部分内容。</mark>
<!-- TOC -->

- [ATAC-seq pipeline](#atac-seq-pipeline)
  - [（一）技术背景](#一技术背景)
  - [（二）与Chip-Seq的异同](#二与chip-seq的异同)
  - [（三）开放染色质的研究方法](#三开放染色质的研究方法)
  - [（四）分析流程](#四分析流程)
    - [1. Quality control 质量控制](#1-quality-control-质量控制)
      - [1.1 FastQC](#11-fastqc)
      - [1.2 Filter data 数据过滤(Trimmomatic)](#12-filter-data-数据过滤trimmomatic)
    - [2. Alignment to Genome 比对至基因组](#2-alignment-to-genome-比对至基因组)
      - [2.1 建立索引](#21-建立索引)
      - [2.2 比对](#22-比对)
    - [3. Filter BAM(Samtools)](#3-filter-bamsamtools)
    - [4. Peak Calling (MACS2)](#4-peak-calling-macs2)
    - [5.  Quality assessment of ChIP 质量评估[Optional]](#5--quality-assessment-of-chip-质量评估optional)
      - [（1）phantompeakqualtools](#1phantompeakqualtools)
      - [（2）ChIPQC](#2chipqc)
    - [6. Combine Peak Calls  [Optional]](#6-combine-peak-calls--optional)
      - [(1) Overlapping peaks using bedtools](#1-overlapping-peaks-using-bedtools)
      - [(2) Irreproducibility Discovery Rate (IDR)重复样本的处理-IDR](#2-irreproducibility-discovery-rate-idr重复样本的处理-idr)
    - [7. Annotation 对峰进行注释](#7-annotation-对峰进行注释)
      - [（1）Homer](#1homer)
      - [（2）ChIPseeker](#2chipseeker)
    - [8. 差异peaks分析——DiffBind](#8-差异peaks分析diffbind)

<!-- /TOC -->

## （一）技术背景
哺乳动物的DNA通过三个主要的层次尺度进行高度压缩：

* a. 第一层是核小体，DNA缠绕在核小体单体上；

* b. 第二层是染色质，核小体单体互相缠绕形成染色质；

* c. 第三层是染色体，染色质进一步压缩为染色体。

DNA压缩的三个尺度及其相互作用共同造就了基因表达调控。

而DNA的复制转录是需要将DNA的紧密结构打开，从而允许一些调控因子结合（转录因子或其他调控因子）。
这部分打开的染色质，就叫<mark>开放染色质，打开的染色质允许其他调控因子结合的特性称为染色质的可及性（chromatin accessibility）</mark>。因此，认为染色质的可及性与转录调控密切相关。  

![](Epigenomics/images/ATAC-seq1.png)

上图比较直观地展示了ATAC-seq的技术原理，其中NFR fragments是开放染色质中两个核小体之间的Linker DNA片段，由Peak Calling分析得到的Peak来鉴定；蓝色的Footprint则是转录因子的足迹，对应的是转录因子足迹分析；核小体单体（Mononucleosome)上结合的DNA片段则反映出核小体的位置信息，对应的是核小体占位分析。

开放染色质的研究方法有ATAC-seq以及传统的DNase-Seq及FAIRE-seq等，ATAC-Seq由于所需细胞量少，实验简单，可以在全基因组范围内检测染色质的开放状态，目前已经成为研究染色质开放性的首选技术方法。

## （二）与Chip-Seq的异同
ATAC-Seq与ChIP-Seq的不同的是：
1. ATAC-Seq是全基因组范围内检测染色质的开放程度，可以得到全基因组范围内的蛋白质可能结合的位点信息，一般用于不知道特定的转录因子，用此方法与其他方法结合筛查感兴趣的特定调控因子；
2. 但是ChIP-Seq是明确知道感兴趣的转录因子是什么，根据感兴趣的转录因子设计抗体去做ChIP实验拉DNA，验证感兴趣的转录因子是否与DNA存在相互作用。
## （三）开放染色质的研究方法
**ATAC-Seq、ChIP-Seq、Dnase-Seq、MNase-Seq、FAIRE-Seq整体的分析思路一致，找到富集区域，对富集区域进行功能分析。**
-   **ChIP-Seq**是揭示特定转录因子或蛋白复合物的结合区域，实际是研究DNA和蛋白质的相互作用，利用抗体将蛋白质和DNA一起富集，并对富集到的DNA进行测序。
-   **DNase-Seq、ATAC-Seq、FAIRE-Seq**都是用来研究开放染色质区域。
	- DNase-Seq是用的DNase I内切酶识别开放染色质区域；
	- ATAC-seq是用的Tn5转座酶，随后进行富集和扩增；
	- FAIRE-Seq是先进行超声裂解，然后用酚-氯仿富集
-   **MNase-Seq**是用来鉴定核小体区域

## （四）分析流程
![](Epigenomics/images/ATAC-seq2.png)
### 1. Quality control 质量控制
#### 1.1 FastQC
```bash
module load conda/miniconda3
module load Java

# Set the directory
Input_data=/home/User/ATAC-seq/cleandata/
Output_directory=/home/User/ATAC-seq
mkdir $Work_directory/result
mkdir $Output_directory/result/1.QC

# QC  
find $Input_data -name *fq.gz| xargs fastqc -t 20 -o $Output_directory/result/1.QC/
```

#### 1.2 Filter data 数据过滤(Trimmomatic)
```bash
module load conda/miniconda3
module load Java

# Set the directory
Input_data=/home/User/cleandata/
Output_directory=/home/User/ATAC-seq
mkdir $Output_directory/result/2.Trim

# trimmomatic clean data  
cd $Input_data

for i in *_1.clean.fq.gz; 
do
	i=${i%_1.clean.fq.gz*}; 
	trimmomatic PE -baseout $Output_directory/result/2.Trim/${i}.fq.gz \
	${i}_1.clean.fq.gz ${i}_2.clean.fq.gz \
	ILLUMINACLIP:/share/old_apps/Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:8:true \
	LEADING:3 \
	TRAILING:3 \
	SLIDINGWINDOW:4:15 \
	MINLEN:8
done
```

### 2. Alignment to Genome 比对至基因组
#### 2.1 建立索引
[bowtie2 官网](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
<mark>**推荐使用！**</mark>
![](Epigenomics/images/bowtie2_install_ref.jpg)
bowtie2官网（mm10) :
```bash
# install
conda install bowtie2
pkg install bowtie2
cd refdata
unzip -o /home/User/refdata/mm10 mm10.zip
rm mm10.zip make_mm10.sh

# index directory
#/home/User/refdata/mm10
```
NCBI官网：

1. 人类参考基因组
```bash
wget -c ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz
gunzip GCF_000001635.26_GRCm38.p6_genomic.fna.gz
```
2. 小鼠参考基因组 
```bash
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

```bash
# Homo sapiens
bowtie2-build -f /home/User/refdata/GCF_000001635.26_GRCm38.p6_genomic.fna --threads 24 GRCh38
# Mus musculus
bowtie2-build -f /home/User/refdata/GCF_000001405.40_GRCh38.p14_genomic.fna --threads 24 mm10
```

#### 2.2 比对
Bowtie2是一个快速精确的比对工具，基于Burrows-Wheeler Transform 构建基因组的FM 索引，比对过程所耗内存少。Bowtie2支持局部、双端、缺口比对模式，对大于50bp的reads比对效果更好（小于50bp的reads用Bowtie1）。

>比对这一步骤中，ATAC-seq使用的bowtie2参数为其默认参数，或者可以加入参数 `-X 2000` ,允许长达2 Kb的片段。

该步骤将得到sam文件，由于sam文件所占空间较大，再运行过程中放入管道转换为bam文件。
```bash
module load conda/miniconda3
module load Java
module load bowtie/v2.4.4

# Set the directory
bowtie2_index='/home/User/refdata/mm10'
Output_directory=/home/User/ATAC-seq/result
mkdir $Output_directory/3.Align
mkdir $Output_directory/3.Align/mapped_fragments_bam

# Alignment
# Alignment
cd $Output_directory/2.Trim
for i in *_1P.fq.gz; 
do
  i=${i%_1P.fq.gz*}; 
  (bowtie2 --no-mixed -p 12 -x $bowtie2_index -1 $Output_directory/2.Trim/${i}_1P.fq.gz -2 $Output_directory/2.Trim/${i}_2U.fq.gz) 2> $Output_directory/3.Align/${i}.bowtie2 | samtools view -bS - >$Output_directory/3.Align/${i}_aligned_reads.bam  
  
  # Filter and keep the mapped read pairs
  samtools view -bh -f 3 -F 4 -F 8 $Output_directory/3.Align/${i}_aligned_reads.bam > $Output_directory/3.Align/mapped_fragments_bam/${i}.step1.bam
  
  ## Convert into bed file format
  bedtools bamtobed -i $Output_directory/3.Align/mapped_fragments_bam/${i}.step1.bam -bedpe >$Output_directory/3.Align/mapped_fragments_bam/${i}_bowtie2.bed
  
  ## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
  awk '$1==$4 && $6-$2 < 1000 {print $0}' $Output_directory/3.Align/mapped_fragments_bam/${i}_bowtie2.bed > $Output_directory/1.Align/mapped_fragments_bam/${i}_bowtie2.clean.bed
  
  ## Only extract the fragment related columns
  cut -f 1,2,6 $Output_directory/3.Align/mapped_fragments_bam/${i}_bowtie2.clean.bed | sort -k1,1 -k2,2n -k3,3n  > $Output_directory/3.Align/mapped_fragments_bam/${i}_bowtie2.fragments.bed

done
```


### 3. Filter BAM(Samtools)
**去除PCR重复**  
PCR扩增和一些重复序列（如微卫星、着丝粒）会产生重复，干扰真实的富集信号，所以在call peaks前需要先去除重复，这里先用picard去除PCR重复。picard去除PCR重复时要加上参数`REMOVE_DUPLICATES=true`,否则只是标记了duplicates，并没有去除。
```bash
module load conda/miniconda3
module load Java
module load GATK/v4.2.2.0

# Set the directory
Output_directory=/home/User/ATAC-seq/result/3.Align
## 所有样本名
File=/home/User/ATAC-seq/SRR_Acc_List.txt
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

### 4. Peak Calling (MACS2)
peaks calling 有不同的方法，MACS2是最常用的call peaks工具。 [MACS全称Model-based Analysis of ChIP-Seq](https://github.com/taoliu/MACS)，最初的设计是用来鉴定转录因子的结合位点，但是它也可以用于其他类型的富集方式测序。  
MACS通过整合序列标签位置信息和方向信息提高结合位点的空间分辨率。MACS的工作流如下所示：

![](Epigenomics/images/ATAC-seq4.png)

**参数解析：**（参考：[callpeak](https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md)）
1. **输入文件选项**
-   `-t`：IP 数据文件（这是 MACS 的唯一必需参数）
-   `-c`：控制或模拟数据文件
-   `-f`: 输入文件的格式；默认为“AUTO”，这将允许 MACS 自动决定格式。
-   `-g`：可映射的基因组大小，定义为可以测序的基因组大小；提供了一些预编译的值。
2. **输出参数**
-   `--outdir`: MACS2 会将所有输出文件保存到此选项的指定文件夹中
-   `-n`: 输出文件的前缀字符串
-   `-B/--bdg`：将片段堆积、控制 lambda、-log10pvalue 和 -log10qvalue 分数存储在 bedGraph 文件中
3. **Shifting model arguments 转移模型参数**
-   `-s`：测序标签的大小。默认情况下，MACS 将使用输入处理文件中的前 10 个序列来确定它
-   `--bw`：仅用于模型构建的基因组扫描带宽。可以设置为预期的超声处理片段大小。
-   `--mfold`：模型构建的上限和下限
4. **Peak calling arguments 峰值调用参数**
-   `-q`: q 值（最小 FDR）截止
-   `-p`: p 值截止（而不是 q 值截止）
-   `--nolambda`：不考虑峰值候选区域的局部偏差/lambda
-   `--broad`: 广泛的峰值呼叫

**<mark>ATAC-seq关心的是在哪切断，断点才是peak的中心，所以使用shift模型，–shift -75或-100。</mark>** 

```bash
module load conda/miniconda3
source activate macs2

# Set the directory
Input_directory=/home/User/ATAC-seq/3.Align/rmdup
Output_directory=/home/User/ATAC-seq/result
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
  macs2 callpeak -t $Input_directory/${i}.bam -g mm -f BAMPE -n ${i} --outdir $Output_directory/4.peak/filter_q0.05/${i} -q 0.05 -B --SPMR --nomodel --shift -75 --extsize 250 --keep-dup all 2> $Output_directory/4.peak/filter_q0.05/${i}/${i}.macs2
  ### add control "-c"
  # p_value=0.05
  mkdir $Output_directory/4.peak/filter_p0.05/${i}
  macs2 callpeak -t $Input_directory/${i}.bam -g mm -f BAMPE -n ${i} --outdir $Output_directory/4.peak/filter_p0.05/${i} -p 0.05 -B --SPMR --nomodel --shift -75 --extsize 250 --keep-dup all 2> $Output_directory/4.peak/filter_p0.05/${i}/${i}.macs2
  
done

# Finished
conda deactivate
```

* **输出文件：**
	-   `_peaks.narrowPeak`: BED6+4 格式文件，包含峰值位置以及峰值、pvalue 和 qvalue
	-   `_peaks.xls`：一个表格文件，其中包含有关被调用峰的信息。附加信息包括堆积和折叠富集
	-   `_summits.bed`：每个峰的峰顶位置。要在结合位点查找基序，建议使用此文件
	-   `_model.R`：一个 R 脚本，可用于根据您的数据和互相关图生成有关模型的 PDF 图像
	-   `_control_lambda.bdg`：输入样本的bedGraph格式
	-   `_treat_pileup.bdg`：处理样本的bedGraph格式
_____
**MACS输出结果解读：**
* 1.  `NAME_peaks.xls`是一个表格文件，其中包含有关称为峰的信息。您可以在 excel 中打开它并使用 excel 函数进行排序/过滤。信息包括：
  
    -   染色体名称
    -   峰值起始位置
    -   峰的结束位置
    -   峰区长度
    -   绝对峰顶位置
    -   峰顶堆积高度
    -   -log10(pvalue) 用于峰顶（例如 pvalue =1e-10，那么这个值应该是 10）
    -   针对局部 lambda 的随机泊松分布，该峰峰的倍数富集，
    -   -log10(qvalue) 在峰顶
  
    XLS 中的坐标是基于 1 的，这与 BED 格式不同。当`--broad`启用宽峰调用时，XLS 文件中的堆积、p 值、q 值和倍数变化将是整个峰区域的平均值，因为在宽峰调用模式下不会调用峰顶.
  
2.  `NAME_peaks.narrowPeak`是 BED6+4 格式文件，其中包含峰值位置以及峰值、p 值和 q 值。您可以将其加载到 UCSC 基因组浏览器。一些特定列的定义是：
    
    -   5th：整数分数显示。它的计算方式 `int(-10*log10pvalue)`或`int(-10*log10qvalue)`取决于`-p`(pvalue) 或`-q`(qvalue) 是否用作分数截止值。请注意，当前此值可能超出[UCSC ENCODE 窄峰格式](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)中定义的 [0-1000] 范围。您可以使用以下 1-liner awk 让值在 1000 处饱和（即 p/q-value = 10^-100）：`awk -v OFS="\t" '{$5=$5>1000?1000:$5} {print}' NAME_peaks.narrowPeak`
    -   7th：峰顶的倍数变化
    -   8th: -log10pvalue at peak Summit
    -   9th：-log10qvalue 在峰顶
    -   10th：相对峰顶位置到峰顶开始
    
    该文件可以直接加载到 UCSC 基因组浏览器。如果您想通过其他工具对其进行分析，请删除起始轨迹线。
    
3.  `NAME_summits.bed`是 BED 格式，其中包含每个峰的峰顶位置。此文件中的第 5 列与文件中的内容相同`narrowPeak`。如果您想在绑定位置找到图案，建议使用此文件。该文件可以直接加载到 UCSC 基因组浏览器。如果您想通过其他工具对其进行分析，请删除起始轨迹线。
    
4.  `NAME_peaks.broadPeak`是 BED6+3 格式，类似于 `narrowPeak`文件，除了缺少注释峰顶的第 10 列。此文件和`gappedPeak`文件只有在启用时才可用`--broad`。由于在宽峰调用模式下，不会调用峰顶，因此第 5 列和第 7-9 列的值是峰区所有位置的平均值。`narrowPeak`如果要修复第 5 列中的值问题，请参阅。
    
5.  `NAME_peaks.gappedPeak`采用 BED12+3 格式，包含宽峰和窄峰。第 5 列是在 UCSC 浏览器上显示灰度的分数，如`narrowPeak`. 第 7 列是该区域第一个窄峰的起点，第 8 列是终点。第 9 列应该是 RGB 颜色键，但是，我们在此处保留 0 以使用默认颜色，因此如果需要，请更改它。第 10 列表示有多少块，包括广泛区域的起始 1bp 和结束 1bp。第 11 列显示每个块的长度，第 12 列显示每个块的开始。第 13 位：倍数变化，第 14 位：_-log10pvalue_，第 15 位：_-log10qvalue_。该文件可以直接加载到 UCSC 基因组浏览器。`narrowPeak`如果要修复第 5 列中的值问题，请参阅 。
    
6.  `NAME_model.r`是一个 R 脚本，可用于根据数据生成模型的 PDF 图像。通过以下方式将其加载到 R：
    
    `$ Rscript NAME_model.r`
    
    `NAME_model.pdf`然后将在您的当前目录中生成一个pdf文件。请注意，绘制此图需要 R。
    
    <mark>注意！如果结果文件没有R脚本，原因是在 BAMPE 模式下，无需估计片段长度，因为将考虑每个读取对的实际插入长度。这就是为什么没有生成model.R。（参考：[issue#281](https://github.com/macs3-project/MACS/issues/281)）</mark>
    
7.  和文件是 bedGraph 格式`NAME_treat_pileup.bdg`，`NAME_control_lambda.bdg`可以导入到 UCSC 基因组浏览器或转换成更小的 bigWig 文件。包含来自 ChIP/处理样品的堆积信号（根据选项标准化 `NAME_treat_pielup.bdg`） 。包含从对照样本或当对照样本不存在时从处理样本估计的每个基因组位置的局部偏差`--scale-to`。 `NAME_control_lambda.bdg`该子命令`bdgcmp`可用于比较这两个文件并制作一个包含 p 值、q 值、对数似然和对数倍数变化等分数的 bedGraph 文件。

### 5.  Quality assessment of ChIP 质量评估[Optional]
#### （1）[phantompeakqualtools](https://www.plob.org/article/24658.html)
<mark>然而，该包的使用过程中spp的R包无法装上。</mark>
* **链交叉相关（Strand cross-correlation）**
链交叉相关是一个有效的评估ChIP-Seq质量的方法，它不依赖于peak calling，而是基于ChIP-Seq实验。如果ChIP-Seq实验成功，DNA富集序列标签（蛋白质相互作用的序列）会在reads的双峰富集中产生显著的聚集。 

>**产生reads的双峰富集的原因如下：**  
在ChIP-Seq实验中，DNA被片段化，蛋白质结合的片段会被免疫沉淀，所以产生了有蛋白质结合的DNA片段（fragments ）。  
DNA的正链从5’端开始被测序（如下图红色reads），DNA负链也从5’末端被测序产生如下图所示的蓝色reads。

![](Epigenomics/images/ATAC-seq5.png)

由于从DNA片段的5′末端测序，使链reads的富集（下图中的蓝色部分）与负链reads的富集（下图红色部分）有少量的相互抵消区域。我们**需要确定峰位移多少碱基数目可以在两个峰间产生最大的相关性。** 我们可以用**交叉相关的度量值（cross-correlation metric）** 计算产生最大相关的位移。

![](Epigenomics/images/ATACvsChIP.jpeg)

```bash
wget https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/phantompeakqualtools/ccQualityControl.v.1.1.tar.gz
tar -xzf ccQualityControl.v.1.1.tar.gz
cd phantompeakqualtools
# 查看README
less README.txt
# Linux下安装
module load R
R
install.packages("caTools", lib="~/R/library")


mkdir -p logs qual
bam_dir=/home/User/ATAC-seq/Align/rmdup
for bam in $bam_dir/A4_FKDL220091111-1a.bam
  do
  bam2=`basename $bam .final.bam`
  Rscript run_spp_nodups.R -c=$bam -savp -out=qual/${bam2}.qual > logs/${bam2}.Rout
done

```

#### （2）[ChIPQC](https://www.plob.org/article/24666.html)
ChIPQC是一个Bioconductor包，输入文件包括BAM和peak文件，可以自动计算一些质量评估值，并产生质量报告。
```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ChIPQC")

library(ChIPQC)

sample = ChIPQCsample("chip.bam")
sample
ChIPQCreport(sample)
```



### 6. Combine Peak Calls  [Optional]
ATAC-Seq要求必须有2次或更多次生物学重复（十分珍贵或者稀有样本除外，但必须做至少2次技术重复）。理论上重复样本的peaks应该有高度的一致性，实际情况并不完全与预期一致。如何评价重复样本的重复性的好坏？如何得到一致性的peaks?
#### (1) Overlapping peaks using bedtools
用`bedtools`计算peaks的overlaps。  
用法：`bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>`
-   `-a`: 参数后加重复样本1（A）
-   `-b`：参数后加重复样本2（B），也可以加多个样本


![](Epigenomics/images/ATAC-seq7.png)


```bash
bedtools intersect -a macs2/~.narrowPeak -b macs2/~.narrowPeak -wo > bedtools/Nanog-overlaps.bed
```

#### (2) Irreproducibility Discovery Rate (IDR)[重复样本的处理-IDR](https://www.plob.org/article/24676.html)
IDR是通过比较一对经过排序的regions/peaks 的列表，然后计算反映其重复性的值。  
IDR在[ENCODE](https://www.encodeproject.org/atac-seq/)和modENCODE项目中被广泛使用，也是[ChIP-seq指南和标准](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3431496/)中的一部分。

**使用IDR的注意事项：**
-   建议使用IDR时，MACS2 call peaks的步骤参数设置不要过于严格，以便鉴定出更多的peaks。
-   使用IDR需要先对MACS2的结果文件narrowPeak根据`-log10(p-value)`进行排序。
```bash
# Call peaks
macs2 callpeak -t sample.final.bam -n sample --shift -100 --extsize 200 --nomodel -B --SPMR -g hs --outdir Macs2_out 2> sample.macs2.log
#Sort peak by -log10(p-value)
sort -k8,8nr NAME_OF_INPUT_peaks.narrowPeak > macs/NAME_FOR_OUPUT_peaks.narrowPeak
# After sort
idr --samples sample_Rep1_sorted_peaks.narrowPeak sample_Rep2_sorted_peaks.narrowPeak --input-file-type narrowPeak --rank p.value --output-file sample-idr --plot --log-output-file sample.idr.log
```

>结果文件详细见[IDR](https://github.com/nboley/idr#output-file-format)



### 7. Annotation 对峰进行注释
**其中，功能注释见GSEA和Transcriptomics的Scripts中GO_KEGG.R。**
#### （1）Homer
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
* **进行注释**
```bash
# a sample
## mm10 为参考基因组
annotatePeaks.pl ./peak/~_peaks.narrowPeak mm10 >~_summits_peak_anno.txt
```

```bash
# for many samples
# Set the directory
Work_directory=/home/User/ATAC-seq
mkdir $Work_directory/Homer

cd $Work_directory/Align/rmdup
for i in *.bam; 
do
  i=${i%.bam*}; 
  annotatePeaks.pl $Work_directory/peak/filter_p0.05/${i}/${i}_peaks.narrowPeak mm10 >$Work_directory/Homer/${i}_summits_peak_anno.txt
done
```


>**快速使用 Deeptools 对 ChIP-seq 数据画图！**

**Deeptools 的用途：**	
- 处理 bam 文件 或者 bam 转化的 bigwig 文件；
- 数据质量控制；
-  作图，比如热图、折线图；
-  其他。

____
转录起始位点（Transcription Start Site, TSS）：是指**一个基因**的5'端转录的第一个碱基，它是与**新生RNA链**第一个核苷酸**相对应DNA链上的碱基**，通常为一个嘌呤（A或G）。

**通常**把**转录起始位点**前即5’末端的序列称为**上游**，而把其后即3‘末端的序列称为**下游**。
**知识点：转录起始位点是指一个基因的5'端转录的第一个碱基。**
____
```bash
# for many samples
module load homor/v4.11
# Set the directory
Work_directory=/home/User/ATAC-seq
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
  --samplesLabel ${i}  \
  --heatmapWidth 6 \
  --dpi 720
done
```
___
得到可视化和数据分析的峰矩阵（需要时再补充相应的路径和文件名）：
```bash
# Generate peak matrix for visualization and statistics
cat A_summits.bed|awk 'BGEIN{OFS="\t"}{print $1,$2-250,$3+250}' > A_summits_250bp.bed
cat A_summits_250bp.bed B_summits_250bp.bed... | awk '{if($2<0){print $1 "\t" 0 "\t" $3} else {print $0}}' |sortBed -i - |bedtools merge -i - > atlas
```


#### （2）ChIPseeker
[ChIPseeker](https://www.plob.org/article/24683.html)虽然最初是为了ChIP-seq注释而写的一个R包，但它不只局限于ChIP-seq，也可用于ATAC-Seq等其他富集peaks注释，也可用于lincRNA注释、DNA breakpoints的断点位置注释等所有genomic coordination的注释，另外提供了丰富的可视化功能。

**使用方法：**

使用ChIPseeker需要准备两个文件：
* 要注释的peaks的文件，需满足BED格式。
* 注释参考文件，即需要一个包含注释信息的TxDb对象

* **从UCSC下载参考基因信息做注释，从UCSC生成TxDb**
```R
# 1. Download packages from Bioconductor
## mm10
BiocManager::install("ChIPseeker")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")
## Find interested species and txdb type here https://bioconductor.org/packages/3.11/data/annotation/
require("TxDb.Mmusculus.UCSC.mm10.knownGene")
library(clusterProfiler)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
```

* **进行注释**
```R
library(ChIPseeker)
library(org.Mm.eg.db)
library(ggplot2)
## ChIPseeker v1.28.3

#2. Load peak file
system("ls")
a4<-readPeakFile("/home/User_summits.bed")

#3. Annotate peaks
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix<-getTagMatrix(a4, windows=promoter)
anno_peak<-annotatePeak(a4,TxDb = txdb,tssRegion = c(-3000, 3000), verbose = FALSE,addFlankGeneInfo = TRUE,flankDistance = 5000, annoDb= "org.Mm.eg.db")
as.GRanges(anno_peak) %>% write.table(file = "anno_peak.tsv", quote = FALSE, sep = "\t", header = 1)

# 4. Visualization
plotAnnoBar(anno_peak)
ggsave(file = "anno_bar.pdf", height = 4)
pdf("anno_pie.pdf", height = 5,)
plotAnnoPie(anno_peak)
dev.off()


pdf("peak_heatmap.pdf", height = 10)
peakHeatmap(a4, weightCol = "V5", TxDb=txdb, upstream=3000, downstream=3000, color =rainbow(length(a4)))
dev.off()

pdf("plotAvgProf.pdf",width = 10)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")
dev.off()
pdf("plotAvgProf_conf0.95.pdf",width = 10)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            conf=0.95,resample=500,
            xlab="Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")
dev.off()
```

### 8. 差异peaks分析——DiffBind
ATAC-Seq下游分析的另一个重点是差异peaks分析。如分析不同的实验条件、多个时间节点、不同的发育时期等的差异区域。
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
