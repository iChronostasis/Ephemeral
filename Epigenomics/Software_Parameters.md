# Software[Epigenomics]
<!-- TOC -->

- [Software[Epigenomics]](#softwareepigenomics)
    - [1.  Quality Control](#1--quality-control)
        - [1.1 [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)](#11-fastqchttpsdnacoremissouriedupdffastqc_manualpdf)
        - [1.2 [Multiqc](https://multiqc.info/docs/)](#12-multiqchttpsmultiqcinfodocs)
    - [2. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)](#2-trimmomatichttpwwwusadellaborgcmspagetrimmomatic)
    - [3.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-I)](#3--bowtie2httpbowtie-biosourceforgenetbowtie2manualshtmlbowtie2-options-i)
    - [4.  [Samtools](http://www.htslib.org/)](#4--samtoolshttpwwwhtsliborg)
    - [5.  [bedtools](https://bedtools.readthedocs.io/en/latest/)](#5--bedtoolshttpsbedtoolsreadthedocsioenlatest)
    - [6. [GATK](https://gatk.broadinstitute.org/hc/en-us)](#6-gatkhttpsgatkbroadinstituteorghcen-us)
    - [7.  [deeptools](https://deeptools.readthedocs.io/en/develop/)](#7--deeptoolshttpsdeeptoolsreadthedocsioendevelop)
    - [8.  [MACS2 callpeak](https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md)](#8--macs2-callpeakhttpsgithubcommacs3-projectmacsblobmasterdocscallpeakmd)

<!-- /TOC -->
## 1.  Quality Control
### 1.1 [FastQC](https://dnacore.missouri.edu/PDF/FastQC_Manual.pdf)
```bash
Usage:
fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam]
           [-c contaminant file] seqfile1 .. seqfileN
Options:
  -t 线程数
  -o 输出目录
  # 其余参数选项为default 
```
### 1.2 [Multiqc](https://multiqc.info/docs/)
```bash
#MultiQC aggregates results from bioinformatics analyses across many samples into a single report.
Usage: multiqc [OPTIONS] [ANALYSIS DIRECTORY]
Options:
  -o 输出目录
```
## 2. [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
```bash
Usage:
       PE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] [-validatePairs] [-basein <inputBase> | <inputFile1> <inputFile2>] [-baseout <outputBase> | <outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] <trimmer1>...
   or:
       SE [-version] [-threads <threads>] [-phred33|-phred64] [-trimlog <trimLogFile>] [-summary <statsSummaryFile>] [-quiet] <inputFile> <outputFile> <trimmer1>...
   or:
       -version
Options:
  PE 双端测序
  SE 单端测序
  -baseout 输出目录
  -threads 线程数
  ILUMINACLIP 找到并去除Illumina adapters（find and remove Illumina adapters）
  Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
  Remove leading low quality or N bases (below quality 3) (LEADING:3)
  Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
  Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
  Drop reads below the 36 bases long (MINLEN:36)
```

## 3.  [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#bowtie2-options-I)
```bash
1. bowtie2 
Usage:
  bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | -b <bam>} [-S <sam>]
Options:
  -end-to-end [Defalut] That is, it searches for alignments involving all of the read characters. This is also called an "untrimmed" or "unclipped" alignment.
  --very-sensitive 通常以较慢的方式结束但更灵敏、更准确
  --no-mixed 默认情况下，当`bowtie2`无法为一对找到一致或不一致的对齐时，它会尝试为各个配合找到对齐。此选项禁用该行为。
  --no-discordant 默认情况下，`bowtie2`如果找不到任何一致的对齐方式，则查找不一致的对齐方式。不一致对齐是两个配对唯一对齐但不满足配对末端约束 ( [`--fr`/ `--rf`/`--ff`]的对齐。此选项禁用该行为。
  --phred33 输入质量是等于[Phred 质量]加 33 的 ASCII 字符。
  -I 有效配对末端比对的最小片段长度。
  -X 有效配对末端比对的最大片段长度。

2. bowtie2-build
Usage: bowtie2-build [options]* <reference_in> <bt2_index_base>
    reference_in            comma-separated list of files with ref sequences
    bt2_index_base          write bt2 data to files with this dir/basename

  -f reference files are Fasta (default)
  --threads 线程数
```

## 4.  [Samtools](http://www.htslib.org/)
```bash
1. Samtools
Usage:   samtools <command> [options]
Commands:
  -- Indexing
     dict           create a sequence dictionary file
     faidx          index/extract FASTA
     index          index alignment

  -- Editing
     calmd          recalculate MD/NM tags and '=' bases
     fixmate        fix mate information
     reheader       replace BAM header
     rmdup          remove PCR duplicates
     targetcut      cut fosmid regions (for fosmid pool only)
     addreplacerg   adds or replaces RG tags

  -- File operations
     collate        shuffle and group alignments by name
     cat            concatenate BAMs
     merge          merge sorted alignments
     mpileup        multi-way pileup
     sort           sort alignment file
     split          splits a file by read group
     quickcheck     quickly check if SAM/BAM/CRAM file appears intact
     fastq          converts a BAM to a FASTQ
     fasta          converts a BAM to a FASTA

  -- Statistics
     bedcov         read depth per BED region
     depth          compute the depth
     flagstat       simple stats
     idxstats       BAM index stats
     phase          phase heterozygotes
     stats          generate stats (former bamcheck)

  -- Viewing
     flags          explain BAM flags
     tview          text alignment viewer
     view           SAM<->BAM<->CRAM conversion
     depad          convert padded BAM to unpadded BAM


2. samtools view
# 在未指定选项或区域的情况下，将指定输入对齐文件（SAM、BAM 或 CRAM 格式）中的所有对齐打印到 SAM 格式（无标题）的标准输出。还可以将SAM格式文件放入管道转变为BAM文件格式。
Usage: samtools view [options] <in.bam>|<in.sam>|<in.cram> [region ...]
Options:
  -b output BAM
  -h include header in SAM output
  -o FILE  output file name [stdout]
  -S ignored (input format is auto-detected)
  -h, --with-header Include header in SAM output
  -f, --require-flags FLAG have all of the FLAGs present
  -F, --exclude-flags FLAG ...have none of the FLAGs present

3. samtools index
# 索引坐标排序的 BGZIP 压缩 SAM、BAM 或 CRAM 文件，以实现快速随机访问。
Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
Options:
  -b       Generate BAI-format index for BAM files [default]
  -c       Generate CSI-format index for BAM files
  -m INT   Set minimum interval size for CSI indices to 2^INT [14]
  -@ INT   Sets the number of threads [none]
```
## 5.  [bedtools](https://bedtools.readthedocs.io/en/latest/)
```bash
1. [ Genome arithmetic ]
# intersect Find overlapping intervals in various ways.
Usage: bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>
Options:
  -a 参数后加重复样本1（A）
  -b 参数后加重复样本2（B），也可以加多个样本
2. [ Format conversion ]
# bedtobam Convert intervals to BAM records.
Usage:   bedtools bedtobam [OPTIONS] -i <bed/gff/vcf> -g <genome>
Options:
  -i input
  -bedpe Set the score field based on BAM tags
  ## 根据BAM标签设置score字段将BAM比对转换为BEDPE格式，从而允许在单个文本行上报告双端比对的两端
```
## 6. [GATK](https://gatk.broadinstitute.org/hc/en-us)
```bash
1. gatk SortSam
USAGE: SortSam [arguments]
Required Arguments:
  --INPUT,-I <File> Input BAM or SAM file to sort.
  --OUTPUT,-O <File> Sorted BAM or SAM output file.
  --SORT_ORDER Sort order of output file.
  ## coordinate (Sorts primarily according to the SEQ and POS fields of the record.)
Optional Arguments:
  --VALIDATION_STRINGENCY Validation stringency for all SAM files read by this program.
  ## Default value: STRICT.Possible values: {STRICT, LENIENT, SILENT}

2. gatk MarkDuplicates
USAGE: MarkDuplicates [arguments]
Required Arguments:
  --INPUT,-I <File> Input BAM or SAM file to sort.
  --OUTPUT,-O <File> Sorted BAM or SAM output file.
  --METRICS_FILE,-M <File> File to write duplication metrics to
Optional Arguments:
  --REMOVE_DUPLICATES <Boolean> If true do not write duplicates to the output file instead of writing them with appropriate flags set.
  ## REMOVE_DUPLICATES=true
```
## 7.  [deeptools](https://deeptools.readthedocs.io/en/develop/)
```bash
1. bamCoverage Convert into bigwig file format
usage: An example usage is:$ bamCoverage -b reads.bam -o coverage.bw
Required arguments:
  --bam BAM file, -b BAM file BAM file to process (default: None)
Output:
  --outFileName FILENAME, -o FILENAME
                        Output file name. (default: None)
  --outFileFormat {bigwig,bedgraph}, -of {bigwig,bedgraph}
                        Output file type. Either "bigwig" or "bedgraph". (default:
                        bigwig)
2. multiBamSummary bins
usage: multiBamSummary bins --bamfiles file1.bam file2.bam -o results.npz
Required arguments:
  --bamfiles FILE1 FILE2 [FILE1 FILE2 ...], -b FILE1 FILE2 [FILE1 FILE2 ...]
                        List of indexed bam files separated by spaces. (default: None)
  --outFileName OUTFILENAME, -out OUTFILENAME, -o OUTFILENAME
                        File name to save the coverage matrix. This matrix can be
                        subsequently plotted using plotCorrelation or or plotPCA.
                        (default: None)

Optional arguments:
  --labels sample1 sample2 [sample1 sample2 ...], -l sample1 sample2 [sample1 sample2 ...]User defined labels instead of default labels from file names.Multiple labels have to be separated by a space, e.g. --labels sample1 sample2 sample3 (default: None)
  --binSize INT, -bs INT Length in bases of the window used to sample the genome.(Default: 10000)
  --numberOfProcessors INT, -p INT Number of processors to use. Type "max/2" to use half the maximum number of processors or "max" to use all available processors.(Default: 1)

Output optional options:
  --outRawCounts FILE   Save the counts per region to a tab-delimited file. (default:None)
  --scalingFactors FILE Compute scaling factors (in the DESeq2 manner) compatible for use with bamCoverage and write them to a file. The file has tab-separated columns "sample" and "scalingFactor". (default: None)
4. plotCorrelation
usage: plotCorrelation [-h] --corData FILE --corMethod {spearman,pearson} --whatToPlot
                       {heatmap,scatterplot} [--plotFile FILE] [--skipZeros]
                       [--labels sample1 sample2 [sample1 sample2 ...]]
                       [--plotTitle PLOTTITLE] [--plotFileFormat FILETYPE]
                       [--removeOutliers] [--version] [--outFileCorMatrix FILE]
                       [--plotHeight PLOTHEIGHT] [--plotWidth PLOTWIDTH] [--zMin ZMIN]
                       [--zMax ZMAX] [--colorMap] [--plotNumbers]
                       [--xRange XRANGE XRANGE] [--yRange YRANGE YRANGE] [--log1p]
Required arguments:
  --corData FILE, -in FILE
                        Compressed matrix of values generated by multiBigwigSummary or multiBamSummary
  --corMethod {spearman,pearson}, -c {spearman,pearson} Correlation method.
  --whatToPlot {heatmap,scatterplot}, -p {heatmap,scatterplot} Choose between a heatmap or pairwise scatter plots

Output optional options:
  --outFileCorMatrix FILE Save matrix with pairwise correlation values to a tab-separated file.

Heatmap options:
  --plotHeight PLOTHEIGHT Plot height in cm. (Default: 9.5)
  --plotWidth PLOTWIDTH Plot width in cm. The minimum value is 1 cm. (Default: 11)
  --zMin ZMIN, -min ZMIN Minimum value for the heatmap intensities. If not specified, the value is set automatically
  --zMax ZMAX, -max ZMAX Maximum value for the heatmap intensities.If not specified, the value is set automatically
  --colorMap Color map to use for the heatmap. Available values can be seen here:
http://matplotlib.org/examples/color/colormaps_reference.html
  --plotNumbers If set, then the correlation number is plotted on top of the heatmap. This option is only valid when plotting a heatmap.

6. computeMatrix reference-point
usage: An example usage is:
  computeMatrix reference-point -S <bigwig file(s)> -R <bed file(s)> -b 1000
  
Required arguments:
  --regionsFileName File [File ...], -R File [File ...] File name or names, in BED or GTF format, containing the regions to plot.
  --scoreFileName File [File ...], -S File [File ...] bigWig file(s) containing the scores to be plotted.

Output options:
  --outFileName OUTFILENAME, -out OUTFILENAME, -o OUTFILENAME
                        File name to save the gzipped matrix file needed by the
                        "plotHeatmap" and "plotProfile" tools. (default: None)
Optional arguments:
  --skipZeros Whether regions with only scores of zero should be included or not. Default is to include them. (default: False)
  --referencePoint Possible choices: TSS, TES, center
#The reference point for the plotting could be either the region start (TSS), the region end (TES) or the center of the region. Note that regardless of what you specify, plotHeatmap/plotProfile will default to using “TSS” as the label. (Default: “TSS”)
  --beforeRegionStartLength, -b, --upstream Distance upstream of the reference-point selected. (Default: 500)
  --afterRegionStartLength, -a, --downstream Distance downstream of the reference-point selected. (Default: 1500)
  --numberOfProcessors, -p Number of processors to use.
```

## 8.  [MACS2 callpeak](https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md)
```bash
usage: macs2 callpeak [-h] -t TFILE [TFILE ...] [-c [CFILE [CFILE ...]]]
                      [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}]
                      [-g GSIZE] [-s TSIZE] [--keep-dup KEEPDUPLICATES]
                      [--outdir OUTDIR] [-n NAME] [-B] [--verbose VERBOSE]
                      [--trackline] [--SPMR] [--nomodel] [--shift SHIFT]
                      [--extsize EXTSIZE] [--bw BW] [--d-min D_MIN]
                      [-m MFOLD MFOLD] [--fix-bimodal] [-q QVALUE | -p PVALUE]
                      [--scale-to {large,small}] [--ratio RATIO]
                      [--down-sample] [--seed SEED] [--tempdir TEMPDIR]
                      [--nolambda] [--slocal SMALLLOCAL] [--llocal LARGELOCAL]
                      [--max-gap MAXGAP] [--min-length MINLEN] [--broad]
                      [--broad-cutoff BROADCUTOFF] [--cutoff-analysis]
                      [--call-summits] [--fe-cutoff FECUTOFF]
                      [--buffer-size BUFFER_SIZE] [--to-large]
Required arguments:
  -t IP 数据文件（这是 MACS 的唯一必需参数）
  -c 控制或模拟数据文件
  -f 输入文件的格式；默认为“AUTO”，这将允许 MACS 自动决定格式。
  -g 可映射的基因组大小，定义为可以测序的基因组大小；提供了一些预编译的值。

Output options:
  --outdir  MACS2 会将所有输出文件保存到此选项的指定文件夹中
  -n 输出文件的前缀字符串
  -B/--bdg 将片段堆积、控制 lambda、-log10pvalue 和 -log10qvalue 分数存储在 bedGraph 文件中

Shifting model arguments 转移模型参数:
  -s 测序标签的大小。默认情况下，MACS 将使用输入处理文件中的前 10 个序列来确定它
  --bw 仅用于模型构建的基因组扫描带宽。可以设置为预期的超声处理片段大小。
  --mfold 模型构建的上限和下限

Peak calling arguments 峰值调用参数:
  -q q 值（最小 FDR）截止
  -p p 值截止（而不是 q 值截止）
  --nolambda 不考虑峰值候选区域的局部偏差/lambda
  --broad  广泛的峰值呼叫
```
