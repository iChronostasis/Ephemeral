#!/bin/bash
#SBATCH --job-name=featurecount              # 任务名
#SBATCH -w compute04                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=subread_err.log          # 日志文
### 加载软件
module load conda

File1=/home/chengxin/2022.03.29/result/02.align

featureCounts \
	-T 16 \
	-p \
	-t exon \
	-g gene_id \
	-a /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/reads_150bp/GCF_000001635.26_GRCm38.p6_genomic.gtf \
	-o /home/chengxin/2022.03.29/result/out_feature.txt \
	$File1/*Aligned.sortedByCoord.out.bam 