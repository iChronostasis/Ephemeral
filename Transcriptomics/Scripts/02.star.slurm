#!/bin/bash
#SBATCH --job-name=align                # 任务名
#SBATCH -w compute04                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=align_err.log          # 日志文
### 加载软件
module load conda
File1=/home/chengxin/2022.03.29/clean
File2=/home/chengxin/RNA-seq_slurm/List1.txt
mkdir /home/chengxin/2022.03.29/result/02.align
File3=/home/chengxin/2022.03.29/result/02.align
cd $File3
cat $File2 | while read LINE || [[ -n ${LINE} ]]
do
a=$(echo $LINE |sed 's/\n//g')
b=$(echo $a |sed 's/\r//g')
mkdir $b
STAR \
--genomeDir /share/data0/reference/STAR_genome_index/Mus_musculus/GRCm38.p6/reads_150bp \
--readFilesIn $File1/$b/$b\_1.fq.gz $File1/$b/$b\_2.fq.gz \
--readFilesCommand zcat \
--runThreadN 8 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--outFileNamePrefix $File3/$b

done