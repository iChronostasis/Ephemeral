#!/bin/bash
#SBATCH --job-name=index                # 任务名
#SBATCH -w compute07                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=index_err.log          # 日志文
### 加载软件
module load conda
module load bowtie/v2.4.4
module load TrimGalore/v0.6.7
module load Bismark/v0.23.1dev

### 循环
File1=~/bismark_result
cd $File1
mkdir 02.index
bismark_genome_preparation ~/homo_sapiens/hg38 /usr/bin/bowtie2/ --verbose $File1/02.index
