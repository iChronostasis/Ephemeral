#!/bin/bash
#SBATCH --job-name=fastqc                # 任务名
#SBATCH -w compute06                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=fastqc_err.log          # 日志文
### 加载软件
module load conda
module load bowtie/v2.4.4
module load TrimGalore/v0.6.7
module load Bismark/v0.23.1dev

### 循环
File1=~/bismark_result	# 结果文件
File2=~/gsa/CRA002787		# 源数据
File3=~/bismark_pipeline_slurm/SRR_Acc_List.txt	# 样本列表

cat $File3 | while read LINE || [[ -n ${LINE} ]]
do
a=$(echo $LINE |sed 's/\n//g')
b=$(echo $a |sed 's/\r//g')
fastqc --t 10 -o $File1/01.fastqc $File2/$b/$b\_f1.fastq.gz $File2/$b/$b\_r2.fastq.gz
done