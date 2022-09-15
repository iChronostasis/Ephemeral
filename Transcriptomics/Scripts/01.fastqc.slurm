#!/bin/bash
#SBATCH --job-name=fastqc                # 任务名
#SBATCH -w compute04                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=fastqc_err.log          # 日志文
### 加载软件
module load conda

File1=/home/chengxin/2022.03.29/clean
File2=/home/chengxin/RNA-seq_slurm/List.txt
File3=/home/chengxin/2022.03.29/result

cat $File2 | while read LINE || [[ -n ${LINE} ]]
do
a=$(echo $LINE |sed 's/\n//g')
b=$(echo $a |sed 's/\r//g')
mkdir $File3/01.fastqc/$b
fastqc --t 10 -o $File3/01.fastqc/$b $File1/$b/$b\_1.fq.gz $File1/$b/$b\_2.fq.gz
done