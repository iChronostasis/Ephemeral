#!/bin/bash
#SBATCH --job-name=deduplicate               # 任务名
#SBATCH -w compute06                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=deduplicate_err.log          # 日志文
### 加载软件
module load conda
module load bowtie/v2.4.4
module load TrimGalore/v0.6.7
module load Bismark/v0.23.1dev

###循环
File1=~/bismark_result
File2=~/bismark_pipeline_slurm/SRR_Acc_List.txt
File3=~/gsa/CRA002787

cat $File2 | while read LINE || [[ -n ${LINE} ]]
do
a=$(echo $LINE |sed 's/\n//g')
b=$(echo $a |sed 's/\r//g')
deduplicate_bismark \
--bam \
--paired  $File1/04.align/$b\_f1_val_1_bismark_bt2_pe.bam \
--output_dir $File1/05.deduplicate
done

