#!/bin/bash
#SBATCH --job-name=report              # 任务名
#SBATCH -w compute06                  # 申请节点名
#SBATCH --nodes 1                     # 申请使用一个节点
#SBATCH --ntasks 20                    # 申请使用一个CPU
#SBATCH --time=3-12:00:00               # 申请使用时间 时:分:秒
#SBATCH --output=report_err.log          # 日志文
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
bismark2report \
--dir $File1/07.report \
--alignment_report $File1/04.align/$b\_f1_val_1_bismark_bt2_PE_report.txt \
--dedup_report $File1/05.deduplicate/$b\_f1_val_1_bismark_bt2_pe.deduplication_report.txt \
--splitting_report $File1/06.extractor/$b\_f1_val_1_bismark_bt2_pe.deduplicated_splitting_report.txt \
--mbias_report $File1/06.extractor/$b\_f1_val_1_bismark_bt2_pe.deduplicated.M-bias.txt 
done

