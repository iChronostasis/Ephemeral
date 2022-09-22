## ♪( ´▽｀)
# About Ephemeral
<mark>**仅供本实验室**</mark><br>
2022.05.07 暂时地存储一些表观遗传学的RRBS技术和处理方法bismark流程的教程及代码。<br>
2022.06.12-13 更新了Chip-seq的流程内容,Transcriptomics的内容和GSEA的内容。<br>
2022.07.15-21 更新了CUT_Tag_pipeline，ATAC-seq和ChIP-seq pipeline的相关内容，Software_Parameters.md为了未来封装流程整理的各种软件参数。ATAC-seq和Chip-seq分析流程基本一致，具体个性化分析可以参考CUT_Tag_pipeline的部分内容。<br>
2022.09.22 更新了学习宏基因组相关内容。<br>

### To be continued...
# User Guide
1. **Epigenomics** 关于表观遗传学
    * **RRBS.md** 关于RRBS的技术和概念总结
    * **bismark.md** bismark的运行教程
    * **Chip-seq_pipeline.md** 关于Chip-seq的运行教程及代码
    * **ATAC-seq_pipeline.md** 关于ATAC-seq的运行教程及代码
    * **CUT_TAG_pipeline.md** 关于CUT&Tag_data的运行教程及代码
    * **Software_Parameters.md** 为了未来封装流程整理的各种软件参数
    * **Scripts**
        * **Bismark_Pipeline** 存储bismark的运行脚本
2. **Transcriptomics** 关于转录组学
    * **RNA-seq_pipeline.md** RNA-seq的运行教程及代码
    * **转录组分析报告.pdf** 分析结果报告模板
    * **Scripts**
        * **01.fastqc.slurm，02.star.slurm，03.feature_count.slurm** 服务器中使用的脚本
        * **RNA-seq_FPKM_visualization.R** R脚本
        * **test_function.R** 批量处理时：选择阈值，GO-BP，GSEA，volcano火山图，Venn图，获得差异分析后的两组差异基因的交集的部分函数代码合集
        * **venn.R** 使用**test_function.R**的批量处理的例子
        * **volcano.R** 使用**test_function.R**的批量处理的例子
        * **GO_KEGG.R** 使用**test_function.R**的批量处理的例子
3. **GSEA** 基因功能富集(均下载于GSEA官方网站)
    * **h.all.v7.4.symbols.gmt** hallmark数据库symbol的基因注释reference
    * **c2.cp.reactome.v7.4.symbols.gmt** reactome数据库symbol的基因注释reference
    * **MS4A4A-医学统计支持办公室.pptx** GSEA网站使用方法，同时也是医学统计支持办公室结果汇总的模板
4. **Metagenomnic** 关于宏基因组分析的workflow的总结
    * **Metagenomic/Metagenomic_StudyAccount.md** 关于宏基因组的背景知识介绍
    * **Metagenomic/Metagenomic_shotgun-seq_pipeline.md** shotgun-seq 分析流程【未上传，汇报完再说】
5. **Cellphonedb_convertgene.R** 小鼠和人类的Gene Name转换
# w(ﾟДﾟ)w 要好好学习啊！
