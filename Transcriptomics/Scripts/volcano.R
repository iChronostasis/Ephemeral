  library(ggplot2)
  library(ggpubr)
# 读取数据
  DEG_data<- as.data.frame(res)
  DEG_data$logP <- -log10(DEG_data$padj) # 对差异基因矫正后p-value进行log10()转换
  dim(DEG_data)
#将基因分为三类：not-siginficant，up，dowm
#将adj.P.value小于0.05，logFC大于2的基因设置为显著上调基因
#将adj.P.value小于0.05，logFC小于-2的基因设置为显著上调基因
  DEG_data$Group <- "not-siginficant"
  DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange >= 0)] = "up-regulated"
  DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange < 0)] = "down-regulated"
  table(DEG_data$Group)
  ggscatter(DEG_data,x = "log2FoldChange",y = "logP",
            color = "Group",
            palette = c("green","gray","red"),
            repel = T,
            ylab = "-log10(Padj)",
            size = 1) + 
    scale_y_continuous(limits = c(0,8))+
    scale_x_continuous(limits = c(-18,18))+
    geom_hline(yintercept = 1.3,linetype = "dashed")+
    geom_vline(xintercept = c(-2,2),linetype = "dashed")
### 添加label
  DEG_data$Label = ""
#### 按p值大小排序
  DEG_data <- DEG_data[order(DEG_data$padj),]
#### 高表达基因中，选择fdr最小的10个
  up.genes <- head(rownames(DEG_data)[which(DEG_data$Group == "up-regulated")],10)
#### 低表达基因中，选择fdr最小的10个
  down.genes <- head(rownames(DEG_data)[which(DEG_data$Group == "down-regulated")],10)
### merge
  deg_top10_genes<-c(as.character(up.genes),as.character(down.genes))
  DEG_data$Label[match(deg_top10_genes,rownames(DEG_data))]<-deg_top10_genes
#添加特定基因label
  ggsave("MS4A4A-vs-WT_volcano.png",vol)
  vol<-ggscatter(DEG_data,x = "log2FoldChange",y = "logP",
            color = "Group",
            palette = c("#2f5688","#BBBBBB","#CC0000"),
            label = DEG_data$Label,
            font.label = 8,
            repel = T,
            ylab = "-log10(Padj)",
            size = 1,
            title = "MS4A4A-vs-WT_volcano_plot") + 
    theme(element_line(size = 0),element_rect(size = 1.5))+ #坐标轴线条大小设置
    scale_y_continuous(limits = c(0,8))+
    scale_x_continuous(limits = c(-18,18))+
    geom_hline(yintercept = 1.3,linetype = "dashed")+
    geom_vline(xintercept = c(-2,2),linetype = "dashed")
  dev.off()


