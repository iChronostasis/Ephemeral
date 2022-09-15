### 区分结果
  selectQ<-function(data,taskname,Q_value){
    setwd("/Users/hecate/研一/RNA-seq/DESeq2/XMU")
    data1<-data[which(data[,4]<=Q_value),]
    up_data1<-data1[which(data1[,2]>=0),]
    down_data1<-data1[which(data1[,2]<0),]
    write.csv(data1,paste0(taskname,Q_value,".csv"))
    write.csv(up_data1,paste0("up",taskname,Q_value,".csv"))
    write.csv(down_data1,paste0("down",taskname,Q_value,".csv"))
  }

  tableQ<-function(data,taskname,Q_value){
    table(data[,5]<=Q_value)
    table(data[,5]<=Q_value&data[,4]>=0)
    table(data[,5]<=Q_value&data[,4]<0)
  }


### GO-BP 
  GO_BP<-function(taskname,dir_path){
    library(DOSE)
    library(org.Hs.eg.db)
    library(topGO)
    library(clusterProfiler)
    library(pathview)
    library(org.Mm.eg.db)
    gene<-read.csv(paste0(dir_path,taskname,".csv"),sep=",",header=T)[,-1]
    go <- enrichGO(gene[,1], OrgDb = org.Mm.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                   qvalueCutoff = 0.2,keyType = 'SYMBOL')
    dot<-dotplot(go,showCategory=5,title=paste0(taskname,"_GO_BP"))
    SaveDR <- c(dir_path) 
    write.csv(go,paste0(SaveDR,taskname,"_GO_BP",".csv"))
    ggsave(paste0(SaveDR,taskname,"_GO_BP",".png"),dot)
    pdf(file = paste0(SaveDR,taskname,"_GO_BP",".pdf"))
    print(dot)
    dev.off()
  }


### GSEA 
  GSEA_plot<-function(taskname,type,dir_path,reference_path,p_value){
    library(DOSE)
    library(org.Hs.eg.db)
    library(topGO)
    library(clusterProfiler)
    library(pathview)
    library(org.Mm.eg.db)
    library(annotationTools)
    library(biomaRt)
    reference<-read.gmt(reference_path)
    data<-read.csv(paste0(dir_path,taskname,".csv"),sep=",",header=T)[,-1]
    geneList = as.numeric(data[,2])#把foldchange按照从大到小提取出来
    names(geneList) <- as.character(data[,1]) #给上面提取的foldchange加上对应上ENTREZID
    
    FCgenelist <- sort(geneList,decreasing=T) #decreasing order
  #GSEA
    gsea <- GSEA(FCgenelist,
                 TERM2GENE = reference, verbose=FALSE, pvalueCutoff = p_value,pAdjustMethod = "BH") #GSEA分析
    head(gsea)
    SaveDR <- c(dir_path) 
    write.csv(gsea,paste0(SaveDR,taskname,type,".csv"))
    dot<-dotplot(gsea,title=paste0(taskname,type))
    ggsave(paste0(SaveDR,taskname,type,".png"),dot)
    pdf(file = paste0(SaveDR,taskname,type,".pdf"))
    print(dot)
    dev.off()
  }


### 画图
  volcano<-function(DEG_data,taskname,dir_path){
    library(ggplot2)
    library(ggpubr)
    # 读取数据
    DEG_data<- as.data.frame(DEG_data)
    DEG_data$logP <- -log10(DEG_data$padj) # 对差异基因矫正后p-value进行log10()转换
    dim(DEG_data)
    #将基因分为三类：not-siginficant，up，dowm
    #将adj.P.value小于0.05，logFC大于0的基因设置为显著上调基因
    #将adj.P.value小于0.05，logFC小于0的基因设置为显著上调基因
    DEG_data$Group <- "not-siginficant"
    DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange >= 0)] = "up-regulated"
    DEG_data$Group[which((DEG_data$padj < 0.05) & DEG_data$log2FoldChange < 0)] = "down-regulated"
    table(DEG_data$Group)
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
    SaveDR <- c(dir_path)
    vol<-ggscatter(DEG_data,x = "log2FoldChange",y = "logP",
                   color = "Group",
                   palette = c("#2f5688","#BBBBBB","#CC0000"),
                   label = DEG_data$Label,
                   font.label = 8,
                   repel = T,
                   ylab = "-log10(Padj)",
                   size = 1,
                   title = paste0(taskname,"_volcano_plot")) + 
      theme(element_line(size = 0),element_rect(size = 1.5))+ #坐标轴线条大小设置
      scale_y_continuous(limits = c(0,8))+
      scale_x_continuous(limits = c(-18,18))+
      geom_hline(yintercept = 1.3,linetype = "dashed")+
      geom_vline(xintercept = c(-2,2),linetype = "dashed")
    ggsave(paste0(SaveDR,taskname,"_volcano.png"),vol)
    pdf(file = paste0(SaveDR,taskname,"_volcano.pdf"))
    print(vol)
    dev.off()
  }


###venn图
  vennplot<-function(Gene1,Gene2,taskname,Q_value,dir_path){
    setwd(dir_path)
    library(ggplot2)
    library(VennDiagram)
    DEGs1<-Gene1[which(Gene1[,5]<Q_value),]
    DEGs2<-Gene2[which(Gene2[,5]<Q_value),]
    venn_list <- list(D10 = DEGs1$gene, D12 = DEGs2$gene)
    venn.diagram(venn_list, filename = paste0("all",taskname,"_venn.png"), imagetype = 'png', 
                 fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
                 col = 'black', cex = 1.5, fontfamily = 'serif', 
                 cat.cex = 1.5, cat.fontfamily = 'serif',main=paste0("all_",taskname,"_venn"))
    Gene1<-D10
    Gene2<-D12
    Q_value = 0.05
    dir_path = "/Users/hecate/研一/Term 2/2022.03.31/venn"
    up1<-Gene1[which(Gene1[,5]<Q_value&Gene1[,3]>=0),]
    up2<-Gene2[which(Gene2[,5]<Q_value&Gene2[,3]>=0),]
    #write.csv(up1,"upDEGsD10.csv")
    #write.csv(up2,"upDEGsD12.csv")
    venn_list <- list(D10 = up1$gene, D12 = up2$gene)
    venn.diagram(venn_list, filename = paste0("up",taskname,"_venn.png"), imagetype = 'png', 
                 fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
                 col = 'black', cex = 1.5, fontfamily = 'serif', 
                 cat.cex = 1.5, cat.fontfamily = 'serif',main=paste0("up_",taskname,"_venn"))
    
    down1<-Gene1[which(Gene1[,5]<Q_value&Gene1[,3]<0),]
    down2<-Gene2[which(Gene2[,5]<Q_value&Gene2[,3]<0),]
    #write.csv(down1,"downDEGsD10.csv")
    #write.csv(down2,"downDEGsD12.csv")
    venn_list <- list(D10 = down1$gene, D12 = down2$gene)
    venn.diagram(venn_list, filename = paste0("down",taskname,"_venn.png"), imagetype = 'png', 
                 fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
                 col = 'black', cex = 1.5, fontfamily = 'serif', 
                 cat.cex = 1.5, cat.fontfamily = 'serif',main=paste0("down_",taskname,"_venn"))
  }


### intersect
  overlap<-function(Gene1,Gene2,taskname,Q_value,dir_path){
    Gene2<-as.data.frame(Gene2)
    DEGs1<-Gene1[which(Gene1[,5]<Q_value),]
    DEGs2<-Gene2[which(Gene2[,6]<Q_value),]
    num1<-length(intersect(rownames(DEGs2),DEGs1$gene_symbol))
    print(num1/nrow(DEGs1))
    up1<-Gene1[which(Gene1[,5]<Q_value&Gene1[,4]>=0),]
    up2<-Gene2[which(Gene2[,6]<Q_value&Gene2[,2]>=0),]
    num2<-length(intersect(rownames(up2),up1$gene_symbol))
    print(num2/nrow(up1))
    down1<-Gene1[which(Gene1[,5]<Q_value&Gene1[,4]<0),]
    down2<-Gene2[which(Gene2[,6]<Q_value&Gene2[,2]<0),]
    num3<-length(intersect(rownames(down2),down1$gene_symbol))
    print(num3/nrow(down1))
  }









