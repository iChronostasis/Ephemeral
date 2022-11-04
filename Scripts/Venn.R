########################   Chronostasis   ##########################################
####################################################################################
# About Venn and intersection
#--------------------------------------------------------------
## 1. The Simplest version
  library(ggplot2)
  library(VennDiagram)
### get the list of venn
  venn_list <- list(D10 = DEGs1$gene, D12 = DEGs2$gene)
  venn.diagram(venn_list, filename = paste0("all_venn.png"), imagetype = 'png', 
              fill = c('red', 'blue'), alpha = 0.50, cat.col = rep('black', 2), 
              col = 'black', cex = 1.5, fontfamily = 'serif', 
              cat.cex = 1.5, cat.fontfamily = 'serif',main=paste0("all_",taskname,"_venn"))
#--------------------------------------------------------------
## 2. The complicated version[Like a flower]
#从CRAN安装ggVennDiagram包；
  install.packages("ggVennDiagram")
  install.packages("ggsci")

#载入所需的R包；
  library(ggplot2)
  library(ggsci)
  library(sf)
  library(ggVennDiagram)

#自定义颜色；
  color1 <- alpha("#f8766d",0.9)
  color2 <- alpha("#FF99CC",0.7)
  color3 <- alpha("#c77cff",0.5)
  color4 <- alpha("#99CC00",0.5)

#绘制常见的4组数据venn图；
  ggVennDiagram(x[1:4], label_alpha=0) +
    scale_fill_gradient(low="white",high =color4 ,guide="none")
  p <- ggVennDiagram(venn_list,size=1,,lty="longdash") +
    scale_fill_gradient(low="white",high = color4 ,guide="none") +
    labs(title = "Overlap of highly variable features",
         subtitle = "After the standard pre-processing workflow") +
    theme(plot.title = element_text(hjust = 0.5,vjust = .5,color = "black",face = 'bold',
                                    size = 20, margin = margin(t = 1, b = 12)),
          plot.subtitle = element_text(hjust = 0,vjust = .5,size=15),
          plot.caption = element_text(face = 'bold',size = 12),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  ggsave("venn.pdf", plot = p, width = 15, height = 10)  

#--------------------------------------------------------------
### Visualize the venn in another type[4 dataset]
  venn_list <- readRDS(paste0(out_data_dir,"venn_list_QC.rds"))
  seurat_list <- readRDS(paste0(out_data_dir,"seurat_list_QC.rds"))
  install.packages("UpSetR")
  library(UpSetR)
  pdf(paste0(fig_dir,"Visualization/UpSet_QC_4.pdf"), width=12, height=10)
  upset(fromList(venn_list),text.scale = 2,
        point.size = 3.5, 
        line.size = 2,  order.by = "freq")
  # upset(fromExpression(expressionInput), order.by = "freq")
  dev.off()

