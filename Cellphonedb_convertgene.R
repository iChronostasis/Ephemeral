########################   Chronostasis   ############################
######################################################################
# https://cloud.tencent.com/developer/article/1886310
#### 少量基因 ####
# 转换人鼠的基因ID
  library("biomaRt")
  human = useEnsembl(biomart="ensembl", dataset = "hsapiens_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  mouse = useEnsembl(biomart="ensembl", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
# Basic function to convert mouse to human gene names
  convertMouseGeneList <- function(x){
    genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
    return(humanx)
  }
  
  musGenes <- c("Hmmr", "Tlx3", "Cpeb4")
  musGenes <- c("Fcgr4")
  convertMouseGeneList(musGenes)

# Basic function to convert human to mouse gene names
  convertHumanGeneList <- function(x){
    genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
    humanx <- unique(genesV2[, 2])
    
    # Print the first 6 genes found to the screen
    return(humanx)
  }
  
  humGenes <- c("KLRF1", "NCAM1", "NCR1", "GNLY" ,"NKG7", "FCGR3A", "FGFBP2", "CX3CR1")
  humGenes <- c("NCAM1", "GNLY" ,"NKG7")
  humGenes <- c("FCGR3A", "KLRF1")
  
  humGenes <- c("cSF1R", "CD68","CD163", "CD14")
  humGenes <- c("CD40","CD68", "HLA-DR","CD80","CD86","FCGR1A","FCGR2A")
  humGenes <- c("CD163", "CCL18","MRC1")
  
  humGenes <- c("CD3D", "CD3E", "CD3G", "CD4", 'CD8A', 'CD8' )
  
  humGenes <- c("CD122")
  convertHumanGeneList(humGenes)
  
#### 批量+单细胞数据转换 ####
#############################
  
##### 从源头转小鼠基因到人类基因 #####  
  
# DownloadMouse reference dataset required for Cell Ranger
  #wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
# Offer downloadable file function
  library(tidyverse)
  library(DT)
  create_dt <- function(x){
    DT::datatable(x,
                  extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('csv', 'excel', 'pdf'),
                                 lengthMenu = list(c(10,25,50,-1), c(10,25,50,"All"))))
  }

# 读取features.tsv文件
  gene_list <- read.table('features.tsv', header=F, sep='\t')
  create_dt(gene_list)
# 转小鼠基因名为人类基因名
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = gene_list$V2, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  genesV2 <- genesV2[!duplicated(genesV2[,2]),] 
  create_dt(genesV2)
  
# 保存文件
  write.table(genesV2, "mouse_to_human_genes.txt", sep="\t", row.names=F, quote=F)

##### 应用转换后的小鼠到人类基因表 #####    
# 读取文件
  genesV2 <- read.table("mouse_to_human_genes.txt", sep="\t", header=T)
  
# 准备小鼠数据集
  library(Seurat)
  library(SeuratData)
# To see a manifest of all available datasets
  AvailableData()
  
# Choose small mouse dataset
  InstallData('stxKidney')
  data('stxKidney')
  stxKidney
  
## To accelerate, sample 100 cells
  seurat_object <- stxKidney
  sub_cells <- subset(seurat_object, cells = sample(Cells(seurat_object), 100))
  sub_cells@assays$Spatial[1:4,1:4]
## Extract Expression Data
  sp1 <- sub_cells
  sp1_counts <- as.matrix(sp1@assays$Spatial@data) # Notice: 这里应用的是空间数据，常规转录组数据提取用sp1@assays$RNA@data
  sp1_counts <- data.frame(gene=rownames(sp1_counts), sp1_counts, check.names = F)
  dim(sp1_counts)
  
# 转小鼠基因名为人类基因名
  sp1_counts$Gene <- genesV2[match(sp1_counts$gene, genesV2[,1]),2]
  sp1_counts <- subset(sp1_counts, Gene!='NA')
  sp1_counts <- dplyr::select(sp1_counts, Gene, everything())
  sp1_counts <- sp1_counts[, !(colnames(sp1_counts) %in% 'gene')]
  dim(sp1_counts)
  
# 保存文件
  write.table(sp1_counts, "sp1_counts_human.txt", row.names=F, sep='\t', quote=F)
  
  
  
  
  