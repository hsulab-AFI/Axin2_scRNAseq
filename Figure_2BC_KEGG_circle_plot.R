
library(Seurat)
library(tidyverse)
library(rJava)
devtools::install_github('colearendt/xlsx')
library(xlsx)

testseu=readRDS("Axin2+_scRNA-seq_combined.rds")
# Idents(testseu) <- "anno_new"

### Find DEGs #########################################################################
marker_celltype <- FindAllMarkers(testseu,
               only.pos = TRUE,
               min.pct =  0.25, 
               min.diff.pct = 0.25)

# filter
marker_celltype=marker_celltype%>%filter(p_val_adj < 0.05)
marker_celltype$d=marker_celltype$pct.1-marker_celltype$pct.2
marker_celltype=marker_celltype%>%filter(d > 0.2)
marker_celltype=marker_celltype%>%arrange(cluster,desc(avg_log2FC))
marker_celltype=as.data.frame(marker_celltype)
write.xlsx(marker_celltype,file = "markers_log2fc0.25_padj0.05_pctd0.2.xlsx",row.names = F,col.names = T)

### enrichment analysis ###########################################################################
library(clusterProfiler)
library(org.Mm.eg.db)
R.utils::setOption("clusterProfiler.download.method","auto") #https://github.com/YuLab-SMU/clusterProfiler/issues/256

source("syEnrich.R")
syEnrich(marker_celltype,outpath = "markers_log2fc0.25_padj0.05_pctd0.2")

### find one cell type #######################################################
kegg.res=read.xlsx("markers_log2fc0.25_padj0.05_pctd0.2.KEGG.xls",sheetIndex = 1,as.data.frame = T,header = T)
kegg.res=kegg.res%>%filter(p.adjust < 0.05)
kegg.res=kegg.res%>%filter(cluster == "SC")

# import kegg_info
kegg_info=read.xlsx("kegg_info.xlsx",sheetIndex = 1,startRow = 3)
kegg_info=kegg_info[,c("ID","Pathway","big.annotion")]

# merge
kegg.res$ID=str_replace(kegg.res$ID,"mmu","")
kegg.res=kegg.res%>%inner_join(kegg_info,by = "ID")

# top20 term
kegg.res=kegg.res%>%arrange(p.adjust)
kegg.res=head(kegg.res,20)

saveRDS(kegg.res, file= "kegg.res.rds")

write.table(kegg.res,file = "kegg.res.txt",quote = F,sep = "\t",row.names = F,col.names = T)
write.xlsx(kegg.res,file = "kegg.res.xlsx",col.names = T,row.names = F)

### use kegg_loop as function ###############################################################
source("kegg_loop.R")
kegg_loop(enrich.res = kegg.res,filename = "SC_kegg_id")
