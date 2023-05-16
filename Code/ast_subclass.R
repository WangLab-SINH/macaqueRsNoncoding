.libPaths("D:/Program Files/R/R-4.2.0/library")
#install.packages("qs")
#BiocManager::install("SingleR")

#BiocManager::install("scran")


library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(SingleR)
library(scran)
library(reshape2)

setwd("F:/Lab_info/wanglab/ww_cellreports_comments")

v1_sc_mfas6=qread("./singleCell_jjm_mfas6/singleCell_jjm_mfas6.qs")
DimPlot(v1_sc_mfas6,reduction  = "umap", label = T)+NoLegend()
ggsave(filename = "./figures/v1_sc_mfas6_umap.pdf",width=4.98,height=4.68)
FeaturePlot(v1_sc_mfas6,features = "SLC1A2",cols = c("lightgrey", "brown"))
VlnPlot(v1_sc_mfas6,features = "NGN2")

# 
macaque_brain_sc=readRDS("snRNA-seq_macaque_brain_Seurat.RDS")
meta=macaque_brain_sc@meta.data 
fig1b_df=read.csv("fig1b_rawdata.txt",sep="\t",row.names = 1, header=T,check.names = F)
meta$celltype_6=fig1b_df$celltype_5
meta$celltype_6[meta$celltype_6=="PERI"]="ENDO"
idx=which(meta$celltype_6=="AST")
meta$celltype_6[idx]=meta$celltype_5[idx]
idx=which(meta$celltype_6=="LAMP5")
meta$celltype_6[idx]=meta$celltype_5[idx]
idx=which(meta$celltype_6=="PVALB")
meta$celltype_6[idx]=meta$celltype_5[idx]

macaque_brain_sc@meta.data=meta

table(meta$celltype_6)

DefaultAssay(macaque_brain_sc)="SCT"
macaque_brain_sc@active.assay

macaque_brain_sc_for_SingleR <- GetAssayData(macaque_brain_sc, 
                                             slot="data") 
test <- GetAssayData(v1_sc_mfas6, slot="data")
v1_sc_mfas6.hesc <- SingleR(test=test, 
                            ref=macaque_brain_sc_for_SingleR, 
                            labels=meta$celltype_6, 
                            de.method="wilcox")

table(v1_sc_mfas6.hesc$labels,v1_sc_mfas6@meta.data$cellType) 

aa=table(v1_sc_mfas6.hesc$labels,v1_sc_mfas6@meta.data$cellType)
aa= apply(aa,2,function(x) 100*x/sum(x))
df=as.data.frame(melt(aa))
df$Var2=as.factor(df$Var2)
a <- ggplot(df, aes(Var1,Var2)) + 
  geom_point(aes(size = value), colour = "brown") + 
  scale_size_continuous(range = c(0,6)) +
  theme_bw()+
  labs(x="Cell type (Lei et.al)",y="Cell type") + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
a


macaque_brain_select_genes=c("SLC17A7","GAD1","ADARB2",
                             "LHX6","SLC1A2","PDGFRA","MOG",
                             "APBB1IP","FLT1","LAMP5","RELN",
                             "PVALB","UNC5B","RORA","PDGFRB",
                             "SST","VIP")

b=DotPlot(v1_sc_mfas6, features=c("RORB","HPCAL1",macaque_brain_select_genes))+
  scale_colour_gradient2(low="#2166AC",mid="#F1F4F6",high="#B2182B") +
  scale_size_continuous(range = c(0,6)) +
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1)) +
  labs(x="genes",y="cell type")

b+a
ggsave(filename = "./figures/Celltype_assignment_to_lei et al.pdf",
       width=13.89,height=6.89)

#VlnPlot(macaque_brain_sc,features = macaque_brain_select_genes)
DotPlot(macaque_brain_sc, features =macaque_brain_select_genes,
        group.by = "celltype_6",
        cols = c("lightgrey", "brown"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))

DotPlot(v1_sc_mfas6, features=c("RORB","HPCAL1",macaque_brain_select_genes),
        cols = c("lightgrey", "brown"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))


#save.image("v1_snRNAseq_mfas6_and_macaque_brain_snRNAseq_analysis.RData")



#AST
# AST subtype
ast.sub <- subset(v1_sc_mfas6@meta.data, cellType=="AST")
ast_sc <- subset(v1_sc_mfas6, cells=row.names(ast.sub))
ast_sc@active.assay

#DefaultAssay(ast_sc)="RNA"

##PCA
ast_sc <- FindVariableFeatures(ast_sc,selection.method = "vst",nfeatures = 3000)
scale.genes <-  rownames(ast_sc)
ast_sc <- ScaleData(ast_sc, features = scale.genes)
# runPCA
ast_sc <- RunPCA(ast_sc, features = VariableFeatures(ast_sc))
ElbowPlot(ast_sc, ndims=50, reduction="pca")
ggsave("./figures/ast_subclass_ElbowPlot.pdf",width=5.89,height=3.69 )
pc.num=1:20
ast_sc <- FindNeighbors(ast_sc, dims = pc.num)
ast_sc <- FindClusters(ast_sc, resolution = 0.5)
table(ast_sc@meta.data$seurat_clusters)
# tSNE
ast_sc = RunTSNE(ast_sc, dims = pc.num)
DimPlot(ast_sc, reduction = "tsne",label = T)
#umap
ast_sc <- RunUMAP(ast_sc, dims = pc.num)
a=DimPlot(ast_sc, reduction = "umap",label=T)
b=DimPlot(ast_sc, reduction = "umap",group.by = "samp")
a+b
ggsave(filename = "./figures/ast subclass_umap.pdf",width=9.98,height=3.98)

FeaturePlot(ast_sc,features = c("SLC1A2"))
ggsave(filename = "./figures/ast subclass_umap_SLC1A2.pdf",width=4.98,height=3.98)


ast_sub_markers=FindAllMarkers(ast_sc)
View(ast_sub_markers)
write.table(ast_sub_markers,file = "AST subclass all markers.txt", sep="\t", quote=F)

ast_sub_markers %>% filter(p_val_adj<0.05 & avg_log2FC>0) %>% group_by(cluster) %>% summarise(  n = n() )

top10 <- ast_sub_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10=top10$gene
DoHeatmap(ast_sc, features = top10) #+ scale_fill_gradientn(colors = c("blue","white","red"))
ggsave(filename = "./figures/ast 3 subclass top 10 markers_heatmap.pdf",height=9.89,width=9.89)
VlnPlot(ast_sc, features=top10,ncol = 10)
DotPlot(ast_sc, features = top10,cols = c("lightgrey", "brown"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))

# classic astrocyte markers
VlnPlot(ast_sc, features=c("GFAP","S100B","NDRG2","SLC1A3","ACSS2"))
DoHeatmap(ast_sc, features=c("GFAP","S100B","NDRG2","SLC1A3","ACSS2")) 
DotPlot(ast_sc, features = c("GFAP","S100B","NDRG2","SLC1A3","ACSS2"),cols = c("lightgrey", "brown"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))

# VlnPlot(ast_sc, features=c("KCNIP4","SNCA","ANO3","ASIC2","RBFOX1",
#                            "ENSMFAG00000052452","ATP13A4","PRKCQ",
#                            "MTSS2","FGFR2"),ncol = 5)

FeaturePlot(ast_sc,features=c("GFAP") )
VlnPlot(ast_sc, features=c("GFAP","INTU","LAMA1","MBP","PLP1","LAMA3","LAMA4",
                           "MOG","APBB1IP",
                           "PCDH15","VCAN","SPOCK3"))
DotPlot(ast_sc, features = c("GFAP","INTU","LAMA1","LAMA3","LAMA4",
                             "MOG","APBB1IP",
                             "MBP","PLP1","PCDH15","VCAN","SPOCK3"),
        cols = c("lightgrey", "brown"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))


VlnPlot(ast_sc, features=c("DMRTA2","CHRDL1","MASP1","CRYM","MYBPC1","WIF1"))

# interlaminar astrocytes
# S100B 
VlnPlot(ast_sc, features=c("S100B","CRYAB","HOPX","AQP4","SLC1A3",
                           "SOX2","PAX6","NES"))



#NC 
nc_ast_subclass_markers=c("LUZP2","COL5A3","GPC5","HSD17B5","CACNB2","GFAP","DPP10","L3MBTL4","TXNIP",
                          "SPARC","ID3","SLC38A1","RERG","TNC","LAMA1","MBP","DOCK10","PLP1","ENPP2",
                          "ST18","PLCL1","SPOCK3","ELMO1","SLC24A2","CTNNA3","PCDH15","VCAN","TNR",
                          "NXPH1","LHFPL3","COL9A1","MMP16","DSCAM","CSMD3") 

DoHeatmap(ast_sc, features=intersect(nc_ast_subclass_markers,ast_sub_markers$gene)) 
VlnPlot(ast_sc, features=nc_ast_subclass_markers,ncol = 10)
DotPlot(ast_sc, features = nc_ast_subclass_markers,
        cols = c("lightgrey", "brown"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))


ast_selected=c("CRYAB","HOPX","AQP4","SLC1A3","SOX2","PAX6",
               "INTU","LAMA1","MBP","PLP1","LAMA3","LAMA4",
               "MOG","APBB1IP","PCDH15","VCAN","SPOCK3",
               "ATP13A4","TRPS1","LGI4","FGFR2",
               "PTPRT","PARM1","RUNX2","POU6F2","TRPC3","SPHKAP","ADAMTS17",
               "NKAIN2","CDH9","NDST3","LRRTM4","GRIA1")
ast_selected=intersect(ast_selected,ast_sub_markers$gene)
ast_selected=c("GFAP","S100B",ast_selected)
VlnPlot(ast_sc, features=ast_selected)
# DotPlot(ast_sc, features = ast_selected,cols = c("lightgrey", "brown"))+
#   theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1))
ast_selected_sub=c("GFAP","S100B","SLC1A3","ATP13A4",
                   "RUNX2","PARM1","NKAIN2","CDH9")
VlnPlot(ast_sc, features=ast_selected_sub,ncol=2)
ggsave(filename = "./figures/ast 3 subclass selected markers vlnplot.pdf",width=9.98,height=7.98)
a=VlnPlot(macaque_brain_sc,features = ast_selected_sub,ncol=1 ,group.by = "celltype_6")
ggsave(a,filename = "./figures/macaque_brain_sc ast 3 subclass selected markers vlnplot.pdf",width=9.98,height=20.98)


VlnPlot(ast_sc, intersect(nc_ast_subclass_markers,ast_sub_markers$gene))

VlnPlot(macaque_brain_sc, "RORA",group.by = "celltype_6")
VlnPlot(ast_sc, "RORA")
DotPlot(ast_sc, features = "RORA")

macaque_brain_sc_ast.sub <- subset(macaque_brain_sc@meta.data, celltype=="AST")
macaque_brain_ast_sc <- subset(macaque_brain_sc, cells=row.names(macaque_brain_sc_ast.sub))

macaque_brain_ast_sc@active.assay
ast_sc@active.assay

macaque_brain_ast_sc_for_SingleR <- GetAssayData(macaque_brain_ast_sc, 
                                             slot="data") 
ast_test <- GetAssayData(ast_sc, slot="data")
ast_test.hesc <- SingleR(test=ast_test, 
                            ref=macaque_brain_ast_sc_for_SingleR, 
                            labels=macaque_brain_ast_sc@meta.data$celltype_6, 
                            de.method="wilcox")
table(ast_test.hesc$labels,ast_sc@meta.data$seurat_clusters) 

aa=table(ast_test.hesc$labels,ast_sc@meta.data$seurat_clusters)
aa= apply(aa,2,function(x) x/sum(x))
df=as.data.frame(melt(aa))
df$Var2=as.factor(df$Var2)
g <- ggplot(df, aes(Var2, Var1)) + 
  geom_point(aes(size = value), colour = "blue") + theme_bw()
g


#----------------------------------------------------
# EX
meta1=v1_sc_mfas6@meta.data
id1=grep("Ex",meta1$cellType)
EX_sc <- subset(v1_sc_mfas6, cells=row.names(meta1)[id1])
EX_sc@active.assay
unique(EX_sc@meta.data$cellType)
EX_sc@active.ident
EX_sc
unique(EX_sc@meta.data$cellType)
nc_ex_markers=c("NXPH4","SYT6","FOXP2","FEZF2", "OPRK1","RORB",	"CUX2",	"SLC30A3",	"SLC17A7")
e=DotPlot(EX_sc, features = nc_ex_markers)+
scale_colour_gradient2(low="#2166AC",mid="#F1F4F6",high="#B2182B") +
  scale_size_continuous(range = c(1,8)) +
  labs(x="genes",y="Ex subclass") + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

macaque_brain_sc_EX.sub <- subset(macaque_brain_sc@meta.data, celltype=="EX")
macaque_brain_EX_sc <- subset(macaque_brain_sc, cells=row.names(macaque_brain_sc_EX.sub))

macaque_brain_EX_sc@active.assay
EX_sc@active.assay

macaque_brain_EX_sc_for_SingleR <- GetAssayData(macaque_brain_EX_sc, 
                                                 slot="data") 
EX_test <- GetAssayData(EX_sc, slot="data")
EX_test.hesc <- SingleR(test=EX_test, 
                         ref=macaque_brain_EX_sc_for_SingleR, 
                         labels=macaque_brain_EX_sc@meta.data$celltype_6, 
                         de.method="wilcox")
#EX_sc@meta.data$cellType=factor(EX_sc@meta.data$cellType,levels = paste0("Ex",0:14))
table(EX_test.hesc$labels,EX_sc@meta.data$cellType) 

aa=table(EX_test.hesc$labels,EX_sc@meta.data$cellType)
aa= apply(aa,2,function(x) 100*x/sum(x))
df=as.data.frame(melt(aa))
df=df[df$value!="NaN",]
df$Var2=as.factor(df$Var2)
f <- ggplot(df, aes( Var1,Var2)) + 
  geom_point(aes(size = value), colour = c("brown")) + theme_bw()+
  labs(x="Ex subclass (Lei et.al)",y="Ex subclass") + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

f+e
ggsave(filename = "./figures/EX_subclass_assignment_to_lei et al.pdf",
       width=10.89,height=5.89)

#inh
meta1=v1_sc_mfas6@meta.data
id1=grep("Inh",meta1$cellType)
Inh_sc <- subset(v1_sc_mfas6, cells=row.names(meta1)[id1])
Inh_sc@active.assay
Inh_sc@active.ident

nc_ihn_markers=c("GAD1","ADARB2","LHX6","VIP", "SST","RELN","PVALB","UNC5B","RORA","LAMP5")
h=DotPlot(Inh_sc, features = nc_ihn_markers)+
  scale_colour_gradient2(low="#2166AC",mid="#F1F4F6",high="#B2182B") +
  scale_size_continuous(range = c(1,8)) +
  labs(x="genes",y="Ihn subclass") + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))



macaque_brain_sc_Inh.sub <- subset(macaque_brain_sc@meta.data, celltype=="IN")
macaque_brain_Inh_sc <- subset(macaque_brain_sc, cells=row.names(macaque_brain_sc_Inh.sub))

macaque_brain_Inh_sc@active.assay
Inh_sc@active.assay

macaque_brain_Inh_sc_for_SingleR <- GetAssayData(macaque_brain_Inh_sc, slot="data") 
Inh_test <- GetAssayData(Inh_sc, slot="data")
Inh_test.hesc <- SingleR(test=Inh_test, 
                        ref=macaque_brain_Inh_sc_for_SingleR, 
                        labels=macaque_brain_Inh_sc@meta.data$celltype_6, 
                        de.method="wilcox")
Inh_sc@meta.data$cellType=factor(Inh_sc@meta.data$cellType,
                                 levels = c("Inh0","Inh1","Inh2"))
table(Inh_test.hesc$labels,Inh_sc@meta.data$cellType) 

aa=table(Inh_test.hesc$labels,Inh_sc@meta.data$cellType)
aa= apply(aa,2,function(x) 100*x/sum(x))
df=as.data.frame(melt(aa))
df=df[df$value!="NaN",]
df$Var2=as.factor(df$Var2)
g <- ggplot(df, aes( Var1,Var2)) + 
  geom_point(aes(size = value), colour = c("brown")) + theme_bw()+
  labs(x="Ihn subclass (Lei et.al)",y="Ihn subclass") + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))

g+h
ggsave(filename = "./figures/Ihn_subclass_assignment_to_lei et al.pdf",
       width=10.89,height=5.89)


save.image("v1_snRNAseq_mfas6_and_macaque_brain_snRNAseq_analysis.RData")

save(ast_sc,file="AST_subclass.RData")
DimPlot(ast_sc,reduction = "umap",label = T)+NoLegend()
