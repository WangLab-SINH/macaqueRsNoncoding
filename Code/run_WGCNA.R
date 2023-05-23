
library(reshape2)
library(ggplot2)
library(WGCNA)
library(pheatmap)
library(stringr, quietly = T, warn.conflicts = F)
library(ggplot2)
library(cowplot)
library(ggplotify)


enableWGCNAThreads(8)


# -------------------------------------------------------------------

# load data
load("vsd4.ncx.mat.rmbatch.region.RData")

# make dir to store results
dir.create("./WGCNA") 
setwd("./WGCNA"")


#----build network using expression data of 23613 genes in 100 ncx regions 
# estimate soft power
powers = c(c(1:10), seq(from = 12, to=30, by=2))
type = "signed"
# Choose a set of soft-thresholding powers
sft = pickSoftThreshold(t(vsd4.ncx.mat.rmbatch.region),
                        powerVector = powers,
                        networkType = type,
                        corFnc=bicor,
                        verbose = 5)
pdf("./estimate_softpower_vsd4_ncx_rmbatch_23613genes100regions.pdf",width=9.42,height = 4.59)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = "Scale independence");
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=0.9,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = "Mean connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()
sft$powerEstimate
# build network
mergingThresh=0.25;softPower=12
vsd4_ncx_rmbatch_100region_net = blockwiseModules(t(vsd4.ncx.mat.rmbatch.region),
                                                  corType="bicor",
                                                  maxBlockSize=24000,networkType="signed",
                                                  power=softPower,minModuleSize=50,
                                                  mergeCutHeight=mergingThresh,numericLabels=TRUE,
                                                  saveTOMs=TRUE,pamRespectsDendro=FALSE,nThreads=8,
                                                  saveTOMFileBase="./vsd4_ncx_rmbnatch_23613genes100regions_singed_power12_wgcna")

vsd4_ncx_rmbatch_100region_net_mColors=vsd4_ncx_rmbatch_100region_net$colors
vsd4_ncx_rmbatch_100region_net_mColors = labels2colors(vsd4_ncx_rmbatch_100region_net_mColors)
table(vsd4_ncx_rmbatch_100region_net_mColors) # 16 modules
blocknumber=1
datColors=data.frame(vsd4_ncx_rmbatch_100region_net_mColors)[vsd4_ncx_rmbatch_100region_net$blockGenes[[blocknumber]],]
# plotDendroAndColors
plotDendroAndColors(vsd4_ncx_rmbatch_100region_net$dendrograms[[blocknumber]],
                    colors=datColors,
                    dendroLabels=FALSE,
                    hang =0.03,
                    addGuide=TRUE,
                    guideHang=0.05,
                    groupLabels=c("Module Colors"),
                    addTextGuide=TRUE,
                    main="vsd4_ncx_rmbatch_23613genes100regions_net Cluster Dendeogram")

# merge modules
MEDissThres=0.2 # correlation 0.8
vsd4_ncx_rmbatch_100region_net_merged_net=mergeCloseModules(t(vsd4.ncx.mat.rmbatch.region), 
                                                                 vsd4_ncx_rmbatch_100region_net_mColors,
                                                                 cutHeight = MEDissThres, 
                                                                 verbose = 3)
vsd4_ncx_rmbatch_100region_net_merged_mColors1=vsd4_ncx_rmbatch_100region_net_merged_net$colors;
table(vsd4_ncx_rmbatch_100region_net_merged_mColors1)
vsd4_ncx_rmbatch_100region_net_merged_ME1=vsd4_ncx_rmbatch_100region_net_merged_net$newMEs;

# KME merge modules
vsd4_ncx_rmbatch_100region_net_merged_mColors2=
  moduleMergeUsingKME(t(vsd4.ncx.mat.rmbatch.region),vsd4_ncx_rmbatch_100region_net_merged_mColors1,
  threshPercent=90, mergePercent=40, omitColors=c("grey"))
vsd4_ncx_rmbatch_100region_net_merged_mColors2=vsd4_ncx_rmbatch_100region_net_merged_mColors2$moduleColors;
table(vsd4_ncx_rmbatch_100region_net_merged_mColors2)
plotDendroAndColors(vsd4_ncx_rmbatch_100region_net$dendrograms[[blocknumber]], 
                    cbind(vsd4_ncx_rmbatch_100region_net_mColors, 
                          vsd4_ncx_rmbatch_100region_net_merged_mColors1,
                          vsd4_ncx_rmbatch_100region_net_merged_mColors2),
                    c("Auto", "Merged","KMEmerged"), dendroLabels = FALSE, hang=0.03,
                    main="vsd4_ncx_rmbatch_23613genes100regions_net Cluster Dendeogram ",
                    addGuide = TRUE, guideHang=0.05)


annoRow=data.frame(row.names = ncx_regionOrder,
                   Lobe= ncx_regionOrder_df[ncx_regionOrder,3])
#bk=c(seq(-0.31,0,by=0.01),seq(0.01,0.83,by=0.01))
removeM=which(colnames(vsd4_ncx_rmbatch_100region_net_merged_ME1)=="MEgrey")

pheatmap(t(vsd4_ncx_rmbatch_100region_net_merged_ME1[ncx_regionOrder,-removeM]),
         cluster_cols = FALSE,
         #color = c(colorRampPalette(colors = c("blue","lightyellow"))(40),colorRampPalette(colors = c("lightyellow","red"))(80)),
         #breaks = bk,
         cellwidth = 10, cellheight = 10, 
         color = colorRampPalette(colors = c("blue","lightyellow","red"))(100),
         annotation_col =annoRow)
