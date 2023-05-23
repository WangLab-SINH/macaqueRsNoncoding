library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggnewscale)
library(ggplot2)
library(stingr)



# define function
myFisher=function(mGenes,fGenes,allgs){
  df1<-intersect(mGenes,fGenes)
  df2<-intersect(setdiff(allgs,mGenes),fGenes)
  df3<-setdiff(mGenes,fGenes)
  df4<-length(setdiff(allgs,union(union(df1,df2),df3)))
  df<-matrix(c(length(df1),length(df2),length(df3),df4),ncol = 2)
  ft<-fisher.test(df)
  
  odds=ft$estimate
  pp=ft$p.value
  res=c(length(mGenes),length(fGenes), length(df1),odds, pp)
  names(res)=c("set1","set2","intersect","oddsratio","p")
  return(res)
}



# load crab expressed genes
load("exprSet4_ncx_phenSet4_757samples.RData")
load("vsd4.ncx.mat.rmbatch.RData")
keep=apply(exprSet4.ncx, 1, function(x){
  ratio=sum(x>=1)
  ratio=ratio/length(x)
  state=F
  if(ratio>=0.1){
    state=T
  }
  return(state)
})

crab_expressedGenes=row.names(exprSet4.ncx)[keep]


# human and crab ortholog
human_crab_ortholog=read.csv("human_crabeating_ortholog.txt", sep="\t", header=T)
human_crab_ortholog=human_crab_ortholog[human_crab_ortholog$Crab.eating.macaque.gene.stable.ID!="",]

crab_expressedGenes_human=human_crab_ortholog[match(crab_expressedGenes,human_crab_ortholog$Crab.eating.macaque.gene.stable.ID ),]
crab_expressedGenes_human=na.omit(crab_expressedGenes_human) 
crab_expressedGenes_human=unique( crab_expressedGenes_human$Gene.name ) 

# crab functional genes
ww_crab_functional_genes=read.csv("funcGenes.csv", header=T)
ww_crab_functional_coding=ww_crab_functional_genes[ww_crab_functional_genes[,2]=="coding",1]
ww_crab_functional_noncoding=ww_crab_functional_genes[ww_crab_functional_genes[,2]=="non-coding",1]
 
ww_crab_functional_coding_human=na.omit( human_crab_ortholog[match(ww_crab_functional_coding, human_crab_ortholog$Crab.eating.macaque.gene.stable.ID ), ] )  # 49
ww_crab_functional_noncoding_human=na.omit( human_crab_ortholog[match(ww_crab_functional_noncoding, human_crab_ortholog$Crab.eating.macaque.gene.stable.ID ), ] ) # 37

ww_crab_functional_coding_human=unique(ww_crab_functional_coding_human$Gene.name) 
ww_crab_functional_noncoding_human=unique(ww_crab_functional_noncoding_human$Gene.name)



#------- calculate crab fMRI genes related genes 
vsd4.ncx.mat.rmbatch=vsd4.ncx.mat.rmbatch[crab_expressedGenes,]

ww_crab_functional_coding_pccRESlist=lapply(ww_crab_functional_coding, function(x){
  v1=as.numeric(vsd4.ncx.mat.rmbatch[x,])
  tmp=vsd4.ncx.mat.rmbatch[setdiff(crab_expressedGenes,x),]
  tmpres=t(apply(tmp,1,function(y){
    corRes=cor.test(v1,y,methods="pearson")
    return(c(corRes$estimate, corRes$p.value))
  }) )
  tmpres=as.data.frame(tmpres)
  tmpres$FDR=p.adjust(tmpres[,2],"fdr")
  tmpres$gene1=rep(x,nrow(tmpres))
  tmpres$gene2=row.names(tmpres)
  tmpres=tmpres[,c(4,5,1,2,3)]
  colnames(tmpres)=c("Gene1","Gene2","R","P","FDR")
  return(tmpres)
})

ww_crab_functional_coding_pccSigGenesNum=c()
for(i in seq(0,1,by=0.1)){
  ww_crab_functional_coding_pccRESlist_subset=lapply(ww_crab_functional_coding_pccRESlist,function(x){
    res=subset(x,FDR<0.05 & (R>=i | R<=(-i)))
  })
  ww_crab_functional_coding_pccSigGenesNum=c(ww_crab_functional_coding_pccSigGenesNum,
                                             length(Reduce(union,(lapply(ww_crab_functional_coding_pccRESlist_subset,function(x) return(x[,2])))))
  )
}
dat=data.frame(cbind(FDR=seq(0,1,by=0.1), 
                     GeneNum=ww_crab_functional_coding_pccSigGenesNum))

ggplot(dat, aes(x=FDR,y=GeneNum))+geom_bar(stat="identity")+scale_x_continuous(breaks = seq(0,1,by=0.1))
  



ww_crab_functional_coding_pccRESlist_subset=lapply(ww_crab_functional_coding_pccRESlist,function(x){
  res=subset(x,FDR<0.05 & (R>=0.4 | R<=(-0.4)))
})
length(Reduce(union,(lapply(ww_crab_functional_coding_pccRESlist_subset,function(x) return(x[,2]))))) # 5004
ww_crab_functional_coding_SigGenes=Reduce(union,(lapply(ww_crab_functional_coding_pccRESlist_subset,function(x) return(x[,2])))) # 5004
ww_crab_functional_coding_PosGenes=Reduce(union,(lapply(ww_crab_functional_coding_pccRESlist_subset,function(x) return(x[x[,3]>=0.4,2])))) # 2809
ww_crab_functional_coding_NegGenes=Reduce(union,(lapply(ww_crab_functional_coding_pccRESlist_subset,function(x) return(x[x[,3]<=(-0.4),2])))) # 2214



ww_crab_functional_coding_PosGenes_human=na.omit(human_crab_ortholog$Gene.name[
  match(ww_crab_functional_coding_PosGenes,human_crab_ortholog$Crab.eating.macaque.gene.stable.ID)]) # 2658
ww_crab_functional_coding_NegGenes_human=na.omit(human_crab_ortholog$Gene.name[
  match(ww_crab_functional_coding_NegGenes,human_crab_ortholog$Crab.eating.macaque.gene.stable.ID)]) # 2120

intersect(ww_crab_functional_coding_PosGenes_human, ww_crab_functional_coding_NegGenes_human) # 20


ww_crab_functional_coding_PosGenes_human_GOBP=enrichGO(gene=ww_crab_functional_coding_PosGenes_human,
                                                     OrgDb = "org.Hs.eg.db",
                                                     keyType = "SYMBOL",
                                                     ont = "BP",
                                                     universe = crab_expressedGenes_human,
                                                     pAdjustMethod = "BH",
                                                     pvalueCutoff = 0.05,
                                                     qvalueCutoff = 0.05)
ww_crab_functional_coding_PosGenes_human_GOBP=simplify(ww_crab_functional_coding_PosGenes_human_GOBP)
barplot(ww_crab_functional_coding_PosGenes_human_GOBP)
x2 <- pairwise_termsim(ww_crab_functional_coding_PosGenes_human_GOBP)
treeplot(x2, showCategory = 60)


ww_crab_functional_coding_NegGenes_human_GOBP=enrichGO(gene=ww_crab_functional_coding_NegGenes_human,
                                                       OrgDb = "org.Hs.eg.db",
                                                       keyType = "SYMBOL",
                                                       ont = "BP",
                                                       universe = crab_expressedGenes_human,
                                                       pAdjustMethod = "BH",
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 0.05)
ww_crab_functional_coding_NegGenes_human_GOBP=simplify(ww_crab_functional_coding_NegGenes_human_GOBP)
x2 <- pairwise_termsim(ww_crab_functional_coding_NegGenes_human_GOBP)
treeplot(x2, showCategory = 60)
#emapplot(x2, showCategory = 60)



crab_expressedGenes_human_entrez=bitr(crab_expressedGenes_human,
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID",
                                      OrgDb = "org.Hs.eg.db")
crab_expressedGenes_human_entrez = crab_expressedGenes_human_entrez[,2] # 17826

ww_crab_functional_coding_PosGenes_human_entrez=bitr(ww_crab_functional_coding_PosGenes_human,
                                                     fromType = "SYMBOL",
                                                     toType = "ENTREZID",
                                                     OrgDb = "org.Hs.eg.db")
ww_crab_functional_coding_PosGenes_human_entrez=unique( ww_crab_functional_coding_PosGenes_human_entrez$ENTREZID) # 2614


ww_crab_functional_coding_NegGenes_human_entrez=bitr(ww_crab_functional_coding_NegGenes_human,
                                                     fromType = "SYMBOL",
                                                     toType = "ENTREZID",
                                                     OrgDb = "org.Hs.eg.db")
ww_crab_functional_coding_NegGenes_human_entrez=unique( ww_crab_functional_coding_NegGenes_human_entrez$ENTREZID) # 2089





ww_crab_functional_coding_PosGenes_human_KEGG=enrichKEGG(gene=ww_crab_functional_coding_PosGenes_human_entrez,
                                                       organism = "hsa",
                                                       keyType = "kegg",
                                                       universe = crab_expressedGenes_human_entrez,
                                                       pAdjustMethod = "BH",
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 0.05)
edox2 <- pairwise_termsim(ww_crab_functional_coding_PosGenes_human_KEGG)
treeplot(edox2,showCategory = 60)


ww_crab_functional_coding_NegGenes_human_KEGG=enrichKEGG(gene=ww_crab_functional_coding_NegGenes_human_entrez,
                                                         organism = "hsa",
                                                         keyType = "kegg",
                                                         universe = crab_expressedGenes_human_entrez,
                                                         pAdjustMethod = "BH",
                                                         pvalueCutoff = 0.05,
                                                         qvalueCutoff = 0.05)
barplot(ww_crab_functional_coding_NegGenes_human_KEGG)
edox2 <- pairwise_termsim(ww_crab_functional_coding_NegGenes_human_KEGG)
treeplot(edox2)



#---- noncoding correlated genes
ww_crab_functional_noncoding_pccRESlist=lapply(ww_crab_functional_noncoding, function(x){
  v1=as.numeric(vsd4.ncx.mat.rmbatch[x,])
  tmp=vsd4.ncx.mat.rmbatch[setdiff(crab_expressedGenes,x),]
  tmpres=t(apply(tmp,1,function(y){
    corRes=cor.test(v1,y,methods="pearson")
    return(c(corRes$estimate, corRes$p.value))
  }) )
  tmpres=as.data.frame(tmpres)
  tmpres$FDR=p.adjust(tmpres[,2],"fdr")
  tmpres$gene1=rep(x,nrow(tmpres))
  tmpres$gene2=row.names(tmpres)
  tmpres=tmpres[,c(4,5,1,2,3)]
  colnames(tmpres)=c("Gene1","Gene2","R","P","FDR")
  return(tmpres)
})

ww_crab_functional_noncoding_pccSigGenesNum=c()
for(i in seq(0,1,by=0.1)){
  ww_crab_functional_noncoding_pccRESlist_subset=lapply(ww_crab_functional_noncoding_pccRESlist,function(x){
    res=subset(x,FDR<0.05 & (R>=i | R<=(-i)))
  })
  ww_crab_functional_noncoding_pccSigGenesNum=c(ww_crab_functional_noncoding_pccSigGenesNum,
                                             length(Reduce(union,(lapply(ww_crab_functional_noncoding_pccRESlist_subset,function(x) return(x[,2])))))
  )
}
dat=data.frame(cbind(FDR=seq(0,1,by=0.1), 
                     GeneNum=ww_crab_functional_noncoding_pccSigGenesNum))
ggplot(dat, aes(x=FDR,y=GeneNum))+geom_bar(stat="identity")+scale_x_continuous(breaks = seq(0,1,by=0.1))





ww_crab_functional_noncoding_pccRESlist_subset=lapply(ww_crab_functional_noncoding_pccRESlist,function(x){
  res=subset(x,FDR<0.05 & (R>=0.5 | R<=(-0.5)))
})
length(Reduce(union,(lapply(ww_crab_functional_noncoding_pccRESlist_subset,function(x) return(x[,2]))))) # 4532

ww_crab_functional_noncoding_SigGenes=Reduce(union,(lapply(ww_crab_functional_noncoding_pccRESlist_subset,function(x) return(x[,2])))) # 4532
ww_crab_functional_noncoding_PosGenes=Reduce(union,(lapply(ww_crab_functional_noncoding_pccRESlist_subset,function(x) return(x[x[,3]>=0.5,2])))) # 3442
ww_crab_functional_noncoding_NegGenes=Reduce(union,(lapply(ww_crab_functional_noncoding_pccRESlist_subset,function(x) return(x[x[,3]<=(-0.5),2])))) # 2078


ww_crab_functional_noncoding_PosGenes_human=na.omit(human_crab_ortholog$Gene.name[
  match(ww_crab_functional_noncoding_PosGenes,human_crab_ortholog$Crab.eating.macaque.gene.stable.ID)]) # 3127
ww_crab_functional_noncoding_NegGenes_human=na.omit(human_crab_ortholog$Gene.name[
  match(ww_crab_functional_noncoding_NegGenes,human_crab_ortholog$Crab.eating.macaque.gene.stable.ID)]) # 1931]

intersect(ww_crab_functional_noncoding_PosGenes_human, ww_crab_functional_noncoding_NegGenes_human) # 910


length(intersect(union(ww_crab_functional_noncoding_PosGenes,ww_crab_functional_noncoding_NegGenes ),
                 union(ww_crab_functional_coding_PosGenes,ww_crab_functional_coding_NegGenes))) # 2300

ww_crab_functional_noncoding_PCCGenes=union(ww_crab_functional_noncoding_PosGenes,ww_crab_functional_noncoding_NegGenes)
ww_crab_functional_coding_PCCGenes=union(ww_crab_functional_coding_PosGenes,ww_crab_functional_coding_NegGenes)

venn(list(CodingPCCGenes=ww_crab_functional_coding_PCCGenes, noncodingPCCGenes=ww_crab_functional_noncoding_PCCGenes),ilabels = T,zcolor = "style")
myFisher(ww_crab_functional_coding_PCCGenes,ww_crab_functional_noncoding_PCCGenes, crab_expressedGenes )




venn(list(CodingPCCGenes=intersect(ww_crab_functional_coding_PCCGenes,mfas5PCG), 
          noncodingPCCGenes=intersect(ww_crab_functional_noncoding_PCCGenes,mfas5PCG)),ilabels = T,zcolor = "style")

myFisher(intersect(ww_crab_functional_coding_PCCGenes,mfas5PCG),
         intersect(ww_crab_functional_noncoding_PCCGenes,mfas5PCG), 
         intersect(crab_expressedGenes,mfas5PCG) )


# functional analysis of noncoding correlated genes
ww_crab_functional_noncoding_PosGenes_human_GOBP=enrichGO(gene=ww_crab_functional_noncoding_PosGenes_human,
                                                       OrgDb = "org.Hs.eg.db",
                                                       keyType = "SYMBOL",
                                                       ont = "BP",
                                                       universe = crab_expressedGenes_human,
                                                       pAdjustMethod = "BH",
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 0.05)
ww_crab_functional_noncoding_PosGenes_human_GOBP=simplify(ww_crab_functional_noncoding_PosGenes_human_GOBP)
barplot(ww_crab_functional_noncoding_PosGenes_human_GOBP)
x2 <- pairwise_termsim(ww_crab_functional_noncoding_PosGenes_human_GOBP)
treeplot(x2,showCategory = 60)


ww_crab_functional_noncoding_NegGenes_human_GOBP=enrichGO(gene=ww_crab_functional_noncoding_NegGenes_human,
                                                       OrgDb = "org.Hs.eg.db",
                                                       keyType = "SYMBOL",
                                                       ont = "BP",
                                                       universe = crab_expressedGenes_human,
                                                       pAdjustMethod = "BH",
                                                       pvalueCutoff = 0.05,
                                                       qvalueCutoff = 0.05)
ww_crab_functional_noncoding_NegGenes_human_GOBP=simplify(ww_crab_functional_noncoding_NegGenes_human_GOBP)
x2 <- pairwise_termsim(ww_crab_functional_noncoding_NegGenes_human_GOBP)
treeplot(x2, showCategory = 60)




ww_crab_functional_noncoding_PosGenes_human_entrez=bitr(ww_crab_functional_noncoding_PosGenes_human,
                                                     fromType = "SYMBOL",
                                                     toType = "ENTREZID",
                                                     OrgDb = "org.Hs.eg.db")
ww_crab_functional_noncoding_PosGenes_human_entrez=unique( ww_crab_functional_noncoding_PosGenes_human_entrez$ENTREZID) # 3069


ww_crab_functional_noncoding_NegGenes_human_entrez=bitr(ww_crab_functional_noncoding_NegGenes_human,
                                                     fromType = "SYMBOL",
                                                     toType = "ENTREZID",
                                                     OrgDb = "org.Hs.eg.db")
ww_crab_functional_noncoding_NegGenes_human_entrez=unique( ww_crab_functional_noncoding_NegGenes_human_entrez$ENTREZID) # 1893





ww_crab_functional_noncoding_PosGenes_human_KEGG=enrichKEGG(gene=ww_crab_functional_noncoding_PosGenes_human_entrez,
                                                         organism = "hsa",
                                                         keyType = "kegg",
                                                         universe = crab_expressedGenes_human_entrez,
                                                         pAdjustMethod = "BH",
                                                         pvalueCutoff = 0.05,
                                                         qvalueCutoff = 0.05)
edox2 <- pairwise_termsim(ww_crab_functional_noncoding_PosGenes_human_KEGG)
treeplot(edox2,showCategory = 30)


ww_crab_functional_noncoding_NegGenes_human_KEGG=enrichKEGG(gene=ww_crab_functional_noncoding_NegGenes_human_entrez,
                                                         organism = "hsa",
                                                         keyType = "kegg",
                                                         universe = crab_expressedGenes_human_entrez,
                                                         pAdjustMethod = "BH",
                                                         pvalueCutoff = 0.05,
                                                         qvalueCutoff = 0.05)
edox2 <- pairwise_termsim(ww_crab_functional_noncoding_NegGenes_human_KEGG)
treeplot(edox2,showCategory = 30)