#Full R script for all anaylses performed in the manuscript, 
#"24-hour age difference causes twice as much gene expression divergence as 100 generations of adaptation to a new novel environment", 
#submitted to genes
##Sheng-Kai Hsu 2018.12.27
##Sheng-Kai Hsu 2019.01.19 updated

rm(list=ls())
####step0: import required libraries and customized function####
library(edgeR)
Library(DESeq2)
library(ExpressionNormalizationWorkflow)
library(VennDiagram)
library(RColorBrewer)
library(ggplot2)
library(pheatmap)
library(topGO)
cont_table=function(query,background,classifyer){
  p1=length(intersect(query,classifyer))
  q1=length(query)-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}
ID_converter=function(ID,db,attributes,filters){
  getBM(attributes=attributes,filters=filters,mart=db,values=ID)
}
setwd("Your working directory")

####step1: input RNASeq count table####
dat_use=read.csv("./Fl_D.sim.RNASeq_readcounts.csv",stringsAsFactors = F,row.names = 1)
counts_use=dat_use[apply(cpm(dat_use),1,function(x) !sum(x<1)>=1),]#filtered for lowly expressed genes
evo=substr(colnames(dat_use),1,1)#translate labels to evolutionary states (B:ancestral; H: evolved)
age=ifelse(substr(colnames(dat_use),9,9)==1,6,5)#translate labels to age (1: 6-d-old; else: 5-d-old)
y=DGEList(counts=counts_use,group = evo)
y=calcNormFactors(y)
####step2: variance decomposition for gene expression####
pca=prcomp(t(log(cpm(y))))
ve=round(pca$sdev^2/sum(pca$sdev^2),4)
mycol_pca=c()
mycol_pca[evo%in%"B"]="forestgreen"
mycol_pca[evo%in%"H"&age==5]="salmon"
mycol_pca[evo%in%"H"&age==6]="maroon"

meta.table=data.frame(evo=evo,age=age,row.names = colnames(y))
annot=data.frame(labelDescription=c("Factor levels","Factor levels"))
annot_factors=AnnotatedDataFrame(meta.table,annot)
expr.set=ExpressionSet(log(cpm(y)),annot_factors)
pvca_res=pvcAnaly(expr.set, 0.75,c("age","evo"))

####Step3: linear modeling and DE analysis####
group=paste0(evo,age)
ModelDesign=model.matrix(~0+group)
DGE=estimateDisp(y,design = ModelDesign,robust = T)
GLM=glmFit(DGE,design = ModelDesign)
mycontrast=makeContrasts("evolution"=groupH5-groupB5,"age"=groupH6-groupH5,"mismatched"=groupH6-groupB5,levels = ModelDesign)
LRT_res_evo=glmLRT(GLM,contrast = mycontrast[,"evolution"])
res_table_evo=LRT_res_evo$table
res_table_evo$padj=p.adjust(res_table_evo$PValue,method = "BH")

LRT_res_age=glmLRT(GLM,contrast = mycontrast[,"age"])
res_table_age=LRT_res_age$table
res_table_age$padj=p.adjust(res_table_age$PValue,method = "BH")

LRT_res_mis=glmLRT(GLM,contrast = mycontrast[,"mismatched"])
res_table_mis=LRT_res_mis$table
res_table_mis$padj=p.adjust(res_table_mis$PValue,method = "BH")

####Step4: identification and comparison of DE genes among different contrasts####
query_ID=list(age_up=rownames(y)[res_table_age$padj<0.05&res_table_age$logFC>log2(1.25)],
                  age_dn=rownames(y)[res_table_age$padj<0.05&res_table_age$logFC< -log2(1.25)],
                  evo_up=rownames(y)[res_table_evo$padj<0.05&res_table_evo$logFC>log2(1.25)],
                  evo_dn=rownames(y)[res_table_evo$padj<0.05&res_table_evo$logFC< -log2(1.25)],
                  false_evo_up=rownames(y)[res_table_mis$padj<0.05&res_table_mis$logFC>log2(1.25)],
                  false_evo_dn=rownames(y)[res_table_mis$padj<0.05&res_table_mis$logFC< -log2(1.25)])
query_ID_cat=list(age_s_up=setdiff(query_ID[[1]],unlist(query_ID[1:4][-1])),
                  age_s_dn=setdiff(query_ID[[2]],unlist(query_ID[1:4][-2])),
                  evo_s_up=setdiff(query_ID[[3]],unlist(query_ID[1:4][-3])),
                  evo_s_dn=setdiff(query_ID[[4]],unlist(query_ID[1:4][-4])),
                  exaggerated_up=rownames(y)[res_table_age$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC>log2(1.25)&res_table_age$logFC>log2(1.25)],
                  exaggerated_dn=rownames(y)[res_table_age$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC< -log2(1.25)&res_table_age$logFC< -log2(1.25)],
                  diminished_up=rownames(y)[res_table_age$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC>log2(1.25)&res_table_age$logFC< -log2(1.25)],
                  diminished_dn=rownames(y)[res_table_age$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC< -log2(1.25)&res_table_age$logFC>log2(1.25)]
)
query_ID_cat2=list(FP_up=setdiff(query_ID[[5]],unlist(query_ID[-1:-2][-3])),
                   FP_dn=setdiff(query_ID[[6]],unlist(query_ID[-1:-2][-4])),
                   FN_up=setdiff(query_ID[[3]],unlist(query_ID[-1:-2][-1])),
                   FN_dn=setdiff(query_ID[[4]],unlist(query_ID[-1:-2][-2])),
                   consistent_up=rownames(y)[res_table_mis$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC>log2(1.25)&res_table_mis$logFC>log2(1.25)],
                   consistent_dn=rownames(y)[res_table_mis$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC< -log2(1.25)&res_table_mis$logFC< -log2(1.25)],
                   reversed_up=rownames(y)[res_table_mis$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC>log2(1.25)&res_table_mis$logFC< -log2(1.25)],
                   reversed_dn=rownames(y)[res_table_mis$padj<0.05&res_table_evo$padj<0.05&res_table_evo$logFC< -log2(1.25)&res_table_mis$logFC>log2(1.25)]
)
new_list=list(age_specific=unlist(query_ID_cat[1:2]),evolution_specific=unlist(query_ID_cat[3:4]),
              exaggerated=unlist(query_ID_cat[5:6]),diminished=unlist(query_ID_cat[7:8]),
              n.s.=rownames(y)[!rownames(y)%in%unlist(query_ID[1:4])])
new_list2=list("False Positive"=unlist(query_ID_cat2[1:2]),"False Negative"=unlist(query_ID_cat2[3:4]),
               Consistent=unlist(query_ID_cat2[5:6]),Reversed=unlist(query_ID_cat2[7:8]))
pie_tab=sapply(new_list2,function(x) sapply(new_list,function(y) sum(x%in%y)))


####Step5: GO analysis####
GO_res_table=list()
for (i in names(query_ID)){
  tmp=factor(as.integer(rownames(y)%in%unlist(query_ID[[i]])))
  names(tmp)=rownames(y)
  tgd1=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd1, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd1, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd1,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  GO_res_table[[i]]=tmp_res
}
GO_res_table=lapply(GO_res_table,function(x) {
  x$Fisher.weight01[x$Fisher.weight01=="< 1e-30"]=1e-30
  return(x)})

GO_res_table2=list()
for (i in names(query_ID_cat)){
  tmp=factor(as.integer(rownames(y)%in%unlist(query_ID_cat[[i]])))
  names(tmp)=rownames(y)
  tgd2=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd2, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd2, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd2,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  GO_res_table2[[i]]=tmp_res
}
GO_res_table2=lapply(GO_res_table2,function(x) {
  x$Fisher.weight01[x$Fisher.weight01=="< 1e-30"]=1e-30
  return(x)})

GO_res_table3=list()
for (i in names(query_ID_cat2)[-7]){ #-7 because error comes when the query is empty
  tmp=factor(as.integer(rownames(y)%in%unlist(query_ID_cat2[[i]])))
  names(tmp)=rownames(y)
  tgd2=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
  resTopGO.classic=runTest(tgd2, algorithm = "classic", statistic = "Fisher")
  resTopGO.weight01=runTest(tgd2, algorithm = "weight01", statistic = "Fisher")
  tmp_res=GenTable(tgd2,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
  #  tmp_res$Fisher.classic=resTopGO.classic@score
  GO_res_table3[[i]]=tmp_res
}
GO_res_table3=lapply(GO_res_table3,function(x) {
  x$Fisher.weight01[x$Fisher.weight01=="< 1e-30"]=1e-30
  return(x)})

####step6: identification of genes of interest based on GO analysis####
GOI_set1=intersect(genesInTerm(tgd2,"GO:0008406")[[1]],query_ID_cat2$FP_up)
GOI_expr_set1=apply(apply(cpm(y)[GOI_set1,],1,function(x) scale(log10(x))),2,function(y) tapply(y,group,mean))

GOI_set2=intersect(genesInTerm(tgd2,"GO:0006270")[[1]],query_ID_cat2$FP_dn)
GOI_expr_set2=apply(apply(cpm(y)[GOI_set2,],1,function(x) scale(log10(x))),2,function(y) tapply(y,group,mean))

GOI_set3=intersect(genesInTerm(tgd2,"GO:0050877")[[1]],query_ID_cat2$FN_dn)
GOI_expr_set3=apply(apply(cpm(y)[GOI_set3,],1,function(x) scale(log10(x))),2,function(y) tapply(y,group,mean))

GOI_set4=intersect(genesInTerm(tgd2,"GO:0070983")[[1]],query_ID_cat2$reversed_dn)
GOI_expr_set4=apply(apply(cpm(y)[GOI_set4,],1,function(x) scale(log10(x))),2,function(y) tapply(y,group,mean))

####step7: PCA analysis on 374 genes which decreased during evolution but up-regulated as flies aged####
pca_goi=prcomp(t(log(cpm(y)[unlist(query_ID_cat[8]),])))
ve_goi=round(pca_goi$sdev^2/sum(pca_goi$sdev^2),4)

####Step8: comparison to Carlson et al., 2015####
carlson_tab=read.csv("./carlson_2015.csv",header = T,stringsAsFactors = F,na.strings = "")#Table S1 Carlson et al. 2015
carlson_goi=carlson_tab$FBGN
fisher.test(cont_table(unlist(query_ID[1:2]),rownames(y),carlson_goi),alternative="greater")

tmp=factor(as.integer(rownames(y)%in%carlson_goi))
names(tmp)=rownames(y)
tgd5=new( "topGOdata", ontology="BP", allGenes = tmp, nodeSize=5,annot=annFUN.org, mapping="org.Dm.eg.db", ID = "ensembl")
resTopGO.classic=runTest(tgd5, algorithm = "classic", statistic = "Fisher")
resTopGO.weight01=runTest(tgd5, algorithm = "weight01", statistic = "Fisher")
GO_res_table5=GenTable(tgd5,Fisher.classic = resTopGO.classic,Fisher.weight01=resTopGO.weight01,orderBy = "Fisher.weight01",ranksOf="Fisher.classic",topNodes=length(resTopGO.classic@score),numChar=100)
fisher.test(cont_table(c(GO_res_table$age_up$GO.ID[GO_res_table$age_up$Fisher.weight01<0.05],
                         GO_res_table$age_dn$GO.ID[GO_res_table$age_dn$Fisher.weight01<0.05]),
                       GO_res_table$age_up$GO.ID,GO_res_table5$GO.ID[GO_res_table5$Fisher.weight01<0.05]),alternative = "greater")


####Step9: Comparison to DESeq2####
coldata=data.frame(condition=group,type=rep("single-read",length(group)))
rownames(coldata)=colnames(y_mon)
dds=DESeqDataSetFromMatrix(countData = monster_counts_use,colData = coldata,design = ~condition)
test=DESeq(dds)
res_DESeq=results(test,contrast = c("condition","H5","B5"))#evolution contrast
res_DESeq2=results(test,contrast = c("condition","H6","H5"))#age contrast
res_DESeq3=results(test,contrast = c("condition","H6","B5"))#mismatched contrast

####visualization####
#Figure 2
png("./figure2.png",width = 24,height = 16,units = "cm",res = 600)
par(mfrow=c(1,2))
plot(pca$x,xlab=paste0("PC1 (", ve[1]*100,"%)"),ylab=paste0("PC2 (", ve[2]*100,"%)"),
     col=mycol_pca,pch=19,asp=1,main="Principal Component Analysis")
legend("topright",pch = 19,col=c("forestgreen","salmon","maroon"),legend = c("anc. 5-day-old","evo. 5-day-old","evo. 6-day-old"))
bp=barplot(pvca_res$dat,  xlab = "Effects",
           ylab = "Weighted average proportion variance", ylim= c(0,1.1),
           col = "grey20", las=2, main="Principal Variance Component Analysis")
axis(1, at = bp, labels = pvca_res$label, xlab = "Effects", cex.axis = 1, las=1)
text(bp,pvca_res$dat,round(pvca_res$dat,3),pos = 3,cex = 0.75)
dev.off()

#Figure 3A
venn.diagram(list(unlist(query_ID_mon[1:2]),unlist(query_ID_mon[3:4]),unlist(query_ID_mon[5:6])),
             "./Figure3A.png",
             width = 12,height = 12,units = "cm",res=600,imagetype = "png",cat.cex=1.25,cex=1.5,
             category.names = c("age","evolution","mismatched"),col=c("forestgreen","salmon","maroon"))

#Figure 3B
png("./figure3AB.png",width = 16,height = 16,units = "cm",res=600)
plot(res_table_evo$logFC,res_table_age$logFC,asp=1,xlab=expression(paste("contrast: ",italic("evolution"))),ylab=expression(paste("contrast: ",italic("age"))),
     pch=19,col="grey",cex=1,cex.lab=1.4,xlim=c(-4,4),ylim=c(-4,4))
abline(v=0,h=0,col="grey70",lty=2)
abline(a=0,b=1,col="grey70",lty=2)
abline(a=0,b=-1,col="grey70",lty=2)
points(res_table_evo$logFC[(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))],
       res_table_age$logFC[(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))],col=brewer.pal(12,"Paired")[3],pch=19,cex=1)
points(res_table_evo$logFC[(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))],
       res_table_age$logFC[(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))],col=brewer.pal(12,"Paired")[5],pch=19,cex=1)
points(res_table_evo$logFC[(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_age$logFC*res_table_evo$logFC>0],
       res_table_age$logFC[(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_age$logFC*res_table_evo$logFC>0],col=brewer.pal(12,"Paired")[4],pch=19,cex=1)
points(res_table_evo$logFC[(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_age$logFC*res_table_evo$logFC<0],
       res_table_age$logFC[(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_age$logFC*res_table_evo$logFC<0],col=brewer.pal(12,"Paired")[8],pch=19,cex=1)
legend("topleft",legend = c("age-specific: 1617","evolution-specific: 500","exaggerated: 38","diminished: 396"),col = brewer.pal(12,"Paired")[c(3,5,4,8)],pch=19,cex = 1.25)
dev.off()

#Figure 3C
png("./figure3C.png",width = 16,height = 16,units = "cm",res=600)
plot(res_table_evo$logFC,res_table_evo2$logFC,asp=1,xlab=expression(paste("contrast: ",italic("evolution"))),ylab=expression(paste("contrast: ",italic("mismatched"))),
     pch=19,col="grey",cex=1,cex.lab=1.4,xlim=c(-4,4),ylim=c(-4,4))
abline(v=0,h=0,col="grey70",lty=2)
abline(a=0,b=1,col="grey70",lty=2)
abline(a=0,b=-1,col="grey70",lty=2)
points(res_table_evo$logFC[(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))],
       res_table_evo2$logFC[(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))],col=brewer.pal(12,"Paired")[9],pch=19,cex=1)
points(res_table_evo$logFC[(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))],
       res_table_evo2$logFC[(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))],col=brewer.pal(12,"Paired")[5],pch=19,cex=1)
points(res_table_evo$logFC[(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_evo2$logFC*res_table_evo$logFC>0],
       res_table_evo2$logFC[(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_evo2$logFC*res_table_evo$logFC>0],col=brewer.pal(12,"Paired")[6],pch=19,cex=1)
points(res_table_evo$logFC[(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_evo2$logFC*res_table_evo$logFC<0],
       res_table_evo2$logFC[(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))&(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))&res_table_evo2$logFC*res_table_evo$logFC<0],col=brewer.pal(12,"Paired")[10],pch=19,cex=1)
legend("topleft",legend = c("false positive: 1260","false negative: 643","consistent: 278","reversed: 13"),col = brewer.pal(12,"Paired")[c(9,5,6,10)],pch=19,cex = 1.25)
dev.off()


#Figure 4
png("./figure4.png",height = 24,width = 24,units = "cm",res = 600,pointsize = 12)
par(mfrow=c(2,2))
boxplot(as.vector(GOI_expr_set1)~rep(sort(unique(group)),length(GOI_set1)),names=c("Anc. 5d-old","Evo. 5d-old","Evo. 6d-old"),
        border=c("forestgreen","salmon","maroon"),lwd=1.5,main="GO:0008406 Gonad development (16)",ylab="Normalized expression",cex.axis=1)
boxplot(as.vector(GOI_expr_set2)~rep(sort(unique(group)),length(GOI_set2)),names=c("Anc. 5d-old","Evo. 5d-old","Evo. 6d-old"),
        border=c("forestgreen","salmon","maroon"),lwd=1.5,main="GO:0006270 DNA replication initiation (11)",ylab="Normalized expression",cex.axis=1)
boxplot(as.vector(GOI_expr_set3)~rep(sort(unique(group)),length(GOI_set3)),names=c("Anc. 5d-old","Evo. 5d-old","Evo. 6d-old"),
        border=c("forestgreen","salmon","maroon"),lwd=1.5,main="GO:0050877 Nervous system process (59)",ylab="Normalized expression",cex.axis=1)
boxplot(as.vector(GOI_expr_set4)~rep(sort(unique(group)),length(GOI_set4)),names=c("Anc. 5d-old","Evo. 5d-old","Evo. 6d-old"),
        border=c("forestgreen","salmon","maroon"),lwd=1.5,main="GO:0070983 Dendrite guidance (3)",ylab="Normalized expression",cex.axis=1)
dev.off()

#Figure 5
png("./figure5.png",height = 16,width = 16,units = "cm",res=600)
plot(pca_goi$x,col=mycol_pca,pch=19,asp=1,xlab = paste0("PC1 (",ve_goi[1]*100,"%)"),ylab = paste0("PC2 (",ve_goi[2]*100,"%)"))
legend("topright",pch = 19,col=c("forestgreen","salmon","maroon"),legend = c("anc. 5-day-old","evo. 5-day-old","evo. 6-day-old"))
dev.off()
anno_col=data.frame(pop=group)
rownames(anno_col)=colnames(y)
png("/Volumes/Temp1/shengkai/ageing_CGE/ageing_result/heatmap_374_diminished_genes.png",height = 16,width = 16,units = "cm",res=600)
pheatmap(log(cpm(y)[unlist(query_ID_cat[8]),]),scale = "row",show_rownames = F,show_colnames = F,annotation_col = anno_col)
dev.off()

#Figure S1
venn.diagram(list(which(res_DESeq$padj<0.05&abs(res_DESeq$log2FoldChange)>log2(1.25)),which(res_table_evo$padj<0.05&abs(res_table_evo$logFC)>log2(1.25))),
             "./FigureS1_evo.png",width = 8,height = 8,units = "cm",res=600,imagetype = "png",
             category.names = c("DESeq2","edgeR"),main="evolution")
venn.diagram(list(which(res_DESeq2$padj<0.05&abs(res_DESeq2$log2FoldChange)>log2(1.25)),which(res_table_age$padj<0.05&abs(res_table_age$logFC)>log2(1.25))),
             "./FigureS1_age.png",width = 8,height = 8,units = "cm",res=600,imagetype = "png",
             category.names = c("DESeq2","edgeR"),main="age")
venn.diagram(list(which(res_DESeq3$padj<0.05&abs(res_DESeq3$log2FoldChange)>log2(1.25)),which(res_table_evo2$padj<0.05&abs(res_table_evo2$logFC)>log2(1.25))),
             "./FigureS1_mismatched.png",width = 8,height = 8,units = "cm",res=600,imagetype = "png",
             category.names = c("DESeq2","edgeR"),main="mismatched")


#Figure S2
png("./FigureS2.png",width = 16,height = 16,units = "cm",res = 600)
hist(abs(res_table_age$logFC),freq = F,xlim=c(0,3),ylim=c(0,5),breaks = 25,col=alpha("forestgreen",0.75),
     main="Magnitude of expression changes",xlab=expression(Abs(paste(log[2],FC))))
hist(abs(res_table_evo$logFC),col=alpha("salmon",0.75),add=T,freq = F,breaks = 25)
legend("topright",legend = c("age","evo"),fill=alpha(c("forestgreen","salmon"),0.75))
dev.off()

#Figure S3
png("./FigureS3.png",width = 16,height = 16,units = "cm",res = 600)
pie_tab_gg=reshape2::melt(t(pie_tab)/colSums(pie_tab))
gg=ggplot(pie_tab_gg,aes(x="",y=value,fill=Var2))+
  geom_bar(width = 1,stat = "identity")+
  coord_polar("y", start=0)+
  facet_wrap(.~Var1,nrow = 2)+
  scale_fill_manual(values=c(brewer.pal(12,"Paired")[c(3,5,4,8)],"grey"))+
  theme_bw()+
  theme(axis.text.x=element_blank(),axis.ticks = element_blank(),
        panel.grid=element_blank(),axis.title = element_blank(),
        legend.title = element_blank(),strip.text.x = element_text(size = 16),
        legend.text=element_text(size=12))
g=ggplot_gtable(ggplot_build(gg))
stript=which(grepl('strip-t', g$layout$name))
fills=brewer.pal(12,"Paired")[c(6,10,9,5)]
k=1
for (i in stript) {
  j=which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill=fills[k]
  k=k+1
}
grid.draw(g)
dev.off()  

