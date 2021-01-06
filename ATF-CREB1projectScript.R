library("affy")
library("annaffy")
library("limma")
library("GenomicRanges")
library("reshape2")
library("hgu133a2.db")
library("hgu133a.db")
library("u133aaofav2cdf")
library("hugene10stv1cdf")
library("hugene10stprobeset.db")
library("hugene10stv1probe")
library("fields")
library("gplots")
library("heatmap.plus")
library("annotationTools")
library(calibrate)
library(lattice)
library(directlabels)
library(plotrix)
library(amap)
library(ggplot2)
################################################################################
target <- "EWS-ATF1vsOther_1FC01Padj"
setwd ("~/Documents/R/Report/ATF1-CREB1project")

################################################################################
u133aID2ucscID <- read.table("/Volumes/HHD/Users/yun-shaosung/HDD_Doc/R/Annotation/u133a2ucsc/u133aID2ucscID.txt", header = TRUE, sep = "\t")
u133aSymbol2ucscID <- read.table("/Volumes/HHD/Users/yun-shaosung/HDD_Doc/R/Annotation/u133a2ucsc/u133aSymbol2ucscID.txt", header = TRUE, sep = "\t")
rownames(u133aID2ucscID) <- u133aID2ucscID$U133A_ID


anno<-read.csv('/Volumes/HHD/Users/yun-shaosung/HDD_Doc/R/Annotation/HG-U133A.na33.annot.csv/HG-U133A.na33.annot.csv',colClasses='character',comment.char='#')
rownames(anno)<- anno$Probe.Set.ID

################################################################################
setwd ("/Volumes/My Passport/PeterAnalysis/GeneSet/")

#Is this code useful? I didn't see allarray in the following script.
allarray <- vector(mode="list", length=20) 

filenames <- list.files(path=".", pattern = "*CEL")
exprArray <- ReadAffy(filenames = as.matrix(filenames))
setwd ("~/Documents/R/Report/ATF1-CREB1project")
normExpr <- expresso (exprArray, bgcorrect.method = "rma", normalize.method = "quantiles", pmcorrect.method = "pmonly", summary.method = "medianpolish")
normExprMatr <-exprs(normExpr)

pdf("boxplot.pdf", width = 15, height = 10)
par(mar=c(20, 7, 20, 2.1))
boxplot(exprs(normExpr), las=3)
boxplot(exprArray, las=3)
dev.off()

################################################################################
exprs.data <- exprs(normExpr)
probe2symb <- aafSymbol(rownames(exprs.data), "hgu133a.db")
symbols <- sapply(probe2symb, function(x) {ifelse(length(x) != 0, x, NA)})  #it is list data type, use symbols[[1]] to call the elements
names(symbols) <- rownames(exprs.data)
aggregate.data <- aggregate(exprs.data, by = list(symbols), FUN = median) #looks like symbol become a variable, which is not good
rownames(aggregate.data) <- aggregate.data[[1]]  #so it add a the symbol to become rownames
exprs.by.symbol <- aggregate.data[, -1]

rownames(Symbol2ID) <- as.matrix(sapply(probe2symb, function(x) {ifelse(length(x) != 0, x, NA)}))  #it is list data type, use symbols[[1]] to call the elements
Symbol2ID <- as.matrix(rownames(exprs.data))

ID2Symbol <- as.matrix(sapply(probe2symb, function(x) {ifelse(length(x) != 0, x, NA)}))  #it is list data type, use symbols[[1]] to call the elements
rownames(ID2Symbol) <- as.matrix(rownames(exprs.data))

# Some useful conversion tools
symb2probe <- as.matrix(rownames(exprs.data))
rownames(symb2probe) <- as.matrix(sapply(probe2symb, function(x) {ifelse(length(x) != 0, x, NA)}))
probe2symb2 <- as.matrix(sapply(probe2symb, function(x) {ifelse(length(x) != 0, x, NA)})) 
rownames(probe2symb2) <- as.matrix(rownames(exprs.data))

################################################################################
mytest=function(x){
  test1=x[grepl("EWS-ATF1", colnames(exprs.by.symbol))] #EWS-CREB1-AFH
  test2=x[!grepl("EWS-ATF1|EWS-CREB1", colnames(exprs.by.symbol))]
  p=t.test(test1, test2)$p.value
  foldchange=mean(test1)-mean(test2)
  c(foldchange, p)
}
results=t(apply(exprs.by.symbol,1,mytest))
adjustedp=p.adjust(results[,2], method = "fdr")  #this consider the FDR to do p-value correction
results <- cbind(results, adjustedp)
rownames(results) <- rownames(exprs.by.symbol)
colnames(results) <- c("logFC", "p", "adjustedp")
exprs.Array.FC <- cbind(exprs.by.symbol, results)

##########################################################################################
#get gene list

#There are two different ways to get the gene list, which one should be used?
#first gene list type?
genelist=exprs.Array.FC[(exprs.Array.FC[,"logFC"]< -1 | exprs.Array.FC[,"logFC"]> 1) & (exprs.Array.FC[,"adjustedp"]<0.01),]
genelist <- genelist[order(genelist[,"logFC"], decreasing = TRUE),]
nrow(genelist)

genes <- c("CREB1", "EWSR1", "ATF1")
genes <- c("NTF3", "POU2AF1", "IRF4", "PAX5", "CCND1", "ALK", "MN1", "NTRK3", "ERBB3", "TACC2")

#second gene list type?
genelist <- as.matrix(exprs.by.symbol[genes,])
genelist <- 2^genelist
##########################################################################################
write.table(genelist, file = sprintf("Table_%s_%sgenes.csv", target, nrow(genelist)) , sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)
genelist <- subset(genelist, select = -c(logFC, p, adjustedp))

#Is the motifResult the same as the one appearing later?
genelist<- 2^(exprs.by.symbol[rownames(motifResult), ])

##########################################################################################
labels_Category <- colnames(genelist)
Sample_color <-c(1:ncol(genelist))
Sample_color[grepl("*", labels_Category)]="gray"
Sample_color[grepl("EWS-ATF1", labels_Category)]="#4aaeee"
Sample_color[grepl("EWS-CREB1-AFH", labels_Category)]="#00ff00"
Sample_color[grepl("EWS-CREB1-GI", labels_Category)]="#ff6666"

filename<- sprintf("Bar_CSTB_CTSB_CTBS_%s_%sgenes.PDF", target, nrow(genelist)) 
pdf(filename, width = 20, height = 10)
for (i in 1:nrow(genelist)){
  par(mfrow=c(1,2))
  #par(mar=c(10, 10, 10, 5))  
  barplot(genelist[i,], main=rownames(genelist)[i],  width=0.1,  ylab="Intensity", col=Sample_color, cex.axis = 1.2, cex.names =1.2, las=3, cex.main=2)
  boxplot(c((genelist[i,][grepl("EWS-ATF1", colnames(genelist))]), 
            (genelist[i,][grepl("EWS-CREB1-AFH", colnames(genelist))]), 
            (genelist[i,][grepl("EWS-CREB1-GI", colnames(genelist))]), 
            (genelist[i,][!grepl("EWS-ATF1|EWS-CREB1-AFH|EWS-CREB1-GI", colnames(genelist))]) )~c(
              rep("EWS-ATF1", length(colnames(genelist)[grepl("EWS-ATF1", colnames(genelist))])), 
              rep("EWS-CREB1-AFH", length(colnames(genelist)[grepl("EWS-CREB1-AFH", colnames(genelist))])), 
              rep("EWS-CREB1-GI", length(colnames(genelist)[grepl("EWS-CREB1-GI", colnames(genelist))])), 
              rep("Rest", length(colnames(genelist)[!grepl("EWS-ATF1|EWS-CREB1-AFH|EWS-CREB1-GI", colnames(genelist))]) )  ) , main=rownames(genelist)[i], ylab="Intensity", col=c("gray"), cex =3, , cex.main=2, cex.axis = 1.2)
  stripchart(c((genelist[i,][grepl("EWS-ATF1", colnames(genelist))]), 
               (genelist[i,][grepl("EWS-CREB1-AFH", colnames(genelist))]), 
               (genelist[i,][grepl("EWS-CREB1-GI", colnames(genelist))]), 
               (genelist[i,][!grepl("EWS-ATF1|EWS-CREB1-AFH|EWS-CREB1-GI", colnames(genelist))]) )~c(
                 rep("EWS-ATF1", length(colnames(genelist)[grepl("EWS-ATF1", colnames(genelist))])), 
                 rep("EWS-CREB1-AFH", length(colnames(genelist)[grepl("EWS-CREB1-AFH", colnames(genelist))])), 
                 rep("EWS-CREB1-GI", length(colnames(genelist)[grepl("EWS-CREB1-GI", colnames(genelist))])), 
                 rep("Rest", length(colnames(genelist)[!grepl("EWS-ATF1|EWS-CREB1-AFH|EWS-CREB1-GI", colnames(genelist))]) )  ),vertical = TRUE, method = "jitter", bg = "bisque",pch = 19, col = c("#4aaeee", "#00ff00", "#ff6666", "gray"), add = TRUE, cex =1.2) 
}
dev.off()

##########################################################################################
labels_Category <- colnames(genelist)
Sample_color <-c(1:ncol(genelist))
Sample_color[grepl("*", labels_Category)]="gray"
Sample_color[grepl("EWS-ATF1", labels_Category)]="#4aaeee"
Sample_color[grepl("EWS-CREB1-AFH", labels_Category)]="#00ff00"
Sample_color[grepl("EWS-CREB1-GI", labels_Category)]="#ff6666"

filename<- sprintf("Heat_%s_%sgenes.PDF", target, nrow(genelist)) 
pdf(filename, width = 10, height = 10)
par(mar=c(5, 5, 5, 5))
#heatmap.2(genelist, col = bluered(128), main=filename, trace="none", density.info="density", ColSideColors=Sample_color, symm=F,symkey=F,symbreaks=T, scale="none")
heatmap.2(genelist, col = bluered(128), main=filename, trace="none", density.info="density", scale="row", ColSideColors=Sample_color)
legend("topright",legend=c("EWS-ATF1", "EWS-CREB1-AFH", "EWS-CREB1-GI", "Other"), 
       fill=c("#4aaeee", "#00ff00", "#ff6666","gray"), border=FALSE, bty="n", y.intersp = 0.7, cex=1)
dev.off()

##########################################################################################

#what is the function of this part?

ews_atf <- read.csv("./Table_EWS-ATF1vsOther_1FC01Padj_386genes.csv")
test <- as.matrix(u133aID2ucscID[symb2probe[rownames(ews_atf),],2])
test <- cbind(u133aID2ucscID[symb2probe[rownames(ews_atf),],2], u133aID2ucscID[symb2probe[rownames(ews_atf),],3], ews_atf)
write.table(genelist, file = sprintf("Table_EWS-ATF1vsOther_1FC01Padj_386genes_motif.tsv", target, nrow(genelist)) , sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)

ews_creb <- read.csv("./Table_EWS-CREB1-AFHvsOther_1FC01Padj_141genes.csv")
mergeGenes <- as.matrix(c(rownames(ews_creb), genes))
rownames(mergeGenes) <- mergeGenes[,1]
ews_creb <- exprs.Array.FC[rownames(mergeGenes),]
test <- as.matrix(u133aID2ucscID[symb2probe[rownames(ews_creb),],2])
test <- cbind(u133aID2ucscID[symb2probe[rownames(ews_creb),],2], u133aID2ucscID[symb2probe[rownames(ews_creb),],3], ews_creb)
write.table(test, file = sprintf("Table_EWS-CREB1-AFHvsOther_1FC01Padj_141genes_plusRnaSeqGenes_motif.tsv", target, nrow(test)) , sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

#What is the file "MotifMatchingResult_EWS-CREB1-AFHvsOther_1FC01Padj_141genes_plusRnaSeqGenes_Cytoscape.txt"? How was it generated?
motifResult <- read.csv("MotifMatchingResult_EWS-CREB1-AFHvsOther_1FC01Padj_141genes_plusRnaSeqGenes_Cytoscape.txt", row.names = NULL)
motifResult[motifResult[,6]=="",6] <- NA 
motifResult <- as.matrix(na.omit(motifResult[,6]))
rownames(motifResult)<- motifResult[,1]

for(i in 1:length(motifResult))
{
  print(rownames(motifResult)[i])
  print (exprs.Array.FC[rownames(motifResult)[i],])
}

##########################################################################################
ews_atf <- read.csv("./Table_EWS-ATF1vsOther_1FC01Padj_386genes.csv")
ews_atf_anno <- cbind(anno[symb2probe[rownames(ews_atf),],"Alignments"], ews_atf)
write.table(ews_atf_anno, file = "Table_EWS-ATF1vsOther_1FC01Padj_386genes_anno.csv" , sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)

ews_creb1 <- read.csv("./Table_EWS-CREB1-AFHvsOther_1FC01Padj_141genes.csv")
ews_creb1_anno <- cbind(anno[symb2probe[rownames(ews_creb1),],"Alignments"], ews_creb1)
write.table(ews_creb1_anno, file = "Table_EWS-CREB1-AFHvsOther_1FC01Padj_141genes.csv" , sep = ",", quote = FALSE, col.names = TRUE, row.names = TRUE)

##########################################################################################
setwd('~/Documents/R/Report/ATF1-CREB1project/')
groupName1 <- "CCS"
groupName2 <- "AFH"
otherName <- "REST"

myTarget1 <- "EWS-ATF1"
myTarget2 <- "AFH"
refNotTarget <- "EWS-ATF1|AFH"



#AFH specific
genes<-c("SGK1", "MXRA5", "CTSB", "SOX10" )
genelist <- as.matrix(exprs.by.symbol[genes,])
genelist <- 2^genelist
mytable <- genelist
mytable <- mytable[, !grepl("GI-CCS", colnames(mytable))]
mytable <- mytable[, !grepl("CCS", colnames(mytable))]
mytable <- mytable[, !grepl("AS_CA_U133A_AS20.CEL|AS_CA_U133A_AS5.CEL|SFT_CA_SFT7_HG-U133A.CEL|MFH_NS_U133A_973.CEL|GISTKIT_CA_U133A_28.CEL|Normal_NS_U133A_kidney.CEL|MFH_NS_U133A_14263rep.CEL|AS_CA_U133A_AS27.CEL|AS_CA_U133A_AS34.CEL|AS_CA_U133A_AS27.CEL|Normal_NS_U133A_adrenal.CEL|SFT_CA_SFT7_HG-U133A.CEL", colnames(mytable))]
titleName <- "AFH_Geneset"

#CCS specific
genes<-c("MITF", "PMEL", "SOX10", "SLC7A5", "DUSP4")
genelist <- as.matrix(exprs.by.symbol[genes,])
genelist <- 2^genelist
mytable <- genelist
mytable <- mytable[, !grepl("GI-CCS", colnames(mytable))]
mytable <- mytable[, !grepl("AFH", colnames(mytable))]
mytable <- mytable[, !grepl("GISTKIT_CA_U133A_28.CEL|GISTKIT_CA_U133A_43.CEL|SBRCT_SY_SBRCT11_HG-U133A.CEL|SBRCT_SY_SBRCT5_HG-U133A.CEL", colnames(mytable))]
titleName <- "CCS_Geneset_test"


colnames(mytable)[!grepl(refNotTarget, colnames(mytable))] <- otherName
colnames(mytable)[grepl(myTarget1, colnames(mytable))] <- groupName1
colnames(mytable)[grepl(myTarget2, colnames(mytable))] <- groupName2

full.m <- melt(as.matrix(mytable),id.vars="Name", measure.vars=c(groupName1, groupName2, otherName))
full.m$X1 <- factor(full.m$X1, levels=genes, labels=genes)
colnames(full.m) <- c("gene", "group", "val")

e <- ggplot(full.m, aes(x = group, y = val), labs(title = titleName)) +
  geom_boxplot(
    aes(color = gene), width = 0.5
    , position = position_dodge(0.9)
  ) 
ggsave(e, device = "pdf", filename = sprintf("%s.pdf", titleName) , dpi = 500, width = 5, height = 8)


##########################################################################################
###### Revise ######

#Is the difference between revise and former part just reserving "AFH" from "mytable"?
setwd('~/Documents/R/Report/ATF1-CREB1project/')
groupName1 <- "CCS"
groupName2 <- "AFH"
otherName <- "REST"

myTarget1 <- "EWS-ATF1"
myTarget2 <- "AFH"
refNotTarget <- "EWS-ATF1|AFH"


#CCS specific
genes<-c("MITF", "PMEL", "SOX10", "SLC7A5", "DUSP4")
genelist <- as.matrix(exprs.by.symbol[genes,])
genelist <- 2^genelist
mytable <- genelist
mytable <- mytable[, !grepl("GI-CCS", colnames(mytable))]
#mytable <- mytable[, !grepl("AFH", colnames(mytable))]
mytable <- mytable[, !grepl("GISTKIT_CA_U133A_28.CEL|GISTKIT_CA_U133A_43.CEL|SBRCT_SY_SBRCT11_HG-U133A.CEL|SBRCT_SY_SBRCT5_HG-U133A.CEL", colnames(mytable))]
titleName <- "CCS_Geneset_test"


colnames(mytable)[!grepl(refNotTarget, colnames(mytable))] <- otherName
colnames(mytable)[grepl(myTarget1, colnames(mytable))] <- groupName1
colnames(mytable)[grepl(myTarget2, colnames(mytable))] <- groupName2

full.m <- melt(as.matrix(mytable),id.vars="Name", measure.vars=c(groupName1, groupName2, otherName))
full.m$X1 <- factor(full.m$X1, levels=genes, labels=genes)
colnames(full.m) <- c("gene", "group", "val")

e <- ggplot(full.m, aes(x = group, y = val), labs(title = titleName)) +
  geom_boxplot(
    aes(color = gene), width = 0.5
    , position = position_dodge(0.9)
  ) 
ggsave(e, device = "pdf", filename = sprintf("%s.pdf", titleName) , dpi = 500, width = 5, height = 8)
