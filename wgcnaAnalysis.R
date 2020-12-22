graphics.off()
rm(list = ls())

library(data.table)
library("anRichment")
library(org.Hs.eg.db)
library('openxlsx')
##########################
#### Read in the data:
##########################

#First, read in the data
miRNA <- fread("miRBase.21.RPKM.csv.gz")
RNA   <- fread("ENSG.GRCh38.82.RPKM.csv.gz")

#Remove the columns containing the labels, to leave just the numeric entries
miRNA_matrix <- as.matrix(miRNA[,-c(1,2)])
RNA_matrix   <- as.matrix(RNA[,-c(1,2)])

#Put back the labels as rownames
rownames(miRNA_matrix) <- miRNA$mirbase_id
rownames(RNA_matrix)   <- RNA$ensembl_gene_id

#Find the IDs of the individuals for which we have both the mRNA and miRNA measurements
commonIDs <- intersect(colnames(RNA_matrix), colnames(miRNA_matrix))

#Extract the data for the common set of individuals:
miRNA_commonIDs <- miRNA_matrix[, commonIDs]
RNA_commonIDs   <- RNA_matrix[, commonIDs]

#Create a single data matrix containing the mRNA and miRNA:
allRNA <- rbind(miRNA_commonIDs, RNA_commonIDs)

# Read in some additional info re: case/ctl status
metaData <- read.csv("metaData.csv", stringsAsFactors = F)
rownames(metaData) <- metaData$Sample.ID
metaData           <- metaData[colnames(allRNA),]

#Filename for WGCNA:
myFilename <- "SuppInfoFile3"


##########################
#### Screen/filter the transcripts:
##########################

#Screen according to the same criteria used elsewhere in the paper:
# For both approaches, we used transcripts passing following filters: 
# (1) FPM ≥ 0.2 for non-circRNA; count ≥5 for circRNA, 
# (2) present in > 10% of the cohort, 

#We suggest removing features whose counts are consistently low (for example, removing all features that have a count of less than say 10 in more than 90% of the samples) because such low-expressed features tend to reflect noise and correlations based on counts that are mostly zero aren't really meaningful. 

transcriptsToKeep     <- which(rowSums(allRNA > 0.2) > 0.1*ncol(allRNA))

allRNA_filtered       <- allRNA[transcriptsToKeep,]
# Keep track of the miRNAs that we have kept:
retained_miRNA        <- intersect(rownames(miRNA_matrix), rownames(allRNA_filtered))
# Log transform:
allRNA_logTransformed <- log2(allRNA_filtered + 1)

##########################
#### Choose power for WGCNA
##########################

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

datExpr <- t(allRNA_logTransformed) # We require columns = genes, rows = samples

# We need to choose a soft-thrsholding power:
# 
# Define a set of soft-thresholding powers to choose from
powers = c(c(1:10), seq(from = 12, to=20, by=1))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
pdf(file = paste0(myFilename,"-choosingSoftThresholdPower.pdf"))
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


##########################
#### Construct WGCNA network:
##########################

net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste0(myFilename,"-TOM"), 
                       verbose = 3)

save.image(paste0(myFilename, "-net100.RData"))
load(paste0(myFilename, "-net100.RData"))


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
pdf(paste0(myFilename, "-clusterDendrogram.pdf"))
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#=====================================================================================
#
#  Save down some specific info
#
#=====================================================================================
load(paste0(myFilename, "-net100.RData"))

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = paste0(myFilename, "-networkConstruction-auto.RData") )



## Create a blank workbook
wb <- createWorkbook()

for(i in sort(unique(net$colors)))
{
    myDF <- names(net$colors)[net$colors == i]
    
    addWorksheet(wb, paste0("membersOfModule_", i))        
    writeData(wb, paste0("membersOfModule_", i), myDF)
    
}

## Save workbook to working directory
saveWorkbook(wb, file = paste0(myFilename, ".xlsx"), overwrite = TRUE) 


#=====================================================================================
#
#  Assess which modules are enriched for miRNAs
#
#=====================================================================================

pvalVector_mirna_enrichment <- vector(mode = "numeric", length = length(unique(net$colors)))
myMatrix <- matrix(nrow = 3, ncol = length(unique(net$colors)))
for(i in 1:length(unique(net$colors)))
{
    myMatrix[1,i] <- sum(net$colors == (i-1))
    myMatrix[2,i] <- sum((names(net$colors) %in% retained_miRNA) & (net$colors == (i-1)) )
    pvalVector_mirna_enrichment[i] <- (fisher.test(x =as.factor(names(net$colors) %in% retained_miRNA), y = as.factor(net$colors == (i-1)), alternative = "greater")$p.value)
}
print(which(pvalVector_mirna_enrichment < 0.05/length(unique(net$colors))) - 1)
myMatrix[3,] <- round(p.adjust(pvalVector_mirna_enrichment, method = "bonferroni"),3)

rownames(myMatrix) <- c("Module size", "Num. miRNA in module", "Bonf. adj. p-value")
colnames(myMatrix) <- paste0("Module", seq(0, max(net$colors)))

wb <- loadWorkbook(file = paste0(myFilename, ".xlsx"))

addWorksheet(wb, "miRNA_enrichedModules")
writeData(wb, "miRNA_enrichedModules", myMatrix, rowNames = T)
## Save workbook to working directory
saveWorkbook(wb, file = paste0(myFilename, ".xlsx"), overwrite = TRUE) 


#=====================================================================================
#
#  Assess which modules are associated with outcome
#
#=====================================================================================
library(pheatmap)
library(qvalue)

caseCtrlStatus <- metaData[rownames(net$MEs),4]
caseCtrlBinary <- caseCtrlStatus != "CTL"
myDf                         <- cbind(net$MEs, as.factor(as.numeric(caseCtrlBinary)))
myDfnames                    <- names(myDf)
myDfnames[length(myDfnames)] <- "y"
names(myDf)                  <- myDfnames

pValVector        <- vector(mode = "numeric", length = length(net$MEs))
names(pValVector) <- names(net$MEs) 
counter         <- 0 
for(i in names(net$MEs) )
{
        counter <- counter + 1
        expr    <- as.formula(paste0("y ~ ", i))
        myglm   <- glm(expr, family = "binomial", data  = myDf)
        pValVector[counter] <- ((coef(summary(myglm))[2,4]))
}

qValVector <- qvalue(pValVector)$qvalues


caseCtrlStatus_noPE <- caseCtrlStatus[caseCtrlStatus != "PE"]
myDf_noPE           <- myDf[caseCtrlStatus != "PE" , ]

pValVector_noPE        <- vector(mode = "numeric", length = length(net$MEs))
names(pValVector_noPE) <- names(net$MEs) 
counter         <- 0 
for(i in names(net$MEs) )
{
        counter <- counter + 1
        expr    <- as.formula(paste0("y ~ ", i))
        myglm   <- glm(expr, family = "binomial", data  = myDf_noPE)
        pValVector_noPE[counter] <- ((coef(summary(myglm))[2,4]))
}

qValVector_noPE <- qvalue(pValVector_noPE)$qvalues




caseCtrlStatus_noFGR <- caseCtrlStatus[caseCtrlStatus != "FGR"]
myDf_noFGR           <- myDf[caseCtrlStatus != "FGR" , ]

pValVector_noFGR        <- vector(mode = "numeric", length = length(net$MEs))
names(pValVector_noFGR) <- names(net$MEs) 
counter         <- 0 
for(i in names(net$MEs) )
{
        counter <- counter + 1
        expr    <- as.formula(paste0("y ~ ", i))
        myglm   <- glm(expr, family = "binomial", data  = myDf_noFGR)
        pValVector_noFGR[counter] <- ((coef(summary(myglm))[2,4]))
}

qValVector_noFGR <- qvalue(pValVector_noFGR)$qvalues




peFgrStatus     <- caseCtrlStatus[caseCtrlStatus == "FGR" | caseCtrlStatus == "PE"]

myDf                         <- cbind(net$MEs[caseCtrlStatus == "FGR" | caseCtrlStatus == "PE",], as.factor(as.numeric(peFgrStatus == "PE")))
myDfnames                    <- names(myDf)
myDfnames[length(myDfnames)] <- "y"
names(myDf)                  <- myDfnames

pValVector_peFgr        <- vector(mode = "numeric", length = length(net$MEs))
names(pValVector_peFgr) <- names(net$MEs) 
counter         <- 0 
for(i in names(net$MEs) )
{
        counter <- counter + 1
        expr    <- as.formula(paste0("y ~ ", i))
        myglm   <- glm(expr, family = "binomial", data  = myDf)
        pValVector_peFgr[counter] <- ((coef(summary(myglm))[2,4]))
}

qValVector_peFgr <- qvalue(pValVector_peFgr)$qvalues



qValMx <- rbind(qValVector, qValVector_noFGR, qValVector_noPE, qValVector_peFgr)
rownames(qValMx) <- c("Case_Ctrl", "PE_Ctrl", "FGR_Ctrl", "PE_FGR")
wb2    <- loadWorkbook(file = paste0(myFilename, ".xlsx"))
addWorksheet(wb2, "eigengeneAssociationTests_qVals")
writeData(wb2, "eigengeneAssociationTests_qVals", round(qValMx,4), rowNames = T)
style010 <- createStyle(bgFill = "#FFFDD0")
style005 <- createStyle(bgFill = "#FCF4A3")
style001 <- createStyle(bgFill = "#FFF200")

conditionalFormatting(wb2, "eigengeneAssociationTests_qVals", cols =1:(ncol(qValMx)+1), rows = 1:(nrow(qValMx)+1), rule = "<0.1", style = style010)
conditionalFormatting(wb2, "eigengeneAssociationTests_qVals", cols =1:(ncol(qValMx)+1), rows = 1:(nrow(qValMx)+1), rule = "<0.05", style = style005)
conditionalFormatting(wb2, "eigengeneAssociationTests_qVals", cols =1:(ncol(qValMx)+1), rows = 1:(nrow(qValMx)+1), rule = "<0.01", style = style001)
saveWorkbook(wb2, paste0(myFilename, ".xlsx"), TRUE)  #


# Visulaise the distributions of the eigengenes for modules 10 and 11

myDfall <- cbind(net$MEs, as.factor(caseCtrlStatus))
myDfnames                    <- names(myDfall)
myDfnames[length(myDfnames)] <- "y"
names(myDfall)               <- myDfnames
myDfall2 <- myDfall[, c(paste0("ME", seq(1,22)), "y")]
suppressWarnings(melted <- melt(myDfall2))
names(melted) <- c("Group", "variable", "Value")



pdf(file = "boxplot_me10squashed.pdf", width = 7, height = 3.5)
p <- ggplot(melted[melted$variable=="ME10",], aes(x=variable, y=Value, fill=Group)) + 
    geom_boxplot() + theme_bw()+ theme( text = element_text(size = 20), axis.text = element_text(size = 20), axis.title.x = element_blank()) +
    labs(title="Module 10 : Eigengene analysis")+ coord_cartesian(ylim = c(-0.15, 0.35))
plot(p)
dev.off()

pdf(file = "boxplot_me11squashed.pdf", width = 7, height = 3.5)
p <- ggplot(melted[melted$variable=="ME11",], aes(x=variable, y=Value, fill=Group)) + 
    geom_boxplot() + theme_bw()+ theme( text = element_text(size = 20), axis.text = element_text(size = 20), axis.title.x = element_blank()) +
    labs(title="Module 11 : Eigengene analysis")+ coord_cartesian(ylim = c(-0.15, 0.35))
plot(p)
dev.off()

#=====================================================================================
#
#  Perform enrichment  analysis using gprofiler2
#
#=====================================================================================

## Create a blank workbook
wb <- createWorkbook()
for(i in 1:length(unique(net$colors)))
{
    print(i)
    evcodesFlag <- FALSE
    currentGeneList   <- names(net$colors[net$colors == (i-1)])
    gostres           <- gost(currentGeneList, evcodes = evcodesFlag)
    newdf             <- gostres$result
    newdf$adj_p_value <- formatC(newdf$p_value, format = "e", digits = 1)
    
    addWorksheet(wb, paste0("enrichmentAnalysisModule_", i-1))  
    if(evcodesFlag)
    {
        writeData(wb, paste0("enrichmentAnalysisModule_", i-1), newdf[,c("source", "term_id", "term_name", "term_size", "intersection_size", "query_size","adj_p_value", "intersection")])
    }
    else
    {
        writeData(wb, paste0("enrichmentAnalysisModule_", i-1), newdf[,c("source", "term_id", "term_name", "term_size", "intersection_size", "query_size","adj_p_value")])
    }
    
}

## Save workbook to working directory
saveWorkbook(wb, file = paste0(myFilename, "_enrichmentAnalysis.xlsx"), overwrite = TRUE) #Supplementary Table yy3

wb_go <- loadWorkbook(file = paste0(myFilename, "_enrichmentAnalysis.xlsx"))

wb_comb <-loadWorkbook(file =  paste0(myFilename, ".xlsx"))


for(i in 1:length(sheets(wb_go)))
{
    currentData <- read.xlsx(paste0(myFilename, "_enrichmentAnalysis.xlsx"), sheet = sheets(wb_go)[i])
    addWorksheet(wb_comb, sheets(wb_go)[i])
    writeData(wb_comb, sheets(wb_go)[i], currentData)
}

saveWorkbook(wb_comb, file = paste0(myFilename, ".xlsx"), overwrite = TRUE) 

#=====================================================================================
#
#  Assess whether modules 10 and 11 are enriched for the targets of the miRNAs they contain
#
#=====================================================================================

library(biomaRt)
library(gprofiler2)
library(clusterProfiler)
library(ggplot2)

myMatrix <- matrix(nrow = 3, ncol = 2)
rownames(myMatrix) <- c("Module size", "For miRNA in module, num. targets also in module", "Bonf. adj. p-value")
colnames(myMatrix) <- paste0("Module", seq(10, 11))

load("filtered.miRTarBase_SE_WR_hsa.RData")
filtered_miRNA_Targets    <- dt.mtb.f

miRNAs_with_known_targets <- intersect(names(net$colors[net$colors == 10]), filtered_miRNA_Targets$miRNA)
targets <- filtered_miRNA_Targets$`Target Gene`[which(filtered_miRNA_Targets$miRNA %in% miRNAs_with_known_targets)]

convertedTargets       <- gconvert(unique(targets))
target_mRNAInDataset   <- intersect(convertedTargets$target, names(net$colors)) 
target_mRNAInModule10  <- intersect(convertedTargets$target, names(net$colors[net$colors == 10])) 
pval_module10          <- fisher.test(y =as.factor(names(net$colors) %in% target_mRNAInDataset), x = as.factor(net$colors == 10), alternative = "greater")$p.value
print(pval_module10)

myMatrix[1,1] <- sum(net$colors == 10)
myMatrix[2,1] <- length(target_mRNAInModule10)


miRNAs_with_known_targets <- intersect(names(net$colors[net$colors == 11]), filtered_miRNA_Targets$miRNA)
targets <- filtered_miRNA_Targets$`Target Gene`[which(filtered_miRNA_Targets$miRNA %in% miRNAs_with_known_targets)]

convertedTargets       <- gconvert(unique(targets))
target_mRNAInDataset   <- intersect(convertedTargets$target, names(net$colors)) 
target_mRNAInModule11  <- intersect(convertedTargets$target, names(net$colors[net$colors == 11])) 
pval_module11          <- fisher.test(y =as.factor(names(net$colors) %in% target_mRNAInDataset), x = as.factor(net$colors == 11), alternative = "greater")$p.value
print(pval_module11)

myMatrix[1,2] <- sum(net$colors == 11)
myMatrix[2,2] <- length(target_mRNAInModule11)

myMatrix[3,]  <- round(p.adjust(c(pval_module10, pval_module11), method = "bonferroni"),4)


wb <- loadWorkbook(file = paste0(myFilename, ".xlsx"))
addWorksheet(wb, "mirnaTargets_enrichedModules")
writeData(wb, "mirnaTargets_enrichedModules", myMatrix, rowNames = T)
## Save workbook to working directory
saveWorkbook(wb, file = paste0(myFilename, ".xlsx"), overwrite = TRUE) 

#=====================================================================================
#
#  For the miRNAs contained in module 10, there are 11 targets also in module 10.
#  Perform enrichment analysis for these 11 mRNAs.
#
#=====================================================================================

gostres <- gost(target_mRNAInModule10)
newdf   <- gostres$result
newdf$p_value <- formatC(newdf$p_value, format = "e", digits = 1)
newdf[,c("source", "term_id", "term_name", "term_size", "intersection_size", "query_size","p_value")]
newdf$adj_p_value <- formatC(newdf$p_value, format = "e", digits = 1)

wb <- loadWorkbook(file = paste0(myFilename, ".xlsx"))
addWorksheet(wb, "module10mirnaTargets_enrichment")
writeData(wb,"module10mirnaTargets_enrichment" , newdf[,c("source", "term_id", "term_name", "term_size", "intersection_size", "query_size","adj_p_value")])
## Save workbook to working directory
saveWorkbook(wb, file = paste0(myFilename, ".xlsx"), overwrite = TRUE) 

#=====================================================================================
#
#  Visualise the enrichment analysis results
#
#=====================================================================================



miRNAs_with_known_targets <- intersect(names(net$colors[net$colors == 11]), filtered_miRNA_Targets$miRNA)
targets <- filtered_miRNA_Targets$`Target Gene`[which(filtered_miRNA_Targets$miRNA %in% miRNAs_with_known_targets)]


level4and5termsBP     <- unique(c(clusterProfiler:::getGOLevel("BP", 4), clusterProfiler:::getGOLevel("BP", 4)))
level12and3termsBP    <- unique(c(clusterProfiler:::getGOLevel("BP", 1), clusterProfiler:::getGOLevel("BP", 2), clusterProfiler:::getGOLevel("BP", 3)))
minlevel4and5termsBP  <- setdiff(level4and5termsBP, level12and3termsBP)

level4and5termsMF     <- unique(c(clusterProfiler:::getGOLevel("MF", 4), clusterProfiler:::getGOLevel("MF", 4)))
level12and3termsMF    <- unique(c(clusterProfiler:::getGOLevel("MF", 1), clusterProfiler:::getGOLevel("MF", 2), clusterProfiler:::getGOLevel("MF", 3)))
minlevel4and5termsMF  <- setdiff(level4and5termsMF, level12and3termsMF)

level4and5termsCC     <- unique(c(clusterProfiler:::getGOLevel("CC", 4), clusterProfiler:::getGOLevel("CC", 4)))
level12and3termsCC    <- unique(c(clusterProfiler:::getGOLevel("CC", 1), clusterProfiler:::getGOLevel("CC", 2), clusterProfiler:::getGOLevel("CC", 3)))
minlevel4and5termsCC  <- setdiff(level4and5termsCC, level12and3termsCC)

minlevel4and5termsALL <- unique(c(minlevel4and5termsBP, minlevel4and5termsCC, minlevel4and5termsMF))



gostres <- gost(names(net$colors[net$colors == 10])[ !names(net$colors[net$colors == 10]) %in% retained_miRNA])

tmpDf               <- gostres$result
tmpDfreduced        <- tmpDf[tmpDf$term_id %in% minlevel4and5termsALL | (tmpDf$source == "REAC" | tmpDf$source == "WP"),]
tmpDfreducedBP      <- tmpDfreduced[tmpDfreduced$source == "GO:BP",]
tmpDfreduced        <- tmpDfreduced[-which(tmpDfreduced$p_value > tmpDfreducedBP$p_value[10] & tmpDfreduced$source == "GO:BP"),]
tmpDfreduced$source[tmpDfreduced$source == "GO:CC"] <- "HO:CC"
tmpDfreduced        <- tmpDfreduced[sort.int(tmpDfreduced$source, index.return = T)$ix, ]
tmpDfreduced$source[tmpDfreduced$source == "HO:CC"] <- "GO:CC"
tmpDfreduced$source <- factor(tmpDfreduced$source, c("GO:BP", "GO:MF", "GO:CC", "REAC",  "WP"))
tmpDfreduced$V4     <- factor(tmpDfreduced$term_name, levels=tmpDfreduced$term_name) # convert V2 into factor
names(tmpDfreduced)[which(names(tmpDfreduced) == "p_value")] <- "adj_P"
library(stringr)



p <- ggplot(tmpDfreduced, aes(x=V4, y=-log10(adj_P), fill=source)) +
    geom_bar(stat="identity")  + scale_fill_manual(values=c(rgb(101,163,198, max = 255), rgb(102,185,118, max = 255), rgb(254,163,101, max = 255), rgb(192,192,192, max = 255), rgb(231,107,243, max = 255),  rgb(248,118,109, max = 255), rgb(231,107,243, max = 255), rgb(192,192,192, max = 255), rgb(231,107,243, max = 255)))+
    #geom_bar(stat="identity")  + scale_fill_manual(values=c(rgb(101,163,198, max = 255), rgb(248,118,109, max = 255), rgb(231,107,243, max = 255)))+
    geom_hline(aes(yintercept = -log10(0.05)), color="black", linetype="dashed") + theme_bw() +
    scale_x_discrete(limits = (levels(tmpDfreduced$V4)), labels = function(x) str_wrap(x, width = 35)) + coord_cartesian(ylim = c(0, 20))+
    #theme(axis.text = element_text(size = 8), axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    theme(axis.text = element_text(size = 18), axis.title.x = element_blank(), axis.text.x = element_text(angle = 60,  hjust=1, size = 8), legend.text = element_text(size = 14), axis.title.y = element_text(size = 20))+ labs(fill = "") + ylab(bquote(~log[10]*'(adj. p-value)'))#expression("-log[10](adj. p-value)"))
#    geom_hline(data = tmpDf, aes(yintercept = 0.5, color = "red")
plot(p)
pdf(file = "enrichment10D_45_35.pdf", width = 13, height =4)
plot(p)
dev.off()

gostres <- gost(names(net$colors[net$colors == 11])[ !names(net$colors[net$colors == 11]) %in% retained_miRNA])


tmpDf               <- gostres$result
tmpDfreduced        <- tmpDf[tmpDf$term_id %in% minlevel4and5termsALL | (tmpDf$source == "REAC" | tmpDf$source == "WP"),]
tmpDfreduced$source <- as.factor(tmpDfreduced$source)
tmpDfreduced$V4     <- factor(tmpDfreduced$term_name, levels=tmpDfreduced$term_name) # convert V2 into factor
names(tmpDfreduced)[which(names(tmpDfreduced) == "p_value")] <- "adj_P"


p <- ggplot(tmpDfreduced, aes(x=V4, y=-log10(adj_P), fill=source)) +
    geom_bar(stat="identity")  + scale_fill_manual(values=c(rgb(101,163,198, max = 255), rgb(102,185,118, max = 255), rgb(254,163,101, max = 255), rgb(192,192,192, max = 255), rgb(231,107,243, max = 255),  rgb(248,118,109, max = 255), rgb(231,107,243, max = 255), rgb(192,192,192, max = 255), rgb(231,107,243, max = 255)))+
    #geom_bar(stat="identity")  + scale_fill_manual(values=c(rgb(101,163,198, max = 255), rgb(248,118,109, max = 255), rgb(231,107,243, max = 255)))+
    geom_hline(aes(yintercept = -log10(0.05)), color="black", linetype="dashed") + theme_bw() +
    scale_x_discrete(limits = (levels(tmpDfreduced$V4)), labels = function(x) str_wrap(x, width = 31)) + coord_cartesian(ylim = c(0, 10))+
    theme(axis.text = element_text(size = 18), axis.title.x = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, size = 8), legend.text = element_text(size = 14), axis.title.y = element_text(size = 20))+ labs(fill = "") + ylab(bquote(~log[10]*'(adj. p-value)'))+ scale_y_continuous(breaks=c(0,2.5,5,7.5,10), labels = c("0","","5","","10"))#expression("-log[10](adj. p-value)"))
#plot(p)

pdf(file = "enrichment11D45_31.pdf", width = 13, height =4)
plot(p)
dev.off()


gostres <- gost(target_mRNAInModule10)

tmpDf               <- gostres$result
tmpDfreduced        <- tmpDf[tmpDf$term_id %in% minlevel4and5termsALL | (tmpDf$source == "REAC" | tmpDf$source == "WP"),]
tmpDfreduced$source <- as.factor(tmpDfreduced$source)
tmpDfreduced$V4     <- factor(tmpDfreduced$term_name, levels=tmpDfreduced$term_name) # convert V2 into factor
names(tmpDfreduced)[which(names(tmpDfreduced) == "p_value")] <- "adj_P"


p <- ggplot(tmpDfreduced, aes(x=V4, y=-log10(adj_P), fill=source)) +
    geom_bar(stat="identity")  + scale_fill_manual(values=c(rgb(101,163,198, max = 255), rgb(254,163,101, max = 255), rgb(102,185,118, max = 255), rgb(248,118,109, max = 255), rgb(231,107,243, max = 255),  rgb(248,118,109, max = 255), rgb(231,107,243, max = 255), rgb(248,118,109, max = 255), rgb(231,107,243, max = 255)))+
    #geom_bar(stat="identity")  + scale_fill_manual(values=c(rgb(101,163,198, max = 255), rgb(248,118,109, max = 255), rgb(231,107,243, max = 255)))+
    geom_hline(aes(yintercept = -log10(0.05)), color="black", linetype="dashed") + theme_bw() +
    coord_flip() + scale_x_discrete(limits = rev(levels(tmpDfreduced$V4))) +
    theme(axis.text = element_text(size = 5), axis.title.y = element_blank() ) 

#    geom_hline(data = tmpDf, aes(yintercept = 0.5, color = "red")
pdf(file = "enrichment10targetsB.pdf")
plot(p)
dev.off()

