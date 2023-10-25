################################################################
#
#                         Configuration
#
#################################################################

# SET THE WORK DIRECTORY
# WORKDIR=""
# setwd(WORKDIR)

################################################################
#
#                         libraries
#
#################################################################

library("ChAMP")
library("RPMM")
library(minfi)
library(ggplot2)
library(DNAmArray)
library(irlba)
library(pheatmap)
library(wateRmelon)
library(ggfortify)
library(BiocManager)
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(shinyMethyl)
library(ggrepel)
#################################################################

#                       DEMO DATASET
# Infinium MethylationEPIC Demo Data Sets can be downloaded
# https://emea.support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html 
# I used Infinium MethylationEPIC v2.0 Demo Data Set (iScan) 821 MB  Nov 17, 2022
# Reference > https://github.com/YuanTian1991/ChAMP-DemoRun/tree/main/EPICv2/illumina_demo_data_iScan 
#################################################################


################################################################
#
#                         ChAMP
#
#################################################################
# .idat or file containing beta values can be loaded with different functions
# ChAMP provides a comprehensive filtering function called champ.filter() which can take any data matrix as input (beta, M, Meth, UnMeth, intensity) and do filtering on it.
# champ.filter() does not call for any specific object or class result, nor does it need an IDAT file, as long as you have a single beta matrix, the filtering can be done. 
# In champ.filter(), since some probes and samples might be removed for reason of low quality, NA might be still exist in filtered data. 
# myImport <- champ.import(testDir)
# myLoad <- champ.filter()


# ADD .idat files to the folder (the RAW data path)
# .idat data can be loaded by pointing the directory to the DataSet like below:

myLoad <- champ.load(directory = getwd(), arraytype = "EPICv2")

# Check the distribution of CpGs on chromosome, CpG island, TSS reagions. e.g. 
# This function can be used for any CpGs list, for example, during your analysis, whenever you get a significant CpG list from DMR, 

CpG.GUI(arraytype="EPICv2")

# Further quality control and exploratory analysis
# Quality Control is an important process to make sure dataset is suitable for downstream analysis. 
# champ.QC() function and QC.GUI() function would draw some plots for user to easily check their data’s quality. 

champ.QC(beta=myLoad$beta, pheno = myLoad$pd$Sample_Group)
QC.GUI(beta=myLoad$beta, pheno = myLoad$pd$Sample_Group, arraytype="EPICv2")

# Normalization: On Illumina bead arrays, probes come in two different designs (called type-I and type-II), with different hybridization chemistries, which means that probes from these two different designs will exhibit different distributions. 
# This is a technical effect and is independent of variations caused by differences in the biological characteristics (e.g. CpG density) of type-I and type-II probes.
# champ.norm() only support BMIQ and PBC method.

myNorm <- champ.norm(myLoad$beta, method = "BMIQ", arraytype = "EPICv2")

myNorm2 <- champ.norm(myLoad$beta, method = "PBC", arraytype = "EPICv2")


# if there is NA in the beta values; fill them with median value > if necessary

# myNorm[is.na(myNorm)] <- median(myNorm, na.rm=TRUE)
# myNorm2[is.na(myNorm2)] <- median(myNorm2, na.rm=TRUE)

write.table(myNorm, "myNorm_16_samples.txt", row.names = TRUE, sep = "\t")  

# 2nd QC after normalization. In the output, it is written before normalization but it is not like that.

champ.QC(beta=myNorm, pheno = myLoad$pd$Sample_Group)
champ.QC(beta=myNorm2, pheno = myLoad$pd$Sample_Group)

# To check the normalized data from QC.GUI()

QC.GUI(beta=myNorm, pheno = myLoad$pd$Sample_Group, arraytype = "EPICv2")
QC.GUI(beta=myNorm2, pheno = myLoad$pd$Sample_Group, arraytype = "EPICv2")


# Assuming you have the 'myNorm' data frame with normalized data

# Find the 1000 most variable CpGs based on standard deviation
most_variable_cpgs <- names(sort(apply(myNorm, 1, sd), decreasing = TRUE))[1:1000]

# Extract the table for the 1000 most variable CpGs
table_1000_cpgs <- myNorm[most_variable_cpgs, ]

# Now you can write this table to a file (e.g., a tab-separated values file)
write.table(table_1000_cpgs, "table_1000_variable_CpGs.txt", row.names = TRUE, sep = "\t")

# If you want to view the table in R, you can print the first few rows using the head() function
head(table_1000_cpgs)


# we noticed that Slide in pd file is numeric, which is not correct, so we firstly manually change it into character.

myLoad$pd$Slide <- as.character(myLoad$pd$Slide)

# champ.SVD() would detect all valid factors to perform analysis, which means the plot contains the two following conditions:
# Covariates are not associated with name, Sample_Name or File_Name
# Covariates contain at least two values(for example, for ‘BeadChip ID’ to be tested as a covariate, samples from at least two different BeadChips must be included in the study)
# If Batch detected, run champ.runCombat() here.

champ.SVD(myNorm, pd=myLoad$pd)

# Cell Type Heterogeneity for blood cells.
# To detect cell proportion and remove cell type influence on WHOLE BLOOD data.

# myRefBase <- champ.refbase(beta=as.data.frame(myNorm) ,arraytype="EPICv2")
# write.csv2(myRefBase$CellFraction,"CellComposition.txt")


# Combat Batch Effect Adjustment; The result shows there are confounding effect for Array, so I want to adjust it with Combat.
# If necessary, run this step
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Array"))

# myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Tissue"))
write.table(myCombat, "myCombat_16_samples.txt", row.names = TRUE, sep = "\t")  


# Differentially Methylated Probes (DMPs) is the most common desired output.
# Here users may use champ.DMP() function to calculate differential methylation, and can use DMP.GUI() to check the result.
# Below code exist because the origin sample group contains space, replace it with `_` works.

tmp_pheno <- myLoad$pd$Sample_Group

tmp_pheno <- gsub(" ", "_", tmp_pheno)

myDMP <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2")

# CpG.GUI() can be used to check the DMP result:
DMP.GUI(DMP=myDMP[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Group)

write.table(myDMP, "001_DMPs_Group_1_vs_Group_2.txt", row.names = TRUE, sep = "\t")


# If you have sample groups, you can compare sample groups also
# myDMP1 <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2", compare.group = c("Group_1", "Group_2"))
# DMP.GUI(DMP=myDMP1[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Group)

# myDMP2 <- champ.DMP(beta = myNorm,pheno=tmp_pheno, arraytype = "EPICv2", compare.group = c("Group_1", "Group_3"))
# DMP.GUI(DMP=myDMP2[[1]], beta=myNorm, pheno=myLoad$pd$Sample_Group)


################################################################
#
#                         Minfi + DNAmArray
#
#################################################################


targets <- read.metharray.sheet(getwd())
RGSet <- read.metharray.exp(targets = targets,force=TRUE,extended=TRUE)

# setting the right array type and annotation will do the work:
RGSet@annotation <- c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38")
RGSet

#Filtering probes
RGSet_filtered=probeFiltering(RGSet)

betas <- getBeta(RGSet, type="Illumina")

# Beta density plots
# Using densityPlot() from minfi (Feinberg, 2014), we can visualize the per sample average beta-value distribution.
# This gives us a global impression of the data and allows us to identify possible anomalous samples. 
# We expect this distribution to be bimodal with the peaks representing methylated and unmethylated signals. 
# Any centre peaks should be further investigated for problems, such as ambiguous mapping.

ggbg <- function() {
  points(0, 0, pch=16, cex=1e6, col="grey90")
  grid(col="white", lty=1)
}

par(mar=c(4,4,3,2), mgp=c(2.5,1,0), 
    cex.main=1.5, font.main="1", 
    fg="#6b6b6b", col.main="#4b4b4b")

densityPlot(RGSet, 
            main="Beta density plot", 
            xlab="Beta values", 
            panel.first=ggbg()) 

# Principal components plot
# Using the prcomp_irlba() function from irlba we can calculate principal components. 
# By assessing the amount of variance explained by these and visualising them, we can better interpret the data. 
# The package ggfortify helps ggplot2 interpret PCA objects, allowing prcomp objects to be passed to the autoplot() function.

pc <- prcomp_irlba(t(betas), n=6)
summary(pc)

# Create the principal components plot
p <- autoplot(pc, data = targets, colour = "Sample_Group", main = "Principal Components Plot")

# Add sample names as labels using ggrepel
p <- p + geom_text_repel(aes(label = Sample_Name), size = 3, show.legend = FALSE, max.overlaps = 20)

# Display the plot
print(p)

# Save the plot as a PDF file
ggsave("PC_plot.pdf", plot = p, width = 8, height = 6)


# shinyMethyl visualization
# shinyMethyl, an interactive visualization R package for exploration and quality control of methylation data. 
# The current version is designed for Illumina Human Methylation 450K arrays.
summarized.data <- shinySummarize(RGSet)
runShinyMethyl(summarized.data)


# QC

MSet <- preprocessRaw(RGSet)
qc <- getQC(MSet)
addQC(MSet,qc)
MSet=fixMethOutliers(MSet, K = -3, verbose = FALSE)

pdf("QC_Report_MSet.pdf")
par(mar=c(3,10,3,10),cex.axis=0.7)
plotQC(qc)
densityPlot(MSet, sampGroups = MSet$Sample_Group)
densityBeanPlot(MSet, sampNames = MSet$Sample_Name)
dev.off()


pdf("QC_Report_RGSet.pdf")
par(mar=c(3,10,3,10),cex.axis=0.7)
plotQC(qc)
densityPlot(MSet, sampGroups = RGSet$Sample_Group)
densityBeanPlot(MSet, sampNames = RGSet$Sample_Name)
dev.off()



#QC report
qcReport(RGSet, pdf= "QC_Report_v2_RGSet.pdf")
par(mar=c(3,10,3,10),cex.axis=0.7)
phenoData <- pData(RGSet)
manifest <- getManifest(RGSet)


#Sex prediction

sex <- getSex(mapToGenome(RGSet))

pdf("Gender.pdf")
ggplot(as.data.frame(sex), aes(x = xMed, y = yMed, fill = predictedSex)) +
  geom_point(shape = 21, size = 3) +
  theme_light() +
  scale_fill_manual(values = c("pink", "blue"), labels = unique(sex$predictedSex)) +
  xlab("X chr, median total intensity (log2)") +
  ylab("Y chr, median total intensity (log2)") +
  guides(fill = guide_legend(title = "Predicted Sex")) +
  theme(legend.position = "top")
dev.off()


# To see sample names
pdf("Gender_sample_names.pdf")
ggplot(as.data.frame(sex), aes(x = xMed, y = yMed, fill = predictedSex)) +
  geom_point(shape = 21, size = 3) +
  geom_text(aes(label = targets$Sample_Name), vjust = 2) +
  theme_light() +
  scale_fill_manual(values = c("pink", "blue"), labels = unique(sex$predictedSex)) +
  xlab("X chr, median total intensity (log2)") +
  ylab("Y chr, median total intensity (log2)") +
  guides(fill = guide_legend(title = "Predicted Sex")) +
  theme(legend.position = "top")
dev.off()


RSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
GRset <- mapToGenome(RSet)


# Normalization Functional Normalization

pc <- screeplot(RGSet)

RGset.funnorm <- preprocessFunnorm(RGSet,sex=sex$predictedSex,nPCs=6)

mybeta <- reduce(RGset.funnorm,RGSet_filtered,what="beta")

write.table(mybeta, "mybeta.txt", row.names = TRUE, sep = "\t")  

most_variable_cpgs_mybeta <- names(sort(apply(mybeta, 1, sd), decreasing = TRUE))[1:1000]

table_1000_cpgs_mybeta <- mybeta[most_variable_cpgs_mybeta, ]

# Now you can write this table to a file (e.g., a tab-separated values file)
write.table(table_1000_cpgs_mybeta, "table_1000_variable_CpGs_mybeta.txt", row.names = TRUE, sep = "\t")


#MDS plot for RAW
pdf("MDS_RAW.pdf", width = 10, height = 8)
mdsPlot(getM(MSet), 
        numPositions = 1000, 
        sampGroups = targets$Sample_Group,
        sampNames = targets$Sample_Name,
        legendPos = "bottomleft",
        main = "RAW")
dev.off()

#MDS plot
pdf("MDS_Funnorm.pdf", width = 10, height = 8)
mdsPlot(getM(RGset.funnorm ), 
        numPositions = 2000, 
        sampGroups = targets$Sample_Group,
        sampNames = targets$Sample_Name,
        legendPos = "bottomleft",
        main = "preprocessFunnorm")
dev.off()


# PCA plot and heatmap

beta <- getBeta(RGset.funnorm)
beta[is.na(beta)] <- median(beta, na.rm=TRUE)
pc <- prcomp_irlba(t(beta), n=6)
summary(pc)
plot(pc$x)

# Making the plot for the correlations of PCs
library(pheatmap)

# Plot the heat map

df <- apply(targets, 2, function(x) as.numeric(factor(x)))

# if there is na in the beta values; fill them with median value
df[is.na(df)] <- median(df, na.rm=TRUE)

keep <- apply(df, 2, sd) > 0    
cxy <- cor(pc$x, scale(df[, keep]))   
pheatmap(cxy[,c(2:3,4,5)], cluster_rows=FALSE)

#Plot dendogram

# Change column names to sequential numbers
colnames(beta) <- 1:ncol(beta)


library(fastcluster)
d <- dist(t(beta),method="euclidean")
fit <- hclust(d, method="average")
Interesting <- sample(colnames(beta),1)

colorLeafs <- function(x) {
  if (is.leaf(x) && attr(x, "label") %in% Interesting) {
    attr(x, "nodePar") <- list(lab.col="red", pch=NA)
  }
  return(x)
} 
dd <- dendrapply(as.dendrogram(fit), colorLeafs)
op <- par(mar=c(10,4,4,2))
plot(dd)



d <- dist(t(beta),method="euclidean")
fit <- hclust(d, method="average")
Interesting <- sample(colnames(beta),1)

colorLeafs <- function(x) {
  if (is.leaf(x)) {
    label <- attr(x, "label")
    if (label %in% colnames(beta)[1:23]) {
      attr(x, "nodePar") <- list(lab.col = "blue", pch = NA)
    } else if (label %in% colnames(beta)[24:41]) {
      attr(x, "nodePar") <- list(lab.col = "red", pch = NA)
    }
  }
  return(x)
}

dd <- dendrapply(as.dendrogram(fit), colorLeafs)
op <- par(mar = c(10, 4, 4, 2))
plot(dd)


#SNPs

snps <- getSnpInfo(GRset)
GRset <- addSnpInfo(GRset)
GRset <- dropLociWithSnps(GRset, snps=c("SBE","CpG"), maf=0)


################################################################
#
#                         Watermelon - Age prediction
#
#################################################################

library(minfi)

familyEPIC <- RGSet
familyEPIC.pf <- pfilter(familyEPIC, logical.return = TRUE)
familyEPIC.dasen.pf <- dasen(familyEPIC.pf)
predictedAges <- agep(familyEPIC.dasen.pf)

pdf("Predicted-ages.pdf")
plot(predictedAges$horvath.age)
dev.off()



################################################################
#
#     Methylations as UCSC tracks to create bigwig files
#
#################################################################
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
data(Locations)

hg38info=seqinfo(BSgenome.Hsapiens.UCSC.hg38)
positions=data.frame(row.names=Locations@rownames,
                     chr=Locations@listData$chr,
                     start=Locations@listData$pos,
                     end=Locations@listData$pos+1,
                     strand=Locations@listData$strand)

for (i in 1:ncol(mybeta)){
  mySample=merge(positions,mybeta[,i],by="row.names")
  names(mySample)[6]="score"
  mySample$Row.names=c()
  name=colnames(mybeta)[i]
  write.table(mySample[,c("chr","start","end","score")],
              sep="\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE,
              file=paste(name,".bedGraph",sep=""))
  myRanges=makeGRangesFromDataFrame(mySample,
                                    keep.extra.columns = TRUE)
  seqinfo(myRanges)=hg38info
  export.bw(myRanges,paste0(name,".bw"))
}

################################################################################


