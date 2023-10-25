
library(minfi)
library(shinyMethyl)
library(ggplot2)
library(DNAmArray)


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


################################################################################
#
#                         Analysing DMPs
#
################################################################################

# Normalization Functional Normalization
RGset.funnorm <- preprocessFunnorm(RGSet,nPCs=6)


#Filtering probes
RGSet_filtered=probeFiltering(RGSet)

mybeta <- reduce(RGset.funnorm,RGSet_filtered,what="beta")

group=data.frame(row.names=row.names(pData(RGset.funnorm)),
                 condition=pData(RGset.funnorm )$Sample_Group,
                 gender=pData(RGset.funnorm)$predictedSex )

dmp <- dmpFinder(mybeta, pheno = group$condition , type = "categorical",qCutoff=1,shrinkVar = TRUE)
write.table(dmp,"DMPs.csv",sep="\t")


################################################################################

################################################################################
#
#                         annotation for the each comparison group
#
################################################################################

library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)

hg38info=seqinfo(BSgenome.Hsapiens.UCSC.hg38)
positions=data.frame(row.names=Locations@rownames,
                     chr=Locations@listData$chr,
                     start=Locations@listData$pos,
                     end=Locations@listData$pos+1,
                     strand=Locations@listData$strand)


library(ChIPseeker)
positions.gr=makeGRangesFromDataFrame(positions)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
positions.annot=annotatePeak(positions.gr,tssRegion=c(-2000, 2000),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
positions.annot.df=as.data.frame(positions.annot,row.names=positions.annot@anno@ranges@NAMES)
positions.annot.df$annotationGroup=ifelse(grepl("Intron",positions.annot.df$annotation),
                                          "Intron",
                                          ifelse(grepl("Exon",positions.annot.df$annotation),"Exon",
                                                 ifelse(grepl("Promoter",positions.annot.df$annotation),"Promoter",
                                                        ifelse(grepl("Downstream",positions.annot.df$annotation),
                                                               "Downstream",positions.annot.df$annotation)
                                                 )
                                          )
)

################################################################################
# violin plots for mybeta
# prepare data frames for each group for comparison
Cell_line <- data.frame(CpG = row.names(mybeta),
                           condition = "Cell Line",
                           methylation = unlist(apply(mybeta, 1, function(x) {mean(x[1:3])})))

Coriell <- data.frame(CpG = row.names(mybeta),
                           condition = "Coriell",
                           methylation = unlist(apply(mybeta, 1, function(x) {mean(x[4:5])})))

Cell_line_Coriell.df <- rbind(Cell_line, Coriell)

# Look up annotations using match instead of sapply
Cell_line_Coriell.df$annotation <- positions.annot.df$annotationGroup[match(Cell_line_Coriell.df$CpG, rownames(positions.annot.df))]

library(ggplot2)
Cell_line_Coriell.df$annotation <- factor(Cell_line_Coriell.df$annotation)

Cell_line_Coriell_violin <- ggplot(Cell_line_Coriell.df, aes(x = annotation, y = methylation, fill = condition, order = annotation)) +
  geom_violin() +
  theme_classic() +
  stat_summary(aes(group = condition), fun.y = mean, geom = "point", shape = 19, size = 2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("white", "grey"))

pdf("violinRegions_Cell_line_Coriell.pdf", width = 15, height = 10)
print(Cell_line_Coriell_violin)
dev.off()

Cell_line_Coriell_violin_v2 <- ggplot(Cell_line_Coriell.df, aes(x = annotation, y = methylation, fill = condition, order = annotation)) +
  geom_violin() +
  theme_classic() +
  stat_summary(aes(group = condition), fun.y = mean, geom = "point", shape = 19, size = 2, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("white", "grey")) +
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5)) +  # Center, bold, and bigger title
  labs(title = "Distribution of DNA methylation for DMPs > Cell_line vs Coriell",  # Adding the title
       x = "Annotation",  # x-axis label
       y = "Methylation")  # y-axis label
pdf("violinRegions_Cell_line_Coriell_v2.pdf", width = 14, height = 7)
print(Cell_line_Coriell_violin_v2)
dev.off()

###################################################################

################################################################################
#
#           Genomic annotations for CHRx and CHRy
#
#################################################################

annEPIC= getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
keepX <- (featureNames(RGset.funnorm) %in% annEPIC$Name[annEPIC$chr %in% 
                                                          c("chrX")])

positionsX= positions[positions$chr=="chrX",]
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

positionsX.gr=makeGRangesFromDataFrame(positionsX)
chrX.anno=annotatePeak(positionsX.gr,tssRegion=c(-2000, 2000),
                       TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("chrX.pdf")
plotAnnoPie(chrX.anno)
dev.off()

positionsY= positions[positions$chr=="chrY",]
positionsY.gr=makeGRangesFromDataFrame(positionsY)
chrY.anno=annotatePeak(positionsY.gr,tssRegion=c(-2000, 2000),
                       TxDb=txdb, annoDb="org.Hs.eg.db")
pdf("chrY.pdf")
plotAnnoPie(chrY.anno)
dev.off()

################################################################################

# To obtain annotations for all autosomes,

library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Assuming you have already loaded the necessary data

# Filter for autosomes (excluding chrX and chrY)
autosomes <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
               "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
               "chr21", "chr22")

keepAutosomes <- (featureNames(RGset.funnorm) %in% annEPIC$Name[annEPIC$chr %in% autosomes])

positionsAutosomes <- positions[positions$chr %in% autosomes,]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

positionsAutosomes.gr <- makeGRangesFromDataFrame(positionsAutosomes)
autosomes.anno <- annotatePeak(positionsAutosomes.gr, tssRegion=c(-2000, 2000),
                               TxDb=txdb, annoDb="org.Hs.eg.db")

pdf("autosomes_annotations.pdf")
plotAnnoPie(autosomes.anno)
dev.off()

################################################################################

# To obtain annotations for 1000 variable CpGs

# Find the 1000 most variable CpGs based on standard deviation
most_variable_cpgs <- names(sort(apply(mybeta, 1, sd), decreasing = TRUE))[1:1000]

# Extract the table for the 1000 most variable CpGs
table_1000_cpgs <- mybeta[most_variable_cpgs, ]

# Now you can write this table to a file (e.g., a tab-separated values file)
write.table(table_1000_cpgs, "table_1000_variable_CpGs.txt", row.names = TRUE, sep = "\t")

# If you want to view the table in R, you can print the first few rows using the head() function
head(table_1000_cpgs)


library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# Assuming you have the 'positions' data frame and 'table_1000_cpgs' data frame

# Extract the CpG names from 'table_1000_cpgs'
cpg_names_1000 <- rownames(table_1000_cpgs)

# Subset 'positions' to include only the rows with CpG names in 'cpg_names_1000'
positions_1000 <- positions[rownames(positions) %in% cpg_names_1000, ]

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

positions_1000.gr <- makeGRangesFromDataFrame(positions_1000)
cpg_names_1000.anno <- annotatePeak(positions_1000.gr, tssRegion=c(-2000, 2000),
                                    TxDb=txdb, annoDb="org.Hs.eg.db")

pdf("1000_variable_CpGs_annotations.pdf")
plotAnnoPie(cpg_names_1000.anno)
dev.off()
################################################################################