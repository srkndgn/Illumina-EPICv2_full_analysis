
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library("ChAMP")
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(minfi)
library(ggplot2)
library(devtools)
library(DNAmArray)
library(ChIPseeker)
library(dplyr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(Gviz)
library(GenomicRanges)
library(ComplexHeatmap)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
library(BSgenome.Hsapiens.UCSC.hg38)


###############################################################################
data(Locations)
positions=data.frame(row.names=Locations@rownames,
                     chr=Locations@listData$chr,
                     start=Locations@listData$pos,
                     end=Locations@listData$pos+1,
                     strand=Locations@listData$strand)


positions.gr=makeGRangesFromDataFrame(positions)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
positions.annot=annotatePeak(positions.gr,tssRegion=c(-2500, 2500),
                             TxDb=txdb, annoDb="org.Hs.eg.db")
positions.annot.df=as.data.frame(positions.annot,row.names=positions.annot@anno@ranges@NAMES)

#select pcdha and promoter
PD <- positions.annot.df[grep("PCDHA|PCDHB|PCDHG", positions.annot.df$SYMBOL), ]
PD2 <- PD[,c('annotation', 'SYMBOL', 'start', 'end')]
PD2 <- PD2[grep("Promoter", PD$annotation), ]
PD2_unique <- distinct(PD2, SYMBOL, .keep_all= TRUE)

#combine individual beta values (beta) with pd2, from UCSC script
combo <- merge(PD2, mybeta, by = 0)

###################################################################################################################
#make average beta value for each pcdha gene promoter (sample number + 5 columns)
sorted_combo <- combo[order(combo$SYMBOL),]
j <- 1
average_table <- data.frame("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21")

for (i in 1:nrow(sorted_combo)){
  s <- sorted_combo[i,3]
  print(s)
  x <- grep(s,sorted_combo$SYMBOL) 
  print(x)
  average <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
  for (k in x){
    for (m in 6:21){
      average[m]=average[m]+sorted_combo[k,m]
    }
  }
  average = average/length(x)
  average[3] = s
  average_table[j,]=average
  j = j + 1
}

average_table <- unique(average_table)
average_table <- average_table[,-c(1,2,4,5)]

# Rename 'X.3.' column to 'SYMBOL'
names(average_table)[names(average_table) == 'X.3.'] <- 'SYMBOL'

# Merge the data frames
data <- merge(PD2_unique, average_table, by = "SYMBOL")

# Rename other columns
names(data)[names(data) == 'X.6.'] <- 'Cell_Line_1'
names(data)[names(data) == 'X.7.'] <- 'Cell_Line_2'
names(data)[names(data) == 'X.8.'] <- 'Cell_Line_3'
names(data)[names(data) == 'X.9.'] <- 'Coriell_1'
names(data)[names(data) == 'X.10.'] <- 'Coriell_2'
names(data)[names(data) == 'X.11.'] <- 'Cell_Line_4'
names(data)[names(data) == 'X.12.'] <- 'Cell_Line_5'
names(data)[names(data) == 'X.13.'] <- 'Cell_Line_6'
names(data)[names(data) == 'X.14.'] <- 'Cell_Line_7'
names(data)[names(data) == 'X.15.'] <- 'Cell_Line_8'
names(data)[names(data) == 'X.16.'] <- 'Cell_Line_9'
names(data)[names(data) == 'X.17.'] <- 'Coriell_3'
names(data)[names(data) == 'X.18.'] <- 'Coriell_4'
names(data)[names(data) == 'X.19.'] <- 'Cell_Line_10'
names(data)[names(data) == 'X.20.'] <- 'Cell_Line_11'
names(data)[names(data) == 'X.21.'] <- 'Cell_Line_12'

write.table(sorted_combo, "sorted_combo.txt",
            row.names = TRUE, sep = "\t")

write.table(data, "PCDH_a_b_g_data.txt",
            row.names = TRUE, sep = "\t")
##################################

##########
#hg38

genome <- "hg38"
chrom <- "chr"
from <- 140770957
to <- 141529546

axTrack <- GenomeAxisTrack()
idxTrack <- IdeogramTrack(genome = "hg38", chromosome = "chr5")

#plotTracks(list(idxTrack, axTrack, knownGenes), from = from, to = to, showTitle = FALSE)

refGenes <- UcscTrack(genome = "hg38", chromosome = "chr5",
                      track = "xenoRefGene", from = from, to = to,
                      trackType = "GeneRegionTrack", rstarts = "exonStarts",
                      rends = "exonEnds", gene = "name", symbol = "name2",
                      transcript = "name", strand = "strand", fill = "#8282d2",
                      stacking = "dense", name = "PCDHA PCDHB PCDHG", background.panel = "white",
                      background.title = "black")

plotTracks(refGenes)


cpgIslands <- UcscTrack(genome = "hg38", chromosome = "chr5", 
                        track = "cpgIslandExt", from = from, to = to,
                        trackType = "AnnotationTrack", 
                        start = "chromStart", end = "chromEnd", 
                        id = "name", shape = "box", fill = "black",
                        name = "CpG Islands", background.panel = "white",
                        background.title = "black",)
plotTracks(cpgIslands)
## make data track
dTrack <- DataTrack(data,chromosome = "chr5", genome = "hg38", name = "<<  M E T H Y L A T I O N  L E V E L  (beta value 0-1)  >>", background.panel = "white",
                    background.title = "black", groups = rep(c("01-Cell Line", "02-Cell Line", "03-Coriell", "04-Coriell", "05-Cell Line", "06-Cell Line", "07-Cell Line", "08-Cell Line", "09-Cell Line", "10-Cell Line", "11-Cell Line", "12-Cell Line", "13-Cell Line", "14-Cell Line", "15-Cell Line", "16-Cell Line", "17-Coriell", "18-Coriell", "19-Cell Line", "20-Cell Line", "21-Cell Line"
)), type = "boxplot", type = c("a", "p"), ylim = c(0, 1), col=c("orange", "orange", "orange", "purple", "purple", "orange", "orange", "orange","orange", "orange", "orange","orange", "orange", "orange","orange", "orange", "purple","purple", "orange", "orange","orange",  ))

plotTracks(list(idxTrack, axTrack, refGenes, cpgIslands, dTrack), from = from, to = to, showTitle = TRUE)
plotTracks(list(idxTrack, axTrack, cpgIslands, dTrack), from = from, to = to, showTitle = TRUE)

#########################
#########################

###### making heatmaps with beta values##################
PD <- positions.annot.df[grep("PCDHA|PCDHB|PCDHG", positions.annot.df$SYMBOL), ]
combo <- merge(PD, mybeta, by = 0)



#select the beta values,
combo2 <- data[,c("Cell_Line_1", "Cell_Line_2", "Cell_Line_3", "Coriell_1", "Coriell_2", "Cell_Line_4", "Cell_Line_5", "Cell_Line_6", "Cell_Line_7", "Cell_Line_8", "Cell_Line_9", "Coriell_3", "Coriell_4", "Cell_Line_10", "Cell_Line_11", "Cell_Line_12"
)]
combo2_matrix <- data.matrix(data)
Heatmap(combo2, name = "beta", column_title = "Heatmap according to 561 CpGs located in the promoter regions(2500-TSS+2500) of PCDHAs, PCDHBs and PCDHGs ",show_row_names = FALSE)
