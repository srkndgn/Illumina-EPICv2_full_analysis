################################################################################
library(minfi)
library(shinyMethyl)
library(ggplot2)
library(DNAmArray)
################################################################################
# Correlation
# If you have many samples, it takes too much time to sa as PDF !!!!

getCorrelation=function(meth.mat,method="spearman",plot=TRUE,nrow=2e6){
  
  print( cor(meth.mat,method=method) )
  panel.cor.pearson <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,method="pearson"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  panel.cor.kendall <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,method="kendall"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  panel.cor.spearman <- function(x, y, digits=2, prefix="", cex.cor, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,method="spearman"))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * r)
  }
  
  panel.my.smooth2<-function(x, y, col = par("col"), 
                             bg = NA, pch = par("pch"), cex = 1, 
                             col.smooth = "darkgreen", span = 2/3, 
                             iter = 3, ...) 
  {
    par(new = TRUE)    #par(usr = c(usr[1:2], 0, 1.5) )
    smoothScatter(x, y,colramp=colorRampPalette(topo.colors(100)), bg = bg)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)) 
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...)
    abline(lm(y[ok]~x[ok]), col="red")
  }
  
  panel.my.smooth<-function(x, y, col = par("col"), bg = NA, 
                            pch = par("pch"), cex = 0.3, 
                            col.smooth = "green", span = 2/3, iter = 3, ...) 
  {
    points(x, y, pch = 20, col = densCols(x,y,
                                          colramp=colorRampPalette(topo.colors(20))), 
           bg = bg, cex = 0.1)
    ok <- is.finite(x) & is.finite(y)
    if (any(ok)){
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),col = col.smooth, ...);
      abline(lm(y[ok]~x[ok]), col="red")}
  }
  panel.hist <- function(x, ...)
  {
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
  }
  
  if(plot)
  {  
    if(method=="spearman")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.spearman,
            diag.panel=panel.hist,main=paste(
              method,"cor.") )
    }
    if(method=="kendall")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.kendall,
            diag.panel=panel.hist,main=paste(
              method,"cor.") )
    }
    if(method=="pearson")
    { pairs(meth.mat, 
            lower.panel=panel.my.smooth2, 
            upper.panel=panel.cor.pearson,
            diag.panel=panel.hist,main=paste(
              method,"cor.") )
    }
  }
}  
data(Locations)
positions=data.frame(row.names=Locations@rownames,
                     chr=Locations@listData$chr,
                     start=Locations@listData$pos,
                     end=Locations@listData$pos+1,
                     strand=Locations@listData$strand)

# Enable parallelization
require(doParallel)
registerDoParallel(cores = 4)

sex=getSex(mapToGenome(RGSet))
RGSet.funnorm <- preprocessFunnorm(RGSet,sex=sex$predictedSex,nPCs=4)

annEPIC= getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
keep <- !(featureNames(RGSet.funnorm) %in% annEPIC$Name[annEPIC$chr %in% 
                                                          c("chrX","chrY")])

################################################################################
pdf("EPPICv2_correlation_Autosomes.pdf")
getCorrelation(meth.mat=getBeta(RGSet.funnorm[keep,]))
dev.off()
################################################################################

cairo_pdf("EPPICv2_correlation_Autosomes_v2.pdf", width = 10, height = 8)  
getCorrelation(meth.mat = getBeta(RGSet.funnorm[keep, ]))
dev.off()
################################################################################
cairo_pdf("EPPICv2_correlation_Autosomes_v3.pdf")  # Output as PDF file
getCorrelation(meth.mat = getBeta(RGSet.funnorm[keep, ]))
dev.off()
################################################################################

# Number of samples per page in the PDF

samples_per_page <- 5

pdf("EPPICv2_correlation_Autosomes_splitted.pdf")  # Output as PDF file
for (i in seq(1, ncol(RGSet.funnorm), samples_per_page)) {
  end_idx <- min(i + samples_per_page - 1, ncol(RGSet.funnorm))
  getCorrelation(meth.mat = getBeta(RGSet.funnorm[keep, i:end_idx]))
}
dev.off()

################################################################################