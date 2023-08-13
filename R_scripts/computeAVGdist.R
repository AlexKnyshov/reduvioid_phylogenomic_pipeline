library(ape)
library(seqinr)
library(phangorn)
args = commandArgs(trailingOnly=TRUE)
aaseqpath <- args[1]
proteinfiles <- list.files(path=aaseqpath, pattern="\\.fas$", full.names=FALSE, recursive=FALSE)

ComputeDistMatrix <- function(alignment) {
	alignmentmx <- as.phyDat.AAbin(alignment)
	dist_mx_dims <- c(length(alignmentmx), length(alignmentmx[[1]]))
	dist_mx <- matrix(NA, nrow =dist_mx_dims[1], ncol=dist_mx_dims[1], dimnames = list(names(alignmentmx),names(alignmentmx)) )
	logicmx <- lower.tri(dist_mx,diag = FALSE) 
	for (row1 in 1:length(dist_mx[,1])){
	  for (col1 in 1:length(dist_mx[1,])){
	    if (logicmx[row1, col1] == TRUE ){
	      dist_mx[row1, col1] <- dist.ml(alignmentmx[c(row1,col1), ], model="LG")
	    }
	  }
	}
	dist_mx[upper.tri(dist_mx, diag = F)] <- t(dist_mx)[upper.tri(dist_mx, diag=F)]
	#diag(dist_mx) <- 0
	return (dist_mx) 
}
outdata <- data.frame(loc=character(),distval = numeric())
for (i in 1:length(proteinfiles)){
	#path
	proteinfile <- paste0(aaseqpath, "/", proteinfiles[i])
	#read
	proteinlocus <- read.alignment(proteinfile, format = "fasta")
	proteinlocus <- as.matrix.alignment(proteinlocus)
	proteinlocus <- toupper(proteinlocus)
	proteinlocus <- replace(proteinlocus, proteinlocus == "X", "-")
	proteinlocus <- as.AAbin.character(proteinlocus)
	proteinlocusdist <- ComputeDistMatrix(proteinlocus)
	meandist <- mean(proteinlocusdist, na.rm=T)
	tempdf <- data.frame(loc=proteinfiles[i],distval = meandist)
	outdata <- rbind(outdata, tempdf)

	write (paste0(proteinfiles[i], "\t", meandist), "")
}
write.csv(outdata, "dist_output.csv")