library(ape)
library(phangorn)

args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript rfFilter.R [path to concat genetrees] [path to text file with gene names] ([path to ref tree])\n")
  cat("Example: Rscript rfFilter.R genetrees.tre genetrees.txt\n")
  cat("Example: Rscript rfFilter.R genetrees.tre genetrees.txt concattree.tre\n")
  quit()
} else if (length(args) >= 2){
	gtrees <- read.tree(args[1])
	names(gtrees) <-  readLines(args[2])
	if (length(args) == 3){
		#this is mode of comparing all vs ref tree
		tree <- read.tree(args[3])
	}
}

FindOutliers <- function(data) {
  #https://stackoverflow.com/questions/12866189/calculating-the-outliers-in-r
  lowerq = quantile(data, na.rm = T)[2]
  upperq = quantile(data, na.rm = T)[4]
  iqr = upperq - lowerq
  extreme.threshold.upper = (iqr * 3) + upperq
  extreme.threshold.lower = lowerq - (iqr * 3)
  result <- which(data > extreme.threshold.upper)
}

if (length(args) == 3){
	#ref tree mode
	tree_dist <- data.frame(locus=character(), rf = numeric())
	for (t in 1:length(gtrees)){
	  temptree <- drop.tip(tree, tree$tip.label[!(tree$tip.label %in% gtrees[[t]]$tip.label)])
	  distval <- wRF.dist(unroot(gtrees[t]), unroot(temptree))
	  tree_dist <- rbind(tree_dist, data.frame(locus=names(gtrees[t]), rf=distval))
	}
} else {
	#all vs all gene tree mode
	#create distance matrix and iterate over non redundant cells
	#getting comparable trees (same tip sets)
	#and measuring their wRF
	tree_dist_mx <- matrix(NA, nrow =length(gtrees), ncol=length(gtrees), dimnames = list(names(gtrees),names(gtrees)) )
	logicmx <- lower.tri(tree_dist_mx,diag = FALSE) 
	for (row1 in 1:length(tree_dist_mx[,1])){
	  for (col1 in 1:length(tree_dist_mx[1,])){
	    if (logicmx[row1, col1] == TRUE ){
	      temptree1 <- drop.tip(gtrees[[row1]], setdiff(gtrees[[row1]]$tip.label, gtrees[[col1]]$tip.label))
	      temptree2 <- drop.tip(gtrees[[col1]], setdiff(gtrees[[col1]]$tip.label, gtrees[[row1]]$tip.label))
	      if (!is.null(temptree1) & !is.null(temptree2)) {
	        if (length(temptree1$tip.label) > 3 & length(temptree2$tip.label) > 3){
	          tree_dist_mx[row1, col1] <- wRF.dist(unroot(temptree1), unroot(temptree2))  
	        }
	      }
	    }
	  }
	}
	#get mean locus wRF
	tree_dist_mx[upper.tri(tree_dist_mx, diag = F)] <- t(tree_dist_mx)[upper.tri(tree_dist_mx, diag=F)] 
	tree_dist <- data.frame(locus=rownames(tree_dist_mx), rf = apply(tree_dist_mx, 1, mean, na.rm = T))
	#in the all vs all mode also save distance matrix
	write.csv(tree_dist_mx, "rf_dist_mx.csv")
}
#sort loci by wRF and output outliers
tree_dist <- tree_dist[order(-tree_dist$rf),]
write(as.character(tree_dist$locus[FindOutliers(tree_dist$rf)]), "rf_loci.txt")