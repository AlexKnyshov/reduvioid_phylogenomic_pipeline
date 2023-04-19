library(ape)
library(phangorn)
library(seqinr)
args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cat("Syntax: Rscript crosscontTree.R [path to ref tree] [path to gene trees] [path to gene tree names] [path to NT alignment] [path to AA alignment] [locus distance threshold] [global distance threshold]\n")
  cat("Example: Rscript crosscontTree.R concat.tre genetrees.tre genetrees.txt ./NT/ ./AA/ 0.01 0.2\n")
  quit()
}

gtrees <- read.tree(args[2])
names(gtrees) <-  readLines(args[3])
reftree <- read.tree(args[1])
dnaseqpath <- args[4]
aaseqpath <- args[5]

distthresh <- as.numeric(args[6])
globaldistthresh <- as.numeric(args[7])

#check all trees
for (t in 1:length(gtrees)){
  gtname <- names(gtrees)[t]
  print(paste("gene tree name", gtname))
  # get pairwise distance matrix out of the trees (GTR dist in case of RAxML trees)
  tree_dist <- cophenetic(gtrees[[t]])
  diag(tree_dist) <- NA
  tree_dist[upper.tri(tree_dist)] <- NA
  # find distances smaller than threshold
  small_dist <- which(tree_dist < distthresh, arr.ind = TRUE)
  final_outliers <- character()
  #check if any found
  if (length(small_dist) > 0){
    print(paste("outliers detected in", names(gtrees)[t]))
    #get the global tree trimmed to the gene tree tips
    tempreftree <- drop.tip(reftree, reftree$tip.label[!(reftree$tip.label %in% gtrees[[t]]$tip.label)])
    tempreftree_dist <- cophenetic(tempreftree)
    outliers <- character()
    #check all small dist pairs
    for (pair in 1:length(small_dist[,1])){
      #get global tree dist for the pair
      globpairwise <- tempreftree_dist[rownames(tree_dist)[small_dist[pair,1]],colnames(tree_dist)[small_dist[pair,2]]]
      treelen <- gtrees[[t]]$edge.length
      if (globpairwise > globaldistthresh){
        #globally the pair has larger distance than threshold
        #try removing one and checking if RF dist to global tree improves
        out1 <- rownames(tree_dist)[small_dist[pair,1]]
        out2 <- colnames(tree_dist)[small_dist[pair,2]]
        drop1_ref <- drop.tip(tempreftree, out1)
        drop1_loc <- drop.tip(gtrees[[t]], out1)
        rf1 <- RF.dist(unroot(drop1_ref), unroot(drop1_loc))
        drop2_ref <- drop.tip(tempreftree, out2)
        drop2_loc <- drop.tip(gtrees[[t]], out2)
        rf2 <- RF.dist(unroot(drop2_ref), unroot(drop2_loc))
        if (rf1 < rf2) {
          outliers <- c(outliers, out1)
        } else if (rf1 > rf2) {
          outliers <- c(outliers, out2)
        } else {
          outliers <- c(outliers, out1, out2)
        }
      }
    } 
    final_outliers <- unique(outliers)
    print (final_outliers)
  }
  #read alignments
  dnalocus <- read.dna(paste0(dnaseqpath, "/", gtname), format="fasta")
  proteinlocus <- read.alignment(paste0(aaseqpath, "/", gtname), format = "fasta")
  proteinlocus <- as.matrix.alignment(proteinlocus)
  proteinlocus <- toupper(proteinlocus)
  proteinlocus <- as.AAbin.character(proteinlocus)
  #if any outliers detected correct alignments and gene trees
  if (length(final_outliers)>0){
    dnalocus <- dnalocus[!(rownames(dnalocus) %in% final_outliers),]
    proteinlocus <- proteinlocus[!(rownames(proteinlocus) %in% final_outliers),]
    gtrees[[t]] <- drop.tip(gtrees[[t]], final_outliers)
  }
  #output
  write.FASTA(dnalocus, paste0(dnaseqpath, "/", gtname,".edited"))
  write.FASTA(proteinlocus, paste0(aaseqpath, "/", gtname,".edited"))
}
#output corrected gene trees
write.tree(gtrees, paste0(args[2], ".edited"))