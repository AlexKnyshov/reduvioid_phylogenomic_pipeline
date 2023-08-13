# OMCL bait pipeline

Procedures to prepare the OMCL target locus set (2471 loci) from the original orthogroup set (32,675 loci of Gordon, 2017):

## 1. Alignment of all original orthogroups

Starting with OMCL orthogroups of Gordon (2017) (Supplemental file `Gordon2017OMCLorthogroups.tar.gz`), run this sed command on each file to rename the sequence headers to remove `lcl`:
```
sed -e '/>lcl.*/d' INPUT_FILE
```

Then we run MAFFT to align all the sequences with E-INS-i algorithm using the following command:
```
mafft --genafpair --maxiterate 1000 --thread 12 --inputorder INPUT_FILE > OUTPUT_FILE
```
It can be combined with the previous step.

## 2. Run an copy number corrector script

We wanted to select single copy orthologs for downstream procedures, however, at the current step copy number is inflated by 2 issues: 1) multiple isoforms of the same transcript 2) multiple sequences / pieces belonging to the same transcript. We thus used a custom script to collapse isoforms for a given taxon at less or equal 1% raw distance and pick the longest isoform. For remainder of the sequences for the same taxon, we keep and annotate other single copy chunks if they are non overlapping (overlap less or equal 10 AA).
```
mkdir deisoformed_fasta
cd deisoformed_fasta
for f in ./aligned/*.fas; do python deisoform.py ${f}; done
cd ..
mkdir loci_stats
mv deisoformed_fasta/*.txt loci_stats/
```
The script outputs both partly corrected alignments (isoforms removed, partial sequences of same transcript remain, sequence names changed) and a data table on each alignment (TXT file), containing counts of each class of events. We pick only the last column (copy number) for each file and combine that info into a single table, using a custom python script:
```
python copy_counter.py ./loci_stats/
```
The script will produce a corrected table of OG counts for all loci and taxa

## 3. Select OGs, present in 10 or more taxa and single copy

For this, we used the following R script:
```
tab1 <- read.csv("output.csv", stringsAsFactors = F)
c1 <- character()
for(row in 1:length(tab1[,1])){
  #print(tab1[row,1])
  rowtab = as.numeric(tab1[row,2:22])
  tab2 <- as.matrix(table(rowtab))
  cls <- as.numeric(rownames(tab2))
  if (sum(cls > 1)==0){
    if (1 %in% rowtab && tab2["1",] >=10) {
      c1 <- c(c1, tab1[row,1])
    }
  }
}
length (c1)
write(c1, "selected.txt",1)
```
Result contains the names of the OGs passing the filter

## 4. Remove AHE loci

The order of this step could be different, but in our case, this is when we check for AHE loci in our OrthoMCL set and remove those from the `selected.txt` set of selected single copy OGs.

First, run BLASTX and get an output of OG proteins matched at 1e-25 e-value with 95% identity:
```
# run BLASTX of AHE target regions (based on Arilus) against Arilus proteins in OGs
# (totalprots.fas) using ALiBaSeq's blast wrapper:
bash blast_wrapper.sh ./PATH_TO_final_ref_cds/ totalprots.fas 1e-10 blastx 12 y

# run ALiBaSeq to get the list of all found proteins at a given threshold
python alibaseq.py -b totalprots.fas.blast --ac tdna-aa -f S -x a -c=-1 -i 95 -e 1e-25
cut -f1 -d, totalprots.fas.blast_default_ttable.tab > all_found_prots.txt
```

Then find names of orthogroups where these proteins are found (below is the cumbersome collection of scripts showing how it was done in our case):
```
# compile list of seq names in each original OG for seq name to file name correspondence table
grep ">" ORIGINAL_OG_FOLDER/OG1.5_* > allog2protP1.txt

# get paired files of file name / seq name for each sequence in the original OGs
# file names
cut -f1 -d: allog2protP1.txt > allog2protP2.txt
# seq names
cut -f2 -d">" allog2protP1.txt | cut -f1 -d" " > allog2protP3.txt

# get line numbers of AHE-blasted proteins in the seq name file from above,
# then use the line numbers to get matching file names of OGs to be removed, unique it
grep -n -w -f all_found_prots.txt allog2protP3.txt | cut -f1 -d: | while read l; do sed -n ${l}p allog2protP2.txt | rev | cut -f1 -d/ | rev; done | uniq > foundogs.txt
```

Then remove these OGs from selected OGs from above
```
grep -v -f foundogs.txt selected.txt | sort > ahe_screen_selected.txt
mkdir selected_fasta
#copy files from isoform correction output (folder named fasta) that pass the AHE screening filter
while read l; do cp fasta/$l selected_fasta/; done < ahe_screen_selected.txt
```

## 5. Run distance assessment script and remove OGs with large pairwise seq distance

First, run average distance assessment R script on all alignments (ML distance under LG model):
```
Rscript computeAVGdist.R ./selected_fasta/
```
Then only keep OGs with average pairwise ML distance under LG model under 0.95 (determined empirically to remove clear artefacts):
```
#filter the list
awk -F, '{if ($3 <= 0.95) print $2}' dist_output.csv > dist_screen_selected.txt

#copy files that pass the AHE screening filter based on the list
mkdir selected_fasta2
while read l; do cp selected_fasta/$l selected_fasta2/; done < dist_screen_selected.txt
```

## 6. Construct HMM profiles

For each file in selected_fasta2 construct a hmmer profile:
```
mkdir selected_fasta2_profiles/
for f in selected_fasta2/*; do hmmbuild --cpu 3 ./selected_fasta2_profiles/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1,2 -d. )".hmm" $f ; done

```
**Resulting hmm files in `selected_fasta2_profiles` was used for the forward HMMER search in the OrthoMCL pipeline (with the exception of reference assemblies, where the longest AA sequence in each locus was used for the forward TBLASTN search).**

## 7. Obtain sequences of the reference taxon

The most complete sample from the initial OrthoMCL set (Pasir = Pasiropsis sp.) was selected as a reference taxon for the RBH check

### 7.1. Run collapser script to combine transcript pieces

This script will merge fragmentary transcripts for each taxon (if any) into a single transcript. This is most important for the reference taxon, as having most complete target region span helps the reciprocal check later on.
```
mkdir selected_fasta2_collapsed
cd selected_fasta2_collapsed
for f in ./selected_fasta2/*.fas; do python get_consensus.py ${f}; done
```

### 7.2. Extract Pasir taxon sequences

Extract Pasir sequences out of collapsed sequence files from above, one sequence per file

### 7.3. Run TBLASTN search of Pasir regions against Pasir transcriptome

Assuming `EXTRACTED_PASIR_LOCI_FILES` contains Pasir only sequences, we used ALiBaSeq's blast wrapper to run TBLASTN against the corresponding transcriptome at 1e-10 e-value threshold:
```
bash blast_wrapper.sh EXTRACTED_PASIR_LOCI_FILES/ Pasir.Trinity.fasta 1e-10 tblastn 32 y
```
**Resulting file `Pasir.Trinity.fasta.blast` is used for the reciprocal blast results in ALiBaSeq in the OrthoMCL pipeline.**
