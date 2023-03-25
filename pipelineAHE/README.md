# AHE pipeline

Pipeline for processing AHE data


## paths setup

```
repo_folder=""
alibaseq_folder=""
blast=""
raxml=""
iqtree=""
macseFolder=""
Nthreads=""
baitfolderCDS=""
baitfolderAA=""
transcriptomicAssembliesFolder=""
transcriptomicSampleList=""
genomicAssembliesFolder=""
genomicSampleList=""
referenceAssembliesFolder=""
referenceSampleList=""
reciprocalReference=""
reciprocalRefBlast=""
```


## 1. Run BLAST (forward search)
Run appropriate BLAST search (in the case of conserved AHE loci dc-megablast did best, minimizing spurious false positive matches due to conserved AA sequence) on all assemblies and group results in folders. We group by data type, since downstream steps slightly differ for different sequencing types. We set e-value to 1e-05.
```
#transcriptomes
bash ${alibaseq_folder}/blast_wrapper.sh \
	${baitfolderCDS} \
	${transcriptomicAssembliesFolder}/ \
	1e-05 \
	dc-megablast \
	${Nthreads} \
	y \
	${transcriptomicSampleList}
mkdir blastT/
mv *.blast blastT/

#genomes and AHE data
bash ${alibaseq_folder}/blast_wrapper.sh \
	${baitfolderCDS} \
	${genomicAssembliesFolder}/ \
	1e-05 \
	dc-megablast \
	${Nthreads} \
	y \
	${genomicSampleList}
mkdir blastG/
mv *.blast blastG

#reference genomes
bash ${alibaseq_folder}/blast_wrapper.sh \
	${baitfolderCDS} \
	${referenceAssembliesFolder}/ \
	1e-05 \
	dc-megablast \
	${Nthreads} \
	y \
	${referenceSampleList}
mkdir blastB/
mv *.blast blastB
```


## 2. Decontamination
### 2.1. Pull out genomic contigs to do decontamination
Pull out all matched contigs from genomic / AHE assemblies to do contamination check for non-arthropod lineages since our data had some of these issues. ALiBaSeq pulls out original contigs grouping them by sample.
```
python ${alibaseq_folder}/alibaseq.py \
	-b ./blastG/ \
	-f M \
	-x n \
	-t ${genomicAssembliesFolder}/ \
	-e 1e-05 \
	-o ./alibaseq_outGens/ \
	-c 0 \
	--om target \
	-s gens
```


### 2.2. Run decontamination megablast against NCBI
This approach was inspired by https://github.com/damurdock/SIDR

While the original approach of inferring unknown contaminants did not work for our data, the part to identify known contaminants worked well. For that, megablast is run against NCBI and taxon ids of matches are recorded
```
mkdir megaGout

blastn \
	-task megablast \
	-query ./alibaseq_outGens/$sample \
	-db $NCBI_DB/nt \
	-outfmt '6 qseqid staxids' \
	-culling_limit 5 \
	-evalue 1e-25 \
	-out "./megaGout/"$sample".blast" \
	-num_threads ${Nthreads}
```


### 2.3. Run decontamination analysis python helper script
This script analyses the blast results from the previous step together with taxon dump file from NCBI (rankedlineage.dmp) to find contaminant contigs
```
python ${repo_folder}/python_scripts/parsedecontblast.py rankedlineage.dmp megaGout/
```


### 2.4. Filter out contaminant matches from blast results
In order to save time, instead of removing contaminants from assemblies and redoing the search, original blast files from step 1 are filtered
```
mkdir decontblastG

for f in blastG/*.blast
do
	fname=$(basename ${f})
	grep -v -f "megaGout/"$fname"_decont_out.txt" $f > "decontblastG/"$fname
done
```


## 3. Run reciprocal search
Output (decontaminated) blast tables are used to check all contigs with matches against the reference transcriptome for a reciprocal best hit (RBH) check which greatly improves orthology prediction accuracy.
```
#transcriptomes
bash ${alibaseq_folder}/reciprocal_search.sh \
	./blastT/ \
	${transcriptomicAssembliesFolder}/ \
	${reciprocalReference} \
	dc-megablast \
	${Nthreads} \
	y \
	${alibaseq_folder}/reciprocal_get_contigs.py \
	${transcriptomicSampleList}

#genomes
bash ${alibaseq_folder}/reciprocal_search.sh \
	./decontblastG/ \
	${genomicAssembliesFolder}/ \
	${reciprocalReference} \
	dc-megablast \
	${Nthreads} \
	y \
	${alibaseq_folder}/reciprocal_get_contigs.py \
	${genomicSampleList}
```
For large chromosome-level assemblies, the contigs are too large, so to speed up the search and lower memory requirements, chromosomes are chunked into smaller pieces that are then blasted. Each sample was submitted as its own job on our compute cluster.
```
#ref genomes
bash ${repo_folder}/TBA/reciprocal_chunky_search.sh \
	${referenceAssembliesFolder}/$sample \
	${reciprocalReference} \
	dc-megablast \
	${Nthreads} \
	${repo_folder}/TBA/chunkify.py \
	y
```


## 4. Run ALiBaSeq to extract loci
This compiles per locus sets of sequences for all samples. We've noticed that for our purposes BLAST search was not accurate to extract all CDS sequence while minimizing intronic / unalignable sequence. Thus, ALiBaSeq was run in b mode to extract all sequence between outer most matched coordinates for a given contig. CDS was then spliced out at the next step.
```
# transcriptomes
python ${alibaseq_folder}/alibaseq.py \
	-b ../blastT/ \
	-f M \
	-x b \
	-t ${transcriptomicAssembliesFolder}/ \
	-o trans_out/ \
	-e 1e-15 \
	-c 1 \
	-r ../blastT/ \
	-R ${reciprocalRefBlast} \
	--is \
	--ac dna-dna \
	-s trans \
	--rm-rec-not-found \
	--amalgamate-hits
mv *.log logs/
mv *.tab logs/

# reference genomes
python ${alibaseq_folder}/alibaseq.py \
	-b ../blastB/ \
	-f M \
	-x b \
	-t ${referenceAssembliesFolder}/ \
	-q ./trans_out/ \
	-o big_out/ \
	-e 1e-10 \
	-c 1 \
	-r ../blastB/ \
	-R ${reciprocalRefBlast} \
	--amalgamate-hits \
	--ac dna-dna \
	-s bigs \
	--rm-rec-not-found \
	--hit-ovlp 50 \
	--max-gap 20000
mv *.log logs/
mv *.tab logs/

# genomes and AHE
python ${alibaseq_folder}/alibaseq.py \
	-b ../decontblastG/ \
	-f M \
	-x b \
	-t ${genomicAssembliesFolder}/ \
	-q ./big_out/ \
	-o gens_out_b_is/ \
	-e 1e-10 \
	-c 1 \
	--is \
	-r ../decontblastG/ \
	-R ${reciprocalRefBlast} \
	--amalgamate-hits \
	--ac dna-dna \
	-s gens \
	--rm-rec-not-found \
	--hit-ovlp 50 \
	--ctg-ovlp 50
mv *.log logs/
mv *.tab logs/
```

## 5. Run exonerate pipeline
Since the previous search was intentionally run to extract as much sequence as possible, but also many sequences were too long to be aligned first and then trimmed (MSA-based method of homologous sequence recovery), we used exonerate to splice CDS of all sequences first based on the reference taxon, before performing a multiple sequence alignment (MSA). In our dataset with many divergent taxa we had to use exact alignment strategy in exonerate 2.2.0 to get maximal sequence retrieval. This step additionally filtered out spurious matches. Obtained CDS were translated in the first frame, as exonerate returned in frame results.
```
mkdir exoner_temp
mkdir exoner_concat
mkdir exoner_out
cd exoner_temp
echo "debug set up" > ../exoner_debug.txt

for x in $(ls ../gens_out_b_is/*.fas)
do
	#get the reference
	echo $x >> ../exoner_debug.txt
	refname=${baitfolderAA}"/"$(basename ${x})
	#splitting
	python ${repo_folder}/python_scripts/split_fasta.py $x
	for f in *
	do
		exonerate \
			--model protein2genome \
			-q $refname \
			-t $f \
			-E \
			-n 1 \
			--showalignment false \
			--showvulgar false \
			--verbose 0 \
			--ryo ">%ti\n%tcs" 1> ../exoner_concat/$f 2>> ../exoner_debug.txt
	done
	cat ../exoner_concat/* > ../exoner_out/$(basename ${x})
	rm *
	rm ../exoner_concat/*
done
cd ..
rmdir exoner_temp
rmdir exoner_concat

## check empty files and delete:

cd exoner_out
find . -type f -exec awk -v x=1 'NR==x{exit 1}' {} \; -exec rm -f {} \;
cd ..

## translate
python ${repo_folder}/python_scripts/translator.py ./exoner_out/ -u
mv translated translated1
```


## 6. Alignment and block trimming
### 6.1. MSA on AA
Align translated alignments. The value of --unalignlevel was determined empirically to reduce spurious alignment and assist trimming later on.
```
## align proteins
mkdir realignedAA

for f in translated1/*
do
	outf=$(basename ${f})
	mafft \
		--globalpair \
		--thread ${Nthreads} \
		--unalignlevel 0.8 \
		$f > realignedAA/$outf
done
```


### 6.2. Align NT based on AA
CDS of loci were aligned based on AA using MACSE
```
## align NT
mkdir realignedNT

for f in realignedAA/*
do
	sample=$(echo $f | rev | cut -d/ -f1 | rev)
	java \
		-Xmx7g \
		-jar ${macseFolder}/macse_v2.03.jar \
		-prog reportGapsAA2NT \
		-align_AA $f \
		-seq nt_translator/$sample \
		-out_NT realignedNT/$sample
done

```


### 6.3. Alignment block trimming
A custom script was used to do per position trimming of AA alignments at 40% missing data or less, removing sequences with less than 10% data. Then samples that were completely removed from AA were also removed from NT using a custom script, and MACSE was used to mask the remaining CDS. Since MACSE does additional base masking (isolated bases), AA alignment is no longer in sync with CDS after this. So after MACSE trimmed CDS, these files were translated to create new AA files.
```
## trim AA alignments
python ${repo_folder}/python_scripts/customtrim.py \
	./realignedAA/ \
	-d prot
mv trimmed trimmedAA


#remove completely trimmed out taxa
python ${repo_folder}/python_scripts/removeTaxa.py \
	./realignedNT/ \
	-m ./trimmedAA/
mv rmtaxaout realignedNTeq


#run MACSE to trim NT based on AA
mkdir trimmedNT
for f in realignedNTeq/*
do
	sample=$(echo $f | rev | cut -d/ -f1 | rev)
	java \
		-Xmx7g \
		-jar ${macseFolder}/macse_v2.03.jar \
		-prog reportMaskAA2NT \
		-align_AA "trimmedAA/"$sample"_masked" \
		-align $f \
		-mask_AA $ \
		-out_NT trimmedNT/$sample \
		-out_mask_detail "trimmedNT/"$sample"_mask" \
		-dist_isolate_AA 4
done


## translate NT to obtain new AA
python ${repo_folder}/python_scripts/translator.py ./trimmedNT/ -t 1.0
```

### 6.4. Remove flanking sequence artefacts
As previous steps were not sensitive enough to remove flanking artefactual sites, a custom end-trimming script was used to mask such sites. MACSE was then run to trim these sites in these sequences, followed by translation to obtain a new set of CDS and AA sequences.
```
#end trimming based on AA alignments
mv translated trimmedNT_translated

for f in trimmedNT_translated/*.fas
do
	echo $f
	python ${repo_folder}/endbite.py $f
done

# mask NT using AA as template and
# translate to get new AA

mkdir trimmedNTa
for f in trimmedNT/*.fas
do
	sample=$(echo $f | rev | cut -d/ -f1 | rev)
	java \
		-Xmx7g \
		-jar ${macseFolder}/macse_v2.03.jar \
		-prog reportMaskAA2NT \
		-align_AA "trimmedNT_translated/"$sample".masked" \
		-align $f \
		-mask_AA $ \
		-out_NT trimmedNTa/$sample \
		-out_mask_detail "trimmedNTa/"$sample"_mask" \
		-dist_isolate_AA 4
done

python ${repo_folder}/translator.py ./trimmedNTa/ -t 1.0
mv translated trimmedNTa_translated
```

## 7. AA Gene tree based filtering
### 7.1. Infer AA gene trees
Create a list of files
```
mkdir gt_filt

ls trimmedNTa_translated/ > gt_filt/locilist.txt

date

```
RAxML was run as an array job on the list of files to efficiently parallelize on many small files
```
${raxml} \
	-T 8 \
	--asc-corr lewis \
	-f a \
	-m PROTGAMMALG \
	-p 12345 \
	-x 12345 \
	-n $(sed -n ${SLURM_ARRAY_TASK_ID}p locilist.txt) \
	-N 50 \
	-s ../trimmedNTa_translated/$(sed -n ${SLURM_ARRAY_TASK_ID}p locilist.txt)
```


### 7.2. Run long branch trimmer
Removes long branch outliers (over 50 times mean edge length, parameter determined empirically). Run twice to trim both NT and AA files.
```
> lbfilt_AA.log
for f in trimmedNTa_translated/*.fas
do
	echo $f
	Rscript ${repo_folder}/R_scripts/lbfilter.R \
		"gt_filt/RAxML_bipartitions."$(basename ${f}) \
		$f \
		prot \
		50 >> lbfilt_AA.log
done

> lbfilt_NT.log
for f in trimmedNTa/*.fas
do
	echo $f
	Rscript ${repo_folder}/R_scripts/lbfilter.R \
		"gt_filt/RAxML_bipartitions."$(basename ${f}) \
		$f \
		dna \
		50 >> lbfilt_NT.log
done


#move edited files to new location

mkdir gtfilt_trimmedNT_translated
for f in trimmedNTa_translated/*.edited
do
	mv $f ./gtfilt_trimmedNT_translated/$(basename ${f} | cut -f1-2 -d.)
done

mkdir gtfilt_trimmedNT
for f in trimmedNTa/*.edited
do
	mv $f ./gtfilt_trimmedNT/$(basename ${f} | cut -f1-2 -d.)
done
```


## 8. Iterative HMMCleaner and global distance filtering
We ran two iterations of HMMCLeaner and a custom script to remove outlier sequences from alignments. Details on the second approach are in 8.2.
## 8.1. Run HMMCleaner I
```
#run HMMCleaner
for f in gtfilt_trimmedNT_translated/*.fas
do
	echo $f
	HmmCleaner.pl $f -v=0
done

#run transferCleaner
for f in gtfilt_trimmedNT/*.fas
do
	sample=$(basename ${f} | cut -f1 -d.)
	echo $sample
	transferCleaner.pl $f -log="gtfilt_trimmedNT_translated/"$sample"_hmm.log"
done

#run ali2fasta
ali2fasta.pl gtfilt_trimmedNT/*_cleaned.ali 

mkdir hmmclean_gtfilt_trimmedNT_translated


# move output files to new location
for f in gtfilt_trimmedNT_translated/*_hmm.fasta
do
	mv $f hmmclean_gtfilt_trimmedNT_translated/$(basename ${f} | cut -f1 -d_).fas
done

mkdir hmmclean_gtfilt_trimmedNT

for f in gtfilt_trimmedNT/*_cleaned.fasta
do
	mv $f hmmclean_gtfilt_trimmedNT/$(basename ${f} | cut -f1 -d_).fas
done
```


### 8.2. Trim missing data and prepare files for distance filter
```
# remove seqs with only 25% data or less
python ${repo_folder}/python_scripts/removeTaxa.py hmmclean_gtfilt_trimmedNT/ -l 0.25
mv rmtaxaout hmmclean_gtfilt_trimmedNT_l025

# remove loci with less than 40 samples
python ${repo_folder}/python_scripts/removeTaxa.py hmmclean_gtfilt_trimmedNT_l025/ -ll 40
mv rmtaxaout hmmclean_gtfilt_trimmedNT_l025_ll40

# remove samples, removed from NT, from AA
python ${repo_folder}/python_scripts/removeTaxa.py \
	./hmmclean_gtfilt_trimmedNT_translated/ \
	-m \
	./hmmclean_gtfilt_trimmedNT_l025_ll40/
mv rmtaxaout hmmclean_gtfilt_trimmedNT_translated_l025_ll40

# get a concat NT matrix, no partition file needed
python ${repo_folder}/python_scripts/concat.py hmmclean_gtfilt_trimmedNT_l025_ll40 -1
mv COMBINED.phy concatNT.phy
rm partitions.prt

# convert to fasta (default output is phylip for phylogenetics but not appropriate for distfilter)
python ${repo_folder}/python_scripts/fconv.py \
	-a concatNT.phy phylip-relaxed fasta .fas

# get a concat AA matrix
python ${repo_folder}/python_scripts/concat.py hmmclean_gtfilt_trimmedNT_translated_l025_ll40 -1
mv COMBINED.phy concatAA.phy
rm partitions.prt

# conver to fasta
python ${repo_folder}/python_scripts/fconv.py \
	-a \
	concatAA.phy \
	phylip-relaxed \
	fasta \
	.fas

# move files around
mkdir distfilt1
mv concat* distfilt1/
```


### 8.3. Run distfilter I
Distance outliers in each locus are determined based on how average pairwise distance of a sample in a NT alignment far from average pairwise distance of a sample in a concat NT matrix. Same is checked for protein data and protein concat file. For NT raw distance was used, while for AA an ML distance under WAG was calculated. Union of outliers (outlier in NT or AA, or both) is taken and sample sequences removed from NT and AA. First pass outlier threshold settings are IRQ\*5 for the upper (distance too large compared to concat) and IRQ\*10 for the lower (distance too small compared to concat), were determined empirically.
```
cd distfilt1

ln -s ../hmmclean_gtfilt_trimmedNT_l025_ll40 prefiltNT

ln -s ../hmmclean_gtfilt_trimmedNT_translated_l025_ll40 prefiltAA

Rscript ${repo_folder}/R_scripts/distfilter.R \
	concatNT.phy.fas \
	./prefiltNT/ \
	concatAA.phy.fas \
	./prefiltAA/ \
	5 \
	10 \
	rescale

mkdir distfiltNT

for f in prefiltNT/*.edited
do
	mv $f ./distfiltNT/$(basename ${f} | cut -f1-2 -d.)
done

mkdir distfiltAA

for f in prefiltAA/*.edited
do
	mv $f ./distfiltAA/$(basename ${f} | cut -f1-2 -d.)
done

```


### 8.4. Run HMMCleaner II
Second iteration of HMMCleaner, since effect was marginal, no additional trimming of missing data was done.
```
for f in distfiltAA/*.fas
do
	echo $f
	HmmCleaner.pl $f -v=0
done

for f in distfiltNT/*.fas
do
	sample=$(basename ${f} | cut -f1 -d.)
	echo $sample
	transferCleaner.pl $f -log="distfiltAA/"$sample"_hmm.log"
done

ali2fasta.pl distfiltNT/*_cleaned.ali 

mkdir hmmclean_distfiltAA

for f in distfiltAA/*_hmm.fasta
do
	mv $f hmmclean_distfiltAA/$(basename ${f} | cut -f1 -d_).fas
done

mkdir hmmclean_distfiltNT

for f in distfiltNT/*_cleaned.fasta
do
	mv $f hmmclean_distfiltNT/$(basename ${f} | cut -f1 -d_).fas
done
```


### 8.5. Run distfilter II
Second iteration of distfilter. Second pass outlier threshold settings are IRQ\*2 for the upper (distance too large compared to concat, more sensitive than in the first pass) and IRQ\*10 for the lower (distance too small compared to concat), were determined empirically.
```
mkdir ../distfilt2/

cd ../distfilt2/

ln -s ../distfilt1/hmmclean_distfiltNT prefiltNT

ln -s ../distfilt1/hmmclean_distfiltAA prefiltAA

python ${repo_folder}/python_scripts/concat.py prefiltNT -1

mv COMBINED.phy concatNT.phy

rm partitions.prt

python ${repo_folder}/python_scripts/fconv.py \
	-a \
	concatNT.phy \
	phylip-relaxed \
	fasta \
	.fas

python ${repo_folder}/python_scripts/concat.py prefiltAA -1

mv COMBINED.phy concatAA.phy

rm partitions.prt

python ${repo_folder}/python_scripts/fconv.py \
	-a \
	concatAA.phy \
	phylip-relaxed \
	fasta \
	.fas

Rscript ${repo_folder}/R_scripts/distfilter.R \
	concatNT.phy.fas \
	./prefiltNT/ \
	concatAA.phy.fas \
	./prefiltAA/ \
	2 \
	10 \
	rescale

mkdir distfiltNT

for f in prefiltNT/*.edited
do
	mv $f ./distfiltNT/$(basename ${f} | cut -f1-2 -d.)
done

mkdir distfiltAA

for f in prefiltAA/*.edited
do
	mv $f ./distfiltAA/$(basename ${f} | cut -f1-2 -d.)
done
```


## 9. NT gene tree based filtering
Several filtering steps based on NT trees
### 9.1. Run concat and gene tree analyses of NT alignments
Infer gene trees using RAxML
```
${raxml} \
	-T 8 \
	--asc-corr lewis \
	-f a \
	-m GTRGAMMA \
	-p 12345 \
	-x 12345 \
	-n $(sed -n ${SLURM_ARRAY_TASK_ID}p locilist.txt) \
	-N 50 \
	-s ../distfilt2/distfiltNT_la50/$(sed -n ${SLURM_ARRAY_TASK_ID}p locilist.txt)
```
Concatenate all alignments and infer NT ML concat trees using IQ-TREE. Model test was done separately as (at least in early versions of IQ-TREE) missing data in outgroup caused abort during model test.
```
#concat distfilter2 results
${iqtree} \
	-nt AUTO \
	-s an1NT.phy \
	--prefix an1NTtest \
	-m MF \
	-spp an1NT.prt

${iqtree} \
	-nt AUTO \
	-s an1NT.phy \
	-alrt 1000 \
	-bb 1000 \
	--prefix an1NTprt \
	-o Megro.CLC.fasta \
	-spp an1NTtest.best_model.nex
```


### 9.2. Run decrosscontamination
As opposed to a more simple algorithm in distfilter (which does check for a too small distance), cross contamination script more carefully checks for near identical sequence similarity, and tries to remove only one sequence (contaminant acceptor) if this can be inferred based on the overall ML tree. First, suspect contaminated groups of taxa are discovered based on pairwise GTR distance of less than 0.01 in a locus given a concat-based GTR distance of over 0.2 (pair is too similar to each other despite too far in the concat tree), parameter choince was determined empirically. This preserves similar sequences for several congeneric taxa in our analysis.

Then we try to rescue at least one of the suspect sequences. Since gene tree error with respect to species tree can be large and a simple check of ML-based distance to the ML tree (like in distfilter) can be misleading, we instead try to drop each sequence and check if RF distance of the resulting tree to ML tree reduced. A sample which when dropped improves RF similarity more is removed. If two tests give indentical RF, both sequences are removed. 
```
#collapse by support
for f in gt_filtAn1/RAxML_bipartitions.L*
do
	echo $f
	Rscript ${repo_folder}/R_scripts/cmd_collapseN.R $f 33
done

mkdir gt_treesAn1

mv gt_filtAn1/*.newick gt_treesAn1/

# concat gene trees
cat gt_treesAn1/* > genetreesAn1.tre

ls gt_treesAn1/* | cut -f2-3 -d. > genetreesAn1.txt

# run decrosscont
Rscript ${repo_folder}/R_scripts/crosscontTree.R \
	an1NT/an1NTprt_scf.cf.tree genetreesAn1.tre \
	genetreesAn1.txt \
	./distfilt2/distfiltNT_la50/ \
	./distfilt2/distfiltAA_la50/ \
	0.01 \
	0.2

mkdir crosscontNT

for f in ./distfilt2/distfiltNT_la50/*.edited
do
	mv $f ./crosscontNT/$(basename ${f} | cut -f1-2 -d.)
done

mkdir crosscontAA

for f in ./distfilt2/distfiltAA_la50/*.edited
do
	mv $f ./crosscontAA/$(basename ${f} | cut -f1-2 -d.)
done

mv genetreesAn1.tre.edited genetreesAn1_crosscont.tre
```


### 9.3. RF-based filtering of entire loci
After cross contamination is addressed, updated gene trees are checked for RF distance to the concat ML distance and outlier loci (with distance over IRQ\*3) are removed.
```
# run rfFilter
Rscript ${repo_folder}/R_scripts/rfFilter.R \
	genetreesAn1_crosscont.tre \
	genetreesAn1.txt \
	an1NT/an1NTprt_scf.cf.tree

cp -r crosscontNT rf_filt

rm rf_filt/*.reduced

cp -r crosscontAA rf_filt_translated

while read l
do
	echo $l
	rm rf_filt/$l
	rm rf_filt_translated/$l
done < rf_loci.txt 
```


## 10. spruceup filtering
### 10.1. Run concatenation analysis on NT alignment
```
#concat rffilter results
${iqtree} \
	-nt AUTO \
	-s an2NT.phy \
	--prefix an2NTtest \
	-m MF \
	-spp an2NT.prt

${iqtree} \
	-nt AUTO \
	-s an2NT.phy \
	-alrt 1000 \
	-bb 1000 \
	--prefix an2NTprt \
	-o Megro.CLC.fasta \
	-spp an2NTtest.best_model.nex \
	-bsam GENESITE
```


### 10.2. Run spruceup

```
python -m spruceup config.conf
```


## 11. Final missing data and locus filtering
### 11.1. Filter missing data
```
#split file
#missing data 0.25
#missing taxa 40
```


### 11.2. Remove loci with PLS outliers
TBA
