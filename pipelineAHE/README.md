# AHE pipeline

Pipeline for processing AHE data


## paths setup

```
repo_folder=""
alibaseq_folder=""
blast=""
raxml=""
iqtree=""
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
	fname=$(echo $f | rev | cut -f1 -d/ | rev)
	echo $fname
	grep -v -f "megaGout/"$fname"_decont_out.txt" $f > "decontblastG/"$fname
done
```


## 3. Run reciprocal search
Output (decontaminated) blast tables are used to check all contigs with matches against the reference transcriptome for a reciprocal best hit (RBH) check which greatly improves orthology prediction accuracy.
```
echo "transcriptomes"
bash ${alibaseq_folder}/reciprocal_search.sh \
	./blastT/ \
	${transcriptomicAssembliesFolder}/ \
	${reciprocalReference} \
	dc-megablast \
	${Nthreads} \
	y \
	${alibaseq_folder}/reciprocal_get_contigs.py \
	${transcriptomicSampleList}

echo "genomes"
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
For large, chromosome-level assemblies, the contigs are too large, so to speed up the search and lower memory requirements, chromosomes are chunked into smaller pieces that are then blasted.
```
bash ${repo_folder}/TBA/reciprocal_chunky_search.sh \
	${referenceAssembliesFolder}/$sample \
	${reciprocalReference} \
	dc-megablast \
	${Nthreads} \
	${repo_folder}/TBA/chunkify.py y
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
	echo "get the reference" $x
	echo $x >> ../exoner_debug.txt
	refname=${baitfolderAA}"/"$(echo $x | rev | cut -d/ -f1 | rev)
	echo "reference" $refname
	echo "splitting" $x
	python ${repo_folder}/python_scripts/split_fasta.py $x
	for f in *
	do
		echo "running exonerate on" $f
		exonerate --model protein2genome -q $refname -t $f -E -n 1 --showalignment false --showvulgar false --verbose 0 --ryo ">%ti\n%tcs" 1> ../exoner_concat/$f 2>> ../exoner_debug.txt
	done
	cat ../exoner_concat/* > ../exoner_out/$(echo $x | rev | cut -d/ -f1 | rev)
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
	outf=$(echo $f | rev | cut -f1 -d/ | rev)
	mafft --globalpair --thread ${Nthreads} --unalignlevel 0.8 $f > realignedAA/$outf
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
		-jar ~/tools/macse_v2.03.jar \
		-prog reportGapsAA2NT \
		-align_AA $f \
		-seq nt_translator/$sample \
		-out_NT realignedNT/$sample
done

```


### 6.3. Alignment block trimming
A custom script was used to do per position trimming of AA alignments at 40% missing data or less, removing sequences with less than 10% data. Then samples that were completely removed from AA were also removed from NT using a custom script, and MACSE was used to mask the remaining CDS. Since MACSE does additional base masking (isolated bases), AA alignment is no longer in sync with CDS after this. So after MACSE trimmed CDS, these files were translated to create new AA files.
```
## trim
python ${repo_folder}/python_scripts/customtrim.py \
	./realignedAA/ \
	-d prot
mv trimmed trimmedAA
python ${repo_folder}/python_scripts/removeTaxa.py \
	./realignedNT/ \
	-m ./trimmedAA/
mv rmtaxaout realignedNTeq
mkdir trimmedNT
for f in realignedNTeq/*
do
	sample=$(echo $f | rev | cut -d/ -f1 | rev)
	java \
		-Xmx7g \
		-jar ~/tools/macse_v2.03.jar \
		-prog reportMaskAA2NT \
		-align_AA "trimmedAA/"$sample"_masked" \
		-align $f \
		-mask_AA $ \
		-out_NT trimmedNT/$sample \
		-out_mask_detail "trimmedNT/"$sample"_mask" \
		-dist_isolate_AA 4
done

## translate NT
python ${repo_folder}/python_scripts/translator.py ./trimmedNT/ -t 1.0
```

### 6.4. Remove flanking sequence artefacts
As previous steps were not sensitive enough to remove flanking artefactual sites, a custom end-trimming script was used to mask such sites. MACSE was then run to trim these sites in these sequences, followed by translation to obtain a new set of CDS and AA sequences.
```
#end trimming
mv translated trimmedNT_translated

for f in trimmedNT_translated/*.fas
do
	echo $f
	python ${repo_folder}/endbite.py $f
done

# mask NT again and retranslate

mkdir trimmedNTa

for f in trimmedNT/*.fas
do
	sample=$(echo $f | rev | cut -d/ -f1 | rev)
	java \
		-Xmx7g \
		-jar ~/tools/macse_v2.03.jar \
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
Run RAxML as an array job on the list of files
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
Removes long branch outliers (over 50 times mean edge length, parameter determined empirically)
```
> lbfilt_AA.log
for f in trimmedNTa_translated/*.fas
do
	echo $f
	Rscript ${repo_folder}/TBA/lbfilter.R \
		"gt_filt/RAxML_bipartitions."$(echo $f | rev | cut -f1 -d/ | rev) \
		$f \
		prot \
		50 >> lbfilt_AA.log
done

> lbfilt_NT.log
for f in trimmedNTa/*.fas
do
	echo $f
	Rscript ${repo_folder}/TBA/lbfilter.R \
		"gt_filt/RAxML_bipartitions."$(echo $f | rev | cut -f1 -d/ | rev) \
		$f \
		dna \
		50 >> lbfilt_NT.log
done

mkdir gtfilt_trimmedNT_translated
for f in trimmedNTa_translated/*.edited
do
	mv $f ./gtfilt_trimmedNT_translated/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done

mkdir gtfilt_trimmedNT
for f in trimmedNTa/*.edited
do
	mv $f ./gtfilt_trimmedNT/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done
```


## 8. Iterative HMMCleaner and global distance filtering
## 8.1. Run HMMCleaner I
```
for f in gtfilt_trimmedNT_translated/*.fas
do
	echo $f
	HmmCleaner.pl $f -v=0
done

for f in gtfilt_trimmedNT/*.fas
do
	sample=$(echo $f | rev | cut -f1 -d/ | rev | cut -f1 -d.)
	echo $sample
	transferCleaner.pl $f -log="gtfilt_trimmedNT_translated/"$sample"_hmm.log"
done

ali2fasta.pl gtfilt_trimmedNT/*_cleaned.ali 

mkdir hmmclean_gtfilt_trimmedNT_translated

for f in gtfilt_trimmedNT_translated/*_hmm.fasta
do
	mv $f hmmclean_gtfilt_trimmedNT_translated/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1 -d_).fas
done

mkdir hmmclean_gtfilt_trimmedNT

for f in gtfilt_trimmedNT/*_cleaned.fasta
do
	mv $f hmmclean_gtfilt_trimmedNT/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1 -d_).fas
done
```


### 8.2. Trim missing data and prepare files for distance filter
```
python ${repo_folder}/python_scripts/removeTaxa.py hmmclean_gtfilt_trimmedNT/ -l 0.25

mv rmtaxaout hmmclean_gtfilt_trimmedNT_l025

python ${repo_folder}/python_scripts/removeTaxa.py hmmclean_gtfilt_trimmedNT_l025/ -ll 40

mv rmtaxaout hmmclean_gtfilt_trimmedNT_l025_ll40

python ${repo_folder}/python_scripts/removeTaxa.py \
	./hmmclean_gtfilt_trimmedNT_translated/ \
	-m \
	./hmmclean_gtfilt_trimmedNT_l025_ll40/

mv rmtaxaout hmmclean_gtfilt_trimmedNT_translated_l025_ll40

python ${repo_folder}/python_scripts/concat.py hmmclean_gtfilt_trimmedNT_l025_ll40 -1

mv COMBINED.phy concatNT.phy

rm partitions.prt

python ${repo_folder}/python_scripts/fconv.py \
	-a concatNT.phy phylip-relaxed fasta .fas

python ${repo_folder}/python_scripts/concat.py hmmclean_gtfilt_trimmedNT_translated_l025_ll40 -1

mv COMBINED.phy concatAA.phy

rm partitions.prt

python ${repo_folder}/python_scripts/fconv.py \
	-a \
	concatAA.phy \
	phylip-relaxed \
	fasta \
	.fas

mkdir distfilt1

mv concat* distfilt1/
```


### 8.3. Run distfilter I
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


```


### 8.4. Run HMMCleaner II
```
mkdir distfiltNT

for f in prefiltNT/*.edited
do
	mv $f ./distfiltNT/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done

mkdir distfiltAA

for f in prefiltAA/*.edited
do
	mv $f ./distfiltAA/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done

for f in distfiltAA/*.fas
do
	echo $f
	HmmCleaner.pl $f -v=0
done

for f in distfiltNT/*.fas
do
	sample=$(echo $f | rev | cut -f1 -d/ | rev | cut -f1 -d.)
	echo $sample
	transferCleaner.pl $f -log="distfiltAA/"$sample"_hmm.log"
done

ali2fasta.pl distfiltNT/*_cleaned.ali 

mkdir hmmclean_distfiltAA

for f in distfiltAA/*_hmm.fasta
do
	mv $f hmmclean_distfiltAA/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1 -d_).fas
done

mkdir hmmclean_distfiltNT

for f in distfiltNT/*_cleaned.fasta
do
	mv $f hmmclean_distfiltNT/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1 -d_).fas
done
```


### 8.5. Run distfilter II
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
	mv $f ./distfiltNT/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done

mkdir distfiltAA

for f in prefiltAA/*.edited
do
	mv $f ./distfiltAA/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done
```


## 9. NT gene tree based filtering
### 9.1. Run concat and gene tree analyses of NT alignments
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
	mv $f ./crosscontNT/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done

mkdir crosscontAA

for f in ./distfilt2/distfiltAA_la50/*.edited
do
	mv $f ./crosscontAA/$(echo $f | rev | cut -f1 -d/ | rev | cut -f1-2 -d.)
done

mv genetreesAn1.tre.edited genetreesAn1_crosscont.tre
```


### 9.3. RF-based filtering of entire loci
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
