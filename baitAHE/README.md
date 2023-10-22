# AHE bait pipeline

Procedures to reorganize the original (short) AHE target set into a modified (long and with paralogs and duplicates removed) phylogenomic target locus set.

## 1. Check the three hemipteran AHE target region sets to be orthologous

**Original locus count: 478**

### 1.1. Re BLAST the original AHE target regions against transcriptomes

Given that the relevant original AHE regions used belonged to Arilus, Phymata, as well as an auchenorrhynchan taxon, we searched back these regions against Arilus transcriptome (to be used as reference) to verify all previously constructed AHE regions are orthologous.

### 1.2. Pullout using ALiBaSeq

Only pull out the best and all suboptimal hits for each bait:
```
python alibaseq.py -b baitX_vs_transcriptomeY.blast -f S -x b -t transcriptomeY_assembly.fasta -o out_baitX_vs_transcriptomeY -s log_baitX_vs_transcriptomeY -e 1e-10 -c=-1 --ac tdna-tdna
```

Then combine per locus pull out resuts (i.e., file L100.fas from `out_baitX_vs_transcriptomeY` gets concatenated with `out_baitZ_vs_transcriptomeY`).

## 1.3. Run MAFFT to align bait and best match

Align concatenated locus files from the previous step with MAFFT

## 1.4. Run checkdist to figure out distance

Run `check_dist.py` on the aligned files to detect paralogs. The script checks raw pairwise distance excluding gaps and N, and splits original loci into paralogous new loci if the distance is above 5%.

For example,
```
python check_dist.py L100.fas
```

**Cluster split locus count: 504**

## 2. Combine (consensus merge) overlapping loci

### 2.1. Self blast
Combine all clustered sequences from the output of checkdist script into a single file, and run blastn to search for loci having hits to other loci. Locate such hits by checking different query and target name (e.g., `awk '{if ($1 != $2) print $0}'`) and split this subset of blast output into a separate file.

### 2.2. Merge loci with hits to other loci
Now run `merge_overlapped.py` script on the subset blast result to get a list of combinations of overlapped loci, which then can be combined into new merged the locus files.
Then align the obtain files with MAFFT and run `get_consensus.py` to get a combined sequence of the overlapped AHE region.

**Clustered, overlap-merged locus count: 499**

## 3. Combine loci located in the same transcript

Besides combining the original AHE regions that were overlapping, we also opted to combine nearby regions (parts of the same transcript). To do that, we first concatenated all regions from the previous step into one fasta file, and searched that against the Arilus transcriptome. We only searched against the ORF transcripts.

### 3.1. Direct blastn against Arilus TRANSDEC CDS

Run blastn of the clustered merged sequences against the Arilus ORF seqs. To produce ORF sequences, transdecoder was used.

### 3.2. awk reverse blast

The resulting blast output was inverted with awk (to change query and target) and then parsed by alibaseq with a contig stitching mode. This way, AHE regions from the same transcript were being stitched together in the appropriate order according to how they were found in the original transcript.

For example:
```
awk -F $'\t' 'BEGIN {OFS=FS} !_[$1][$3][$4][$8][$9][$10][$11][$12]++ {print $2".fas",$1,$3,$4,$5,$6,$9,$10,$7,$8,$11,$12}' original.blast > reversed.blast
```

### 3.3. parse with alibaseq extract to stitch probes

```
python alibaseq.py -b db_vs_arcriORF_blast_reversed.txt -f S -x b -o merge_probes_out -s log_merge_probes -i 100 -t db_arcri_regions_clustered_merged.fas --is
```

### 3.4 run renamer using reverse logs to rename trans loci to AHE stitched loci

Then run `rename_AHE.py` script to rename the obtained files and sequence names according to the original naming scheme.

**Clustered, overlap-merged, transcript-merged locus count: 402**

## 4. Repeat step 2

Re-check and merge any loci with overlaps to other loci

### 4.1. self blast

After merging by transcript, double check locus overlap. For that, combine the obtained transcript-merged loci and self blast them. As before in step 2.1., `awk '{if ($1 != $2) print $0}'` can be used to track down the non-self hits.

### 4.2 merge loci via consensus

As in step 2.2, combine overlapped loci into a single file, align it and run `get_consensus.py` to produce a combine sequence. 

**Clustered, overlap-merged, transcript-merged, overlap-merged locus count: 400**

## 5. Extend loci to 250bp

We also decided to enlarge the AHE regions by 250bp if transcripts CDS was extending that far in either direction.

### 5.1. run blastn against CDS

We first run blast against the longest ORF again as in 3.1.

### 5.2. and alibaseq to get full seq + 250bp flanks

Then we parse the results with alibaseq to extract 250bp flank if possible:

```
python alibaseq.py -b longest_orfs.cds.blast -f S -x b -o merge_probes_out -s log_merge_probes -t Arcri.Trinity.fasta.transdecoder_longest_orfs.cds --fl 250
```

## 6. Check loci against transcriptomes of 4 taxa

We do a test run of locus retrieval as it would be done later on with other transcriptomic and genomic taxa, but only on 4 high quality diverse transcriptomic taxa. We check that results are translatable across the taxa.

### 6.1. Run tblastx of CDS baits against 4 transcriptomes

### 6.2. Run Alibaseq

We pulled out homologs from the four transcriptomes in `-x a` mode.
```
python alibaseq.py -b ./blasts/ -f M -x a -o albskresults -t transcriptomes/ --amalgamate-hits --ac tdna-tdna -q ./merge_probes_out_fixed/
```

### 6.3. Trim to bait bases

### 6.4. Translate to check frame

### 6.5. Adjust / fix detected artifacts

Several loci were excluded, L415 and L490 frames were corrected, and L512aL513 had its flank manually corrected.

**Final locus count: 397**
