The entire analysis pipeline can be run by running ./submit.sh after changing the paths at the top of the script. Here we describe the steps being completed in each section of the script in order.


## Download Data

The chain file, genomes, repeat masker annotations, and segmental duplications can be downloaded using the following commands:

Chain file
```
wget -c https://garrisonlab.s3.amazonaws.com/t2t-primates/wfmash-v0.13.0/chains_fixedSiamang/chm13%231.p70.aln.chain.gz
```
Genomes
```
wget -c https://s3.amazonaws.com/genomeark/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.dip.cur.20231122.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.dip.cur.20231122.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.dip.cur.20231122.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.dip.cur.20231205.fasta.gz
wget -c https://s3.amazonaws.com/genomeark/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.dip.cur.20231122.fasta.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz
```
RepeatMasker annotations (BigBed format)
``` 
wget -c https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker/chm13v2.0_rmsk.bb
```
Segmental duplications
```
wget -c https://genomeark.s3.amazonaws.com/species/Pan_troglodytes/mPanTro3/assembly_curated/repeats/mPanTro3_v2.0.SD_v1.0.bb (don't know if these match what I've been using but the ones I have came from Jess)
wget -c https://genomeark.s3.amazonaws.com/species/Pan_paniscus/mPanPan1/assembly_curated/repeats/mPanPan1_v2.0.SD_v1.0.bb
```

## Process the data
1. The chain file contains the mapping to all other species in the alignment, so we split the chain file based on species
2. Convert the big bed file to a bed file.
3. Made bed files for each L1PA TE of interest by filtering and isolating the TE subfamily.
4. Get data about the TSDs for each element.


## LiftOver
5. Used liftOver to get the mapped location in each species for each TE element. We use a minMatch value of 0.08 because the worst case alignment occurs when there are 500 bp of flanking region on either side of an absent TE. We allow the TE to be of length 12,000 to account for potential nested repeats within it. As a result, the maximum total region size is 13,000 bp meaning the maximum gap percentage 12/13 = 0.92 (minimum match in this case is 0.08). 
6. Add the lifted coordinates into the call set
7. Ignore any elements that failed to map to all species for the rest of the analysis.


## Alignment
8. We isolate the flanker + te + flanker region for all elements and align these regions.
9. Here we create alignments using mafft but we also pull out the alignments from wfmash. We found that in downstream analysis the alignments we produced with mafft gave us more data to work with.

## Gene Tree Estimation
10. We simulate gene trees on the alignments using IQTREE.
    
## Process Data
11. As a last step, we add all the analysis results to csv files stored in call_set.
