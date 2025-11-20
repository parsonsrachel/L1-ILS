#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=250gb
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --qos=highmem
#SBATCH --time=4-00:00:00

# paths
WORKDIR="/fs/cbcb-lab/ekmolloy/rparsons/primates_T2T_short_v_long/final"
CHAINPATH="/fs/cbcb-lab/ekmolloy/group/data/t2t_primate/data/chains/chm13#1.p70.aln.chain.gz"
BIGBEDTOBED="/fs/cbcb-lab/ekmolloy/group/software/ucsc-linux/bigBedToBed"
LIFTOVER="/fs/cbcb-lab/ekmolloy/group/software/ucsc-linux/liftOver"
SEGDUP="/fs/cbcb-lab/ekmolloy/rparsons/primates_T2T_short_v_long/data/SegmentalDuplications"
BLAST="/fs/cbcb-lab/ekmolloy/group/software/ncbi-blast-2.14.1+-src/c++/ReleaseMT/bin"
MAFFT="/fs/cbcb-lab/ekmolloy/group/software/mafft-7.487-with-extensions/install/bin/mafft"
IQTREE="/fs/cbcb-lab/ekmolloy/group/software/iqtree-2.2.2.6-Linux/bin/iqtree2"

# download data
# mkdir $WORKDIR/repeat_masker_out
# wget -c -P $WORKDIR/repeat_masker_out https://hgdownload.soe.ucsc.edu/hubs/GCA/009/914/755/GCA_009914755.4/bbi/GCA_009914755.4_T2T-CHM13v2.0.t2tRepeatMasker/chm13v2.0_rmsk.bb
# $BIGBEDTOBED $WORKDIR/repeat_masker_out/chm13v2.0_rmsk.bb $WORKDIR/repeat_masker_out/chm13v2.0.bed # convert bigBed to bed

# wget -c -P $WORKDIR/genomes https://s3.amazonaws.com/genomeark/species/Gorilla_gorilla/mGorGor1/assembly_curated/mGorGor1.dip.cur.20231122.fasta.gz
# wget -c -P $WORKDIR/genomes https://s3.amazonaws.com/genomeark/species/Pan_paniscus/mPanPan1/assembly_curated/mPanPan1.dip.cur.20231122.fasta.gz
# wget -c -P $WORKDIR/genomes https://s3.amazonaws.com/genomeark/species/Pan_troglodytes/mPanTro3/assembly_curated/mPanTro3.dip.cur.20231122.fasta.gz
# wget -c -P $WORKDIR/genomes https://s3.amazonaws.com/genomeark/species/Pongo_abelii/mPonAbe1/assembly_curated/mPonAbe1.dip.cur.20231205.fasta.gz
# wget -c -P $WORKDIR/genomes https://s3.amazonaws.com/genomeark/species/Pongo_pygmaeus/mPonPyg2/assembly_curated/mPonPyg2.dip.cur.20231122.fasta.gz
# wget -c -P $WORKDIR/genomes https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz

# gunzip $WORKDIR/genomes/mGorGor1.dip.cur.20231122.fasta.gz 
# gunzip $WORKDIR/genomes/mPanPan1.dip.cur.20231122.fasta.gz
# gunzip $WORKDIR/genomes/mPanTro3.dip.cur.20231122.fasta.gz
# gunzip $WORKDIR/genomes/mPonAbe1.dip.cur.20231205.fasta.gz
# gunzip $WORKDIR/genomes/mPonPyg2.dip.cur.20231122.fasta.gz
# gunzip $WORKDIR/genomes/hs1.fa.gz

# split the chain file into species specific files
# mkdir -p $WORKDIR/chains/all_species
# mkdir $WORKDIR/chains/by_spec_haplotype
# gunzip -c $CHAINPATH > $WORKDIR/chains/all_species/chm13#1.p70.aln.chain
# python3 $WORKDIR/scripts/a_split_chain.py $WORKDIR

# perform filtering on TEs and write to format to use to liftover
# mkdir $WORKDIR/te_beds 
# mkdir $WORKDIR/call_set
# python3 $WORKDIR/scripts/b_te_beds.py -i $WORKDIR/repeat_masker_out/chm13v2.0.bed -o1 $WORKDIR/te_beds -o2 $WORKDIR/call_set -r L1PA2 L1PA3 L1PA4

# look for the TSD and count the number of full length elements
# python3 $WORKDIR/scripts/c_tsd_and_len.py $WORKDIR

# perform liftover
# chmod +x $WORKDIR/scripts/d_liftover.sh
# mkdir -p $WORKDIR/lifted_beds/L1PA{2,3,4}
# $WORKDIR/scripts/d_liftover.sh $WORKDIR $LIFTOVER

# add lifted coordinates to call set
# mkdir -p $WORKDIR/coords_to_pull/L1PA{2,3,4}
# python3 $WORKDIR/scripts/e_cand_orth.py $WORKDIR $SEGDUP

# isolate the sequences from each genome to then align
# mkdir -p $WORKDIR/sequences
# python3 $WORKDIR/scripts/f_get_from_genome.py $WORKDIR

# get alignment of sequences and to orfs
# mkdir -p $WORKDIR/alignments
# mkdir -p $WORKDIR/orf
# python3 $WORKDIR/scripts/g_align.py $WORKDIR $MAFFT $BLAST
# python3 $WORKDIR/scripts/g_wfmash_align.py $WORKDIR/wfmash $WORKDIR/call_set

# simulate gene trees
# mkdir -p $WORKDIR/trees/te
# mkdir -p $WORKDIR/trees/te_f1_f2
# mkdir -p $WORKDIR/trees/f1
# mkdir -p $WORKDIR/trees/f2
# python3 $WORKDIR/scripts/h_sim_gene_trees.py $WORKDIR $IQTREE 

# process the data
# mkdir -p $WORKDIR/tsd
# mkdir -p $WORKDIR/alignment_data
# python3 $WORKDIR/scripts/i_process.py $WORKDIR
# python3 $WORKDIR/scripts/i_process_wfmash.py $WORKDIR/wfmash

mkdir $WORKDIR/figures
python3 $WORKDIR/scripts/k_star_stats.py $WORKDIR
