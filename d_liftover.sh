#!/bin/bash

WORKDIR=$1
LO=$2
TEBED=('te_beds/L1PA2.bed' 'te_beds/L1PA3.bed' 'te_beds/L1PA4.bed')
SPECIES=('chm13-to-mPanTro3#1.aln.chain' 'chm13-to-mPonAbe1#1.aln.chain' 'chm13-to-mPonPyg2#1.aln.chain' 'chm13-to-mGorGor1#M.aln.chain' 'chm13-to-mGorGor1#P.aln.chain' 'chm13-to-mPanPan1#M.aln.chain' 'chm13-to-mPanPan1#P.aln.chain')
CHAINPATH='chains/by_spec_haplotype/'

for spec in "${SPECIES[@]}"; do
    name="${spec%%.*}"
    for te in "${TEBED[@]}"; do
        base="${te##*/}"
        te_num="${base%%.*}"
        echo "$spec for $te_num"
        $LO $WORKDIR/$te $WORKDIR/$CHAINPATH/$spec $WORKDIR/lifted_beds/$te_num/$name.bed $WORKDIR/lifted_beds/$te_num/$name-unmapped.bed -minMatch=0.08
    done
done
