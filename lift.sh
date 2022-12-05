#!/bin/bash

bedin=./data/gwas/TRAF1_region_hg18.bed

# hg19
chain19=./data/chain/hg18ToHg19.over.chain.gz
bedout19=./data/gwas/TRAF1_region_hg19.bed
fail19=./data/gwas/failToLift-hg19.txt

liftOver $bedin $chain19 $bedout19 $fail19

# hg38
chain=./data/chain/hg18ToHg38.over.chain.gz
bedout=./data/gwas/TRAF1_region_hg38.bed
fail=./data/gwas/failToLift.txt

#liftOver $bedin $chain $bedout $fail
