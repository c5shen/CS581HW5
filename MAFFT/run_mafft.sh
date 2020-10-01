#!/bin/bash
pre=../..
post=rose.unaln.true.fasta

declare -a datasets=("1000M1" "1000M4")

# for 16S.M
if [ ! -d 16S.M ]; then
        mkdir 16S.M
fi
mafft --thread -1 $pre/16S.M/R0/cleaned.unaln.fasta > 16S.M/mafft.aln.fasta


# for 1000M1 and 1000M4
for d in ${datasets[@]}; do
        if [ ! -d $d ]; then
                mkdir $d
        fi
        for i in $(seq 0 9); do
                if [ ! -d $d/R$i ]; then
                        mkdir $d/R$i
                fi
                
                # run mafft for one replicate, output to corresponding
                # directory
                mafft --thread -1 $pre/$d/$d/R$i/$post > $d/R$i/mafft.aln.fasta
        done
done
