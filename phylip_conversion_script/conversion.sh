#!/bin/bash

# 16S.M conversion
python3 convert_fasta_to_phy.py ../../16S.M/R0/cleaned.alignment.fasta
python3 convert_fasta_to_phy.py ../mafft/16S.M/mafft.aln.fasta
python3 convert_fasta_to_phy.py ../PASTA/16S.M/pastajob.marker001.cleaned.unaln.aln

# 1000M1 conversion
prefix1=$(pwd)/../../1000M1/1000M1/R
rep=9
for i in $(seq 0 $rep); do
        python3 convert_fasta_to_phy.py $prefix1$i/rose.aln.true.fasta
        python3 convert_fasta_to_phy.py ../PASTA/1000M1/R$i/pastajob.marker001.rose.unaln.true.aln
        python3 convert_fasta_to_phy.py ../mafft/1000M1/R$i/mafft.aln.fasta
done


# 1000M4 conversion
prefix2=$(pwd)/../../1000M4/1000M4/R
for i in $(seq 0 $rep); do
        python3 convert_fasta_to_phy.py $prefix2$i/rose.aln.true.fasta
        python3 convert_fasta_to_phy.py ../PASTA/1000M4/R$i/pastajob.marker001.rose.unaln.true.aln
        python3 convert_fasta_to_phy.py ../mafft/1000M4/R$i/mafft.aln.fasta
done

