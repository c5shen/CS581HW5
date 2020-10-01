#!/bin/bash
bin=../FastSP/FastSP.jar
pre=../..
post=rose.aln.true.fasta
name=pastajob.marker001.rose.unaln.true.aln

declare -a datasets=("1000M1" "1000M4")

# 16S.M
java -Xmx4096m -jar $bin -r $pre/16S.M/R0/cleaned.alignment.fasta \
        -e 16S.M/pastajob.marker001.cleaned.unaln.aln -o 16S.M/fastsp.out

# for 1000M1 and 1000M4
for d in ${datasets[@]}; do
        for i in $(seq 0 9); do
                java -Xmx4096m -jar $bin -r $pre/$d/$d/R$i/$post \
                        -e $d/R$i/$name -o $d/R$i/fastsp.out
        done
done
