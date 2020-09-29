#!/bin/bash
bin=../../pasta-code/pasta/run_pasta.py
pre=../..
post=rose.unaln.true.fasta

declare -a datasets=("1000M1" "1000M4")

# for 16S.M
if [ ! -d 16S.M ]; then
        mkdir 16S.M
fi
python $bin -i "$pre/16S.M/R0/cleaned.unaln.fasta" -o "16S.M/"


# for 1000M1 and 1000M4
for d in ${datasets[@]}; do
        if [ ! -d $d ]; then
                mkdir $d
        fi
        for i in $(seq 0 9); do
                if [ ! -d $d/R$i ]; then
                        mkdir $d/R$i
                fi
                
                # run pasta for one replicate, output to corresponding
                # directory
                python $bin -i "$pre/$d/$d/R$i/$post" -o $d/R$i/
        done
done
