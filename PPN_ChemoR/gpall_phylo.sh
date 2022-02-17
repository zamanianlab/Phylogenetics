#!/bin/bash

### Preparation
repo="Phylogenetics"
proj="PPN_ChemoR"

gh_dir="$GIT_PATH/$repo/$proj"

# for local MacOS
local_dir="/Users/njwheeler/Library/CloudStorage/Box-Box/ZamanianLab/Data/Genomics/$proj"

# for server
# local_dir="$GIT_DATA/$repo/$proj"

## Species list
species="${gh_dir}"/aux/species.txt

db="$HOME/Box/ZamanianLab/Data/WBP"
out="$local_dir/gpallida"


# mkdir "${out}"/1/

## misc
linearize="${gh_dir}"/aux/scripts/linearizefasta.awk

## copy FASTA file to new directory
# cp $out/../identification_pipeline/3/globodera_pallida_3.fa $out/1/

## align
# einsi --reorder --thread 8 $out/1/globodera_pallida_3.fa > $out/1/globodera_pallida_1.aln

## trim and filter
# trimal -automated1 -in $out/1/globodera_pallida_1.aln -out $out/1/globodera_pallida_2.aln
trimal -resoverlap 0.50 -seqoverlap 50 -in $out/1/globodera_pallida_2b.aln -out $out/1/globodera_pallida_3.aln

## tree
iqtree2 -s $out/1/globodera_pallida_2.aln -T AUTO -alrt 1000 -bb 1000
