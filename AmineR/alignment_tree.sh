#!/bin/bash

### Preparation
repo="Phylogenetics"
proj="AmineR"

gh_dir="$GIT_PATH/$repo/$proj"

# for local MacOS
# local_dir="$HOME/Box/ZamanianLab/Data/Genomics/$proj"

# for server
local_dir="$GIT_DATA/$repo/$proj"

## Species list
species="${gh_dir}"/aux/species_selected.txt

db="$HOME/Box/ZamanianLab/Data/WBP"
out="$local_dir/phylo"

## HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
HMMTOP_py="$gh_dir/aux/scripts/HMMTOP_extract.py"
HMMTOP_strict_py="$gh_dir/aux/scripts/HMMTOP_extract_strict.py"

## misc
# linearize="${gh_dir}"/scripts/aux/linearizefasta.awk

## copy FASTA files to new directory
# cp $out/../identification_pipeline/3/*_1.fa $out/1/

# count TM domains and extract TM-only sequences for anything with 5-9 TMs -----

## TM count and sequence extract
# while IFS= read -r line; do
#   for f in $db/$line/*.protein.fa.gz ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     cd $gh_dir/aux/scripts/hmmtop_2.1/
#     ./hmmtop -if=$out/1/${species}_1.fa -of=$out/2/${species}_2_hmmtop.txt -sf=FAS
#     python $HMMTOP_py $out/2/${species}_2_hmmtop.txt $out/1/${species}_1.fa $out/2/${species}_2.fa
#   done;
# done <$species

## label each sequence with its species name
# while IFS= read -r line; do
#   for f in $db/$line/*.protein.fa.gz ; do
#     for f in $out/2/${species}_2.fa ; do
#       array=($(echo "$line" | sed 's/\// /g'))
#       species=${array[0]}
#       python $gh_dir/aux/scripts/id_change.py $f
#     done;
#   done;
# done <$species

## NOTE: id_change.py always messes up the Loa loa IDs, so manuaally changed those

# concatenate and align --------------------------------------------------------

# cat $out/2/*_2_label.fa > $out/3/all_aminer_tm.fa

## align (see https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html)
## for algorithm description
# einsi --reorder --thread 8 $out/3/all_aminer_tm.fa > $out/3/all_aminer_tm_1.aln

## trim and filter
# trimal -gt 0.7 -in $out/3/all_aminer_tm_1.aln -out $out/3/all_aminer_tm_2.aln
# trimal -resoverlap 0.70 -seqoverlap 70 -in $out/3/all_aminer_tm_2.aln -out $out/3/all_aminer_tm_3.aln

## change to single-line FASTA
# awk -f $gh_dir/aux/scripts/linearizefasta.awk < $out/3/all_aminer_tm_2.aln > $out/3/all_aminer_tm_2_linear.aln
# awk -f $gh_dir/aux/scripts/linearizefasta.awk < $out/3/all_aminer_tm_3.aln > $out/3/all_aminer_tm_3_linear.aln

## Get IDs and compare lists
# cat $out/3/all_aminer_tm_2_linear.aln | tr '\t' '\n' | awk 'NR%2==1' > $out/3/all_aminer_tm_2_ids.txt
# cat $out/3/all_aminer_tm_3_linear.aln | tr '\t' '\n' | awk 'NR%2==1' > $out/3/all_aminer_tm_3_ids.txt
# grep -v -f $out/3/all_aminer_tm_3_ids.txt $out/3/all_aminer_tm_2_ids.txt > $out/3/all_aminer_3_filtered.txt

## ML tree on server
iqtree -s $out/3/all_aminer_tm_3.aln -nt 4 -alrt 1000 -bb 1000
