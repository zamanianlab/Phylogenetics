#!/bin/bash

### Preparation
repo="Phylogenetics"
proj="AmineR"

gh_dir="$GIT_PATH/$repo/$proj"

# for local MacOS
local_dir="$HOME/Box/ZamanianLab/Data/Genomics/$proj"

# for server
# local_dir="~/WBP"

## Species list
species="${gh_dir}"/aux/species_selected.txt

db="$HOME/Box/ZamanianLab/Data/WBP"
out="$local_dir/identification_pipeline"

## misc
# linearize="${gh_dir}"/scripts/aux/linearizefasta.awk

# HMM approach -----------------------------------------------------------------

#############
# FILTER ONE
#############

# hmmpress $gh_dir/aux/HMM/7tm_1.hmm

### HMMsearch all proteomes against 7tm_1
# while IFS= read -r line; do
#   for f in $db/$line/*.protein.fa.gz ; do
#     curr_dir=$(dirname $f)
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     echo $curr_dir
#     gzcat $f > $curr_dir/protein.tmp.fa
#     hmmsearch --tblout $out/1/${species}_1.out --noali $gh_dir/aux/HMM/7tm_1.hmm $curr_dir/protein.tmp.fa
#     rm $curr_dir/protein.tmp.fa
#   done;
# done <$species

### Get IDs and sequences of hits
# while IFS= read -r line; do
#   for f in $db/$line/*.protein.fa.gz ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     cat $out/1/${species}_1.out | awk '!/#/' | awk '{print $1 " " $3 " " $4 " " $5}' > $out/1/${species}_1.txt
#     cat $out/1/${species}_1.txt | awk '!/#/' | awk '{print $1}' > $out/1/${species}_1_ids.txt
#     curr_dir=$(dirname "${f}")
#     gzcat $f > $curr_dir/protein.tmp.fa
#     seqtk subseq $curr_dir/protein.tmp.fa $out/1/${species}_1_ids.txt > $out/1/${species}_1.fa
#     rm $curr_dir/protein.tmp.fa
#   done;
# done <$species

#############
# FILTER TWO
#############

# gzcat $gh_dir/aux/HMM/Pfam-A.hmm.gz > $gh_dir/aux/HMM/Pfam-A.hmm
# hmmpress $gh_dir/aux/HMM/Pfam-A.hmm
# rm $gh_dir/aux/HMM/Pfam-A.hmm

### Reciprocal HMMSEARCH of extracted sequences against Pfam-a
# while IFS= read -r line; do
#   for f in $db/$line/*.protein.fa.gz ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     hmmsearch --tblout $out/2/${species}_2.out --noali --cpu 4 $gh_dir/aux/HMM/Pfam-A.hmm $out/1/${species}_1.fa
#   done;
# done <$species

### Keep sequences with 7tm_1 as a top hit; get IDs and sequences
# while IFS= read -r line; do
#   for f in $db/$line/*.protein.fa.gz ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     cat $out/2/${species}_2.out | awk '!/#/' | awk '{print $1 " " $3 " " $4 " " $5}' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1/' | sort -k4 -g > $out/2/${species}_2.txt
#     cat $out/2/${species}_2.txt | awk '!/#/' | awk '{print $1}' > $out/2/${species}_2_ids.txt
#     curr_dir=$(dirname "${f}")
#     gzcat $f > $curr_dir/protein.tmp.fa
#     seqtk subseq $curr_dir/protein.tmp.fa $out/2/${species}_2_ids.txt > $out/2/${species}_2.fa
#     rm $curr_dir/protein.tmp.fa
#   done;
# done <$species

##############
# FILTER THREE
##############

# This pipeline was not completed because the homolog approach (below) was sufficient

# WormBase homolog approach ----------------------------------------------------

# run get_homologues.R interactively to get $out/3/species.txt
while IFS= read -r line; do
  for f in $db/$line/*.protein.fa.gz ; do
    array=($(echo "$line" | sed 's/\// /g'))
    species=${array[0]}
    curr_dir=$(dirname "${f}")
    gzcat $f > $curr_dir/protein.tmp.fa
    seqtk subseq $curr_dir/protein.tmp.fa $out/3/${species}_1_ids.txt > $out/3/${species}_1.fa
    rm $curr_dir/protein.tmp.fa
  done;
done <$species
