#!/bin/bash

### Preparation
repo="Phylogenetics"
proj="Caen_ChemoR"

# for server
gh_dir="$repo/$proj"\

## Species list
species="${gh_dir}"/aux/species.txt

db="$proj/db"
out="$proj/phylo"

## HMMTOP parsing script (filter based on TM range and produce sequences based on TM domains)
# HMMTOP_py="${gh_dir}"/aux/scripts/HMMTOP_extract.py
HMMTOP_strict_py="${gh_dir}"/aux/scripts/HMMTOP_extract_strict.py

## misc
linearize="${gh_dir}"/aux/scripts/linearizefasta.awk

# count TM domains and extract TM-only sequences for anything with 3-10 TMs -----
## copy FASTA files to new directory
mkdir $out/1/
cp $out/../identification_pipeline/3/*.fa $out/1/

mkdir $out/2/

## TM count and sequence extract
while IFS= read -r line; do
  for f in $db/$line/*.fa* ; do
    array=($(echo "$line" | sed 's/\// /g'))
    species=${array[0]}
    cd $gh_dir/aux/scripts/hmmtop_2.1/
    ./hmmtop -if=$out/1/${species}_3.fa -of=$out/1/${species}_1_hmmtop.txt -sf=FAS
    python $HMMTOP_strict_py $out/1/${species}_1_hmmtop.txt $out/1/${species}_3.fa $out/2/${species}_2.fa
  done;
done <$species

## label each sequence with its species name
for f in $out/2/*_2.fa; do
  python $gh_dir/aux/scripts/id_change.py $f
done

### Remove Ce
rm $out/2/caen*_label.fa

# concatenate and align --------------------------------------------------------

mkdir $out/3/

### Choose one or more representatives from each clade (5074 sequences)
cat $out/2/*_label.fa > \
  $out/3/caen_1.fa

cp $gh_dir/aux/23.aln $out/3/

### add non-C.elegans representatives to the alignment
mafft --reorder --thread 8 --addfull $out/3/caen_1.fa --keeplength $out/3/23.aln > $out/3/caen_1.aln
### trim and filter
trimal -gt 0.7 -in $out/3/caen_1.aln -out $out/3/caen_2.aln
trimal -resoverlap 0.70 -seqoverlap 70 -in $out/3/caen_2.aln -out $out/3/caen_3.aln
### Change to single-line FASTA
awk -f $linearize < $out/3/caen_2.aln > $out/3/caen_2_linear.aln
awk -f $linearize < $out/3/caen_3.aln | sed '/^$/d' > $out/3/caen_3_linear.aln
### Get IDs and compare lists
cat $out/3/caen_2_linear.aln | tr '\t' '\n' | awk 'NR%2==1' > $out/3/caen_2_ids.txt
cat $out/3/caen_3_linear.aln | tr '\t' '\n' | awk 'NR%2==1' > $out/3/caen_3_ids.txt
ggrep -v -f $out/3/caen_3_ids.txt $out/3/caen_2_ids.txt > $out/3/caen_3_filtered.txt

### ML tree on server
iqtree -s $out/3/caen_3.aln -nt 8 -alrt 1000 -bb 1000
