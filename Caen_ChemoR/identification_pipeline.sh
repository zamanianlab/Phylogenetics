#!/bin/bash

### Preparation
repo="Phylogenetics"
proj="Caen_ChemoR"

gh_dir="$GIT_PATH/$repo/$proj"

# for local MacOS
local_dir="$HOME/Box/ZamanianLab/LabMembers/Nic/project-$proj"

# for server
# local_dir="~/WBP"

## Species list
species="${gh_dir}"/aux/species.txt

# mkdir "$local_dir/identification_pipeline"
out="$local_dir/identification_pipeline"
db="$local_dir/db"

## misc
linearize="${gh_dir}"/scripts/aux/linearizefasta.awk

# HMM approach -----------------------------------------------------------------

## Combined (GRAFS+ / ChemoR)
GRAFS_ChemoR_HMM="${gh_dir}"/aux/GRAFS_ChemoR.hmm

## Prepare Pfam-A HMM db
# wget -nc -O $gh_dir/aux/Pfam-A.hmm.gz ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
# gzcat $gh_dir/aux/Pfam-A.hmm.gz > $gh_dir/aux/Pfam-A.hmm
# hmmpress $gh_dir/aux/Pfam-A.hmm

#############
# FILTER ONE
#############

# mkdir $out/1/

### Mine for nematode ChemoRs
# while IFS= read -r line; do
#   for f in $db/$line/*.fa* ; do
#     curr_dir=$(dirname $f)
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     gzcat -f $f > $curr_dir/protein.tmp.fa
#     #HMMSEARCH all proteomes against db of All GPCR hmms
#     hmmsearch --tblout $out/1/${species}_1.out --noali $GRAFS_ChemoR_HMM $curr_dir/protein.tmp.fa
#     rm $curr_dir/protein.tmp.fa
#   done;
# done <$species

### Parse hmm outputs to filter out those where first hit is not ChemoR or 7tm_1 HMM, extract sequences of surviving hits
# while IFS= read -r line; do
#   for f in $db/$line/*.fa* ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     cat $out/1/${species}_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g > $out/1/${species}_1.txt
#     cat $out/1/${species}_1.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '!/Frizzled|7tm_2|7tm_3|7tm_4|7tm_6|7tm_7/' | sort -k4 -g | awk '{print $1}' > $out/1/${species}_1_ids.txt
#     curr_dir=$(dirname "${f}")
#     gzcat -f $f > $curr_dir/protein.tmp.fa
#     seqtk subseq $curr_dir/protein.tmp.fa $out/1/${species}_1_ids.txt > $out/1/${species}_1.fa
#     rm $curr_dir/protein.tmp.fa
#     done;
# done <$species

#############
# FILTER TWO
#############

# mkdir $out/2/

### Reciprocal HMMSEARCH of extracted sequences against Pfam-a
# while IFS= read -r line; do
#   for f in $db/$line/*.fa* ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     hmmsearch --tblout $out/2/${species}_2.out --noali --cpu 4 $gh_dir/aux/Pfam-A.hmm $out/1/${species}_1.fa
#   done;
# done <$species

### Parse hmm outputs to remove sequences where first hit is not ChemoR or 7tm_1 HMM, get list of surviving unique IDs, extract sequences
# while IFS= read -r line; do
#   for f in $db/$line/*.fa* ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     cat $out/2/${species}_2.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1|7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g > $out/2/${species}_2.txt
#     cat $out/2/${species}_2.out | awk '{print $1 " " $3 " " $4 " " $5}' | awk '!/#/' | sort -k1,1 -k4,4g | sort -uk1,1 | awk '/7tm_1|7TM_GPCR_S|Srg|Sre|Serpentine_r_xa/' | sort -k4 -g  | awk '{print $1}' > $out/2/${species}_2_ids.txt
#     curr_dir=$(dirname "${f}")
#     gzcat -f $f > $curr_dir/protein.tmp.fa
#     seqtk subseq $curr_dir/protein.tmp.fa $out/2/${species}_2_ids.txt > $out/2/${species}_2.fa
#     rm $curr_dir/protein.tmp.fa
#     # compare the list from _1 to the list from _2 and write out the IDs that were removed
#     ggrep -v -f $out/2/${species}_2_ids.txt $out/1/${species}_1_ids.txt > $out/2/${species}_2_filtered.txt
#   done;
# done <$species

##############
# FILTER THREE
##############

# mkdir $out/3/

### Reciprocal blastp of extracted sequences against C. elegans
# while IFS= read -r line; do
#   for f in $db/$line/*.fa* ; do
#     array=($(echo "$line" | sed 's/\// /g'))
#     species=${array[0]}
#     # blast filtered ChemoRs against C. elegans proteome, using E-value cutoff
#     cd $db/c_elegans/
#     blastp -query $out/2/${species}_2.fa -db c_elegans.PRJNA13758.WS280.protein.fa -out $out/3/${species}_3.out -num_threads 4 -evalue 0.01 -outfmt '6 qseqid sseqid pident ppos length mismatch evalue bitscore'
#   done;
# done <$species

### Run json_parser.py to populate list of C. elegans ChemoR paralogs to use a comparison for the reciprocal BLAST output
# python "${query_api}"
## use WormBase SimpleMine to convert Gene_IDs to Sequence_IDs (becauase BLASTp uses Sequence_IDs)
## upload list of Gene_IDs that result from the above command
## choose "allow duplicate genes" and "Transcript"
## in Sublime, remove transcript numbers (e.g. N.1) and put each isoform on a new line
## 2177 transcripts

### Remove hits that aren't most similar to a C. elegans ChemoR, extract sequences of surviving hits
while IFS= read -r line; do
  for f in $db/$line/*.fa* ; do
    array=($(echo "$line" | sed 's/\// /g'))
    species=${array[0]}
    ## Find genes that had no hit at all
    cat $out/3/${species}_3.out | awk '{print $1}' | ggrep -wFv -f -  $out/2/${species}_2_ids.txt > $out/3/${species}_3_nohits.txt
    ## Find genes that had a significant hit to a C. elegans ChemoR
    cat $out/3/${species}_3.out | awk '{print $1 " " $2 " " $7}' | sort -k1,1 -k3,3g | sort -uk1,1 | ggrep -wF -f "${gh_dir}"/aux/simplemine_results.txt | sort -k3 -g > $out/3/${species}_3.txt
    cat $out/3/${species}_3.txt | awk '{print $1}' | cat - $out/3/${species}_3_nohits.txt > $out/3/${species}_3_ids.txt
    curr_dir=$(dirname "${f}")
    gzcat -f $f > $curr_dir/protein.tmp.fa
    seqtk subseq $curr_dir/protein.tmp.fa $out/3/${species}_3_ids.txt > $out/3/${species}_3.fa
    rm $curr_dir/protein.tmp.fa
    ggrep -v -f $out/3/${species}_3_ids.txt $out/2/${species}_2_ids.txt > $out/3/${species}_3_filtered.txt
  done;
done <$species
