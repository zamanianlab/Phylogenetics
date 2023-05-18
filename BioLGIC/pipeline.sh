#!/bin/bash

line_sub=$1

wbp_prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species"

species=Phylogenetics/BioLGIC/parasite.list.txt


# dowload parasite proteomes
mkdir input/proteomes
proteomes=input/proteomes

while IFS= read -r line
do
  species_dl="$wbp_prefix/$line/"
  printf ${species_dl}"\n"
  wget -nc -r -nH --cut-dirs=10 --no-parent --reject="index.html*" -A 'protein.fa.gz' $species_dl -P $proteomes
done <"$species"


# make blast databases
while IFS= read -r line
do
  species_prjn="$(echo $line | sed 's/\//\./g')"
  echo $species_prjn > work/temp.txt
  gunzip -k $proteomes/"$species_prjn".*.protein.fa.gz
  makeblastdb -in $proteomes/"$species_prjn".*.protein.fa -dbtype prot
done <"$species"
rm $proteomes/*.protein*.gz



# set up directories
cd work
mkdir Ce_seeds Ce_targets alignments Para_targets Para_recip Para_final
cd ../

Ce_seeds=work/Ce_seeds
Ce_targets=work/Ce_targets
Para_targets=work/Para_targets
Para_recip=work/Para_recip
Para_final=work/Para_final
alignments=work/alignments
threads=4


# copy Ce list files (Biogenic amine LGIC seq ids)
mv Phylogenetics/BioLGIC/Ce_BLGIC.list.txt $Ce_seeds


# get sequences of seeds
seqtk subseq $proteomes/caenorhabditis_elegans.PRJNA13758.WBPS18.protein.fa $Ce_seeds/Ce_BLGIC.list.txt > $Ce_seeds/Ce_seed."$line_sub".fasta


# blast Ce seed to Ce proteome to expand targets
blastp -query $Ce_seeds/Ce_seed."$line_sub".fasta -db $proteomes/caenorhabditis_elegans.PRJNA13758.WBPS18.protein.fa -out $Ce_targets/"$line_sub".out -outfmt "6 qseqid sseqid pident evalue qcovs" -max_hsps 1 -evalue 1E-3 -num_threads $threads
cat $Ce_targets/"$line_sub".out | awk '$3>30.000 && $4<1E-4 && $5>40.000 {print $2}' | sort | uniq > $Ce_targets/"$line_sub".list.txt
seqtk subseq $proteomes/caenorhabditis_elegans.PRJNA13758.WBPS18.protein.fa $Ce_targets/"$line_sub".list.txt > $Ce_targets/"$line_sub".ext.fasta

#cat $Ce_targets/"$line_sub".ext.fasta | sed 's/>/>caenorhabditis_elegans|/g' > $alignments/"$line_sub".combined.fasta

while IFS= read -r paradb; do
    #blast expanded Ce targets against parasite dbs
    para_name=$(echo "$paradb" | awk 'BEGIN { FS = "." } ; { print $1 }')
    blastp -query $Ce_targets/"$line_sub".ext.fasta -db $proteomes/$paradb -out $Para_targets/"$line_sub"."$para_name".out -outfmt "6 qseqid sseqid pident evalue qcovs" -max_hsps 1 -evalue 1E-1 -num_threads $threads
    cat $Para_targets/"$line_sub"."$para_name".out | awk '$3>30.000 && $4<1E-4 && $5>40.000 {print $2}' | sort | uniq > $Para_targets/"$line_sub"."$para_name".list.txt
    seqtk subseq $proteomes/$paradb $Para_targets/"$line_sub"."$para_name".list.txt > $Para_targets/"$line_sub"."$para_name".fasta
    #blast parasite hits against Ce db
    blastp -query $Para_targets/"$line_sub"."$para_name".fasta -db $proteomes/caenorhabditis_elegans.PRJNA13758.WBPS18.protein.fa -out $Para_recip/"$line_sub"."$para_name".out -outfmt "6 qseqid sseqid pident evalue qcovs" -max_hsps 1 -evalue 1E-3 -num_threads $threads
    cat $Para_recip/"$line_sub"."$para_name".out | awk '$3>30.000 && $4<1E-4 && $5>40.000 {print $1, $2}' | sort | uniq  > $Para_recip/"$line_sub"."$para_name".list.txt
    #compare to original Ce list to find surviving parasite targets
    grep -Ff $Ce_targets/"$line_sub".list.txt $Para_recip/"$line_sub"."$para_name".list.txt | awk '{print $1}' | sort | uniq > $Para_final/"$line_sub"."$para_name".list.txt
    seqtk subseq $proteomes/$paradb $Para_final/"$line_sub"."$para_name".list.txt > $Para_final/"$line_sub"."$para_name".fasta
    cat $Para_final/"$line_sub"."$para_name".fasta | sed 's/>/>'$para_name'|/g' >> $alignments/"$line_sub".combined.fasta
done < Phylogenetics/BioLGIC/parasite_db.list.txt

#align mafft
einsi --reorder --thread $threads $alignments/"$line_sub".combined.fasta > $alignments/"$line_sub".combined.aln
#trim
trimal -gt 0.7 -in $alignments/"$line_sub".combined.aln -out $alignments/"$line_sub".combined_trim.aln
trimal -resoverlap 0.70 -seqoverlap 70 -in $alignments/"$line_sub".combined_trim.aln -out $alignments/"$line_sub".combined_final.aln
#tree-building
iqtree -s $alignments/"$line_sub".combined_final.aln -nt AUTO -alrt 1000 -bb 1000 -redo

mv $alignments output