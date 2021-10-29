#!/bin/bash



wbp_prefix="ftp://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/"

species=Phylogenetics/Tocris/parasite.list.txt

# -N: only download newer versions
# -nc: no clobber; ignore server files that aren't newer than the local version
# -r: recursive
# -nH: don't mimick the server's directory structure
# -cut-dirs=7: ignore everything from pub to species in the recursive search
# --no-parent: don't ascend to the parent directory during a recursive search
# -A: comma-separated list of names to accept
# -P:

mkdir input/proteomes

while IFS= read -r line
do
  species_dl="$wbp_prefix/$line/"
  printf ${species_dl}"\n"
  wget -nc -r -nH --cut-dirs=9 --no-parent --reject="index.html*" -A 'protein.fa.gz' $species_dl -P input/proteomes
done <"$species"

## Get IDs and sequences of hits
# while IFS= read -r line; do
#  	printf '%s\n' "$line"
#   	echo "$line" > temp.line.txt
#  	line_sub=$(echo "$line" | awk 'BEGIN { FS = "|" } ; { print $3 }')
#  	seqtk subseq human_db/HsUniProt_nr.fasta temp.line.txt > Hs_seeds/Hs_seeds.$line_sub.fasta
#  	rm temp.line.txt
#
#  	#blast seed to human proteome to expand targets
# 	blastp -query Hs_seeds/Hs_seeds.$line_sub.fasta -db human_db/HsUniProt_nr.fasta -out Hs_blast/$line_sub.out -outfmt 6 -max_hsps 1 -evalue 1E-3 -num_threads 4
# 	cat Hs_blast/$line_sub.out | awk '$3>50.000 && $11<1E-3 {print $2}' | sort | uniq  > Hs_blast/$line_sub.list.txt
# 	seqtk subseq human_db/HsUniProt_nr.fasta Hs_blast/$line_sub.list.txt >  Hs_targets/$line_sub.ext.fasta
#
# 	cat Hs_targets/$line_sub.ext.fasta | sed 's/>/>Homo_sapiens|/g' > alignments/$line_sub.combined.fasta
# 	while IFS= read -r paradb; do
# 		#blast expanded human targets against parasite dbs
# 		para_name=$(echo "$paradb" | awk 'BEGIN { FS = "." } ; { print $1 }')
# 		echo "$para_name"
# 		blastp -query Hs_targets/$line_sub.ext.fasta -db parasite_db/$paradb -out Para_blast/$line_sub.$para_name.out -outfmt 6 -max_hsps 1 -evalue 1E-1 -num_threads 4
# 		cat Para_blast/$line_sub.$para_name.out | awk '$3>30.000 && $11<1E-3 {print $2}' | sort | uniq  > Para_blast/$line_sub.$para_name.list.txt
# 		seqtk subseq parasite_db/$paradb Para_blast/$line_sub.$para_name.list.txt >  Para_targets_1/$line_sub.$para_name.fasta
# 		#blast parasite hits against human db
# 		blastp -query Para_targets_1/$line_sub.$para_name.fasta -db human_db/HsUniProt_nr.fasta -out Recip_blast/$line_sub.$para_name.out -outfmt 6 -max_hsps 1 -evalue 1E-3 -num_threads 4
# 		cat Recip_blast/$line_sub.$para_name.out | awk '$3>30.000 && $11<1E-3 {print $1, $2}' | sort | uniq  > Recip_blast/$line_sub.$para_name.list.txt
# 		#compare to original human list to find surviving parasite targets
# 		grep -Ff Hs_blast/$line_sub.list.txt Recip_blast/$line_sub.$para_name.list.txt | awk '{print $1}' | sort | uniq > Para_targets/$line_sub.$para_name.list.txt
# 		seqtk subseq parasite_db/$paradb Para_targets/$line_sub.$para_name.list.txt > Para_targets/$line_sub.$para_name.fasta
# 		cat Para_targets/$line_sub.$para_name.fasta | sed 's/>/>'$para_name'|/g' >> alignments/$line_sub.combined.fasta
#
# 		#align mafft
# 		einsi --reorder --thread 2 alignments/$line_sub.combined.fasta > alignments/$line_sub.combined.aln
# 		#trim
# 		trimal -gt 0.7 -in alignments/$line_sub.combined.aln -out alignments/$line_sub.combined_trim.aln
# 		trimal -resoverlap 0.70 -seqoverlap 70 -in alignments/$line_sub.combined_trim.aln -out alignments/$line_sub.combined_final.aln
#
# 		#tree-building
# 		iqtree-2.1.3-MacOSX/bin/iqtree2 -s alignments/$line_sub.combined_final.aln -nt 4 -alrt 1000 -bb 1000
# 	done < parasite_db/parasite_db.list.txt
#
# done < Hs_seeds/Hs_seeds.list.txt
