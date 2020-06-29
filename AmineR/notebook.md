# Notebook for Amine GPCR comparative genomics and phylogenetic analysis

## Daily Notebook

### June 25, 2020

1. Initialized new repo that will hold all phylogenetic pipelines (zamanianlab/Phylogenetics)
2. This repo will contain directories of gene families (i.e., AmineR)
3. All analysis should be performed in modular format, so scripts can be copied/pasted/edited between gene family projects
4. The entire pipeline should be written to `ZamanianLab/Data/Genomics/{gene family}`
5. I decided to put all intermediary files (i.e., BLAST output, and different transcript_id lists for each filtering step) in `ZamanianLab/Data/Genomics/{gene family}/identification_pipeline`
    - when the pipeline is complete, new diretories will be created in the root `{gene_family}` directory to hold the final list of transcript_ids, amino acid/CDS/mRNA FASTA files, filtered GTFs, etc.

- Can't forget to incorporate PacBio data where appropriate (I believe one or two gene models may have been corrected)
- Pfam families:
  - 7tm_1 (PF00001)
  - Subfamily A17 = 5HT (2/6), adrenergic, dopamine, trace amine, histamine (H2)
  - Subfamily A18 = Histamine (H1, H3, H4), adenosine, muscarinic acetylcholine
  - Subfamily A19 = 5HT (1/5/7)
  - Missing = tyramine, octopamine

### June 29, 2020

- In the end, I realized that a complex HMM/BLAST based pipeline was no better than simply pulling orthologs/paralogs from WBP (written in R in get_homologues.R)


## Pipeline

### Identification
1. Start with list of Ce aminergic GPCRs (`aux/celegans_aminer.txt`) and pull all WBP orthologs
2. Use the list resulting from (1) as seeds for a new ortholog pull
3. Use the list resulting from (2) as seeds for a paralog pull
4. Combine all 3 lists and remove duplicates

### Alignment
1. Predict TMs
2. Filter sequences that don't have 5-9 TMs
3. Extract sequence from first to last TM +/- 1 nt
4. Concatenate sequences and align with mafft
