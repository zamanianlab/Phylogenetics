# Phylogenetics

Comparative genomics and phylogenetics pipelines for the Zamanian Lab.

## Usage

Every comparative genomics project should be encapsulated in a separate directory in this repository (i.e., `GitHub/Phylogenetics/AmineR`, `GitHub/Phylogenetics/BetaTubulin`). Each investigator is welcome to organize that directory as they see fit, but generally it should contain the following steps:

  1. Gene identification
      - This includes the entire pipeline for creating a filtered set of related genes.
      - See `AmineR/get_homologues.R` and `PPN_ChemoR/identification_pipeline.sh` for two examples, the first of which uses the WormBase ParaSite orthology calls and the REST API, and the second that uses seed sequences from *C. elegans* and a reciprocal HMM/BLAST approach.
  2. Alignment and tree-building
      - There are a number of potential alignment strategies, and these pipelines should be organized in a single script. The ChemoR pipeline first built a scaffold of HMMs using curated *C. elegans* sequences, and then it used MAFFT to sequentially add parasite genes to this scaffold. The AmineR pipeline trimmed sequences to only include transmembrane regions and aligned these.
      - The multiple sequence alignment probably needs to be trimmed (removing uninformative columns and gap-heavy sequences). This can be performed manually using an editor like Seaview, but it can also be done automatically with [trimal](http://trimal.cgenomics.org/). Note that the parameters used for trimming will need to be heuristically selected.
      - After trimming the alignment to informative regions, the tree is inferred. Two general approaches are preferred, either maximum-likelihood (ML) or Bayesian approaches. For ML, use either [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html) or [IQ-TREE](http://www.iqtree.org/). For Bayesian, use [MrBayes](http://nbisweden.github.io/MrBayes/).
  3. Tree visualization and annotation
      - After finishing the tree, the lab uses the R package [ggtree](https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html) to visualize and annotate the consensus tree. See `AmineR/aminer_tree.R` for an example script.

## A word on organization

It is recommended to include a Markdown file (`notebook.md`) in each project directory that is a running archive of approaches that were developed, attempted, discarded, and eventually selected.

All output files should be written to `/Box/ZamanianLab/Data/Genomics/{Project}`. Files should be separated by species until the required concatentation prior to alignment. It is recommended that the output of each step of the identification pipeline is written to a directory that can be organized in ascending numeric order (see `/Box/ZamanianLab/Data/Genomics/AmineR`) as an example. The output of the tree inference (most likely run on the BRC6 server) should be copied into this directory. It is then recommended to create a separate `ggtree` directory, where output trees can be written in ascending numeric order (see `/Box/ZamanianLab/Data/Genomics/AmineR/phylo/4/ggtree`).
