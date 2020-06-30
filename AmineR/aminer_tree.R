library(tidyverse)
library(treeio)
library(tidytree)
library(ggtree)
library(cowplot)
library(conflicted)
library(here)

conflicted::conflict_prefer("filter", "dplyr")

# import tree -------------------------------------------------------------

iqtree <- read.newick(here('..', 'all_aminer_tm_3.aln.contree'))

t1 <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none")

# unrooted ----------------------------------------------------------------

unrooted_data <- t1$data %>%
  separate(label, into = c("species", "transcript_id"), sep = "-", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(species = case_when(
    isTip == FALSE ~ as.character(NA),
    isTip == TRUE ~ species,
  )) 

# get protein ID's of C. elegans genes
library(biomaRt)

mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

ce_proteins <- getBM(mart = mart, 
                     filters = c("wbps_transcript_id"),
                     value = list((filter(unrooted_data, species == 'celeg') %>% 
                                     dplyr::select(transcript_id) %>% 
                                     mutate(transcript_id = str_c(transcript_id, '.1')))$transcript_id),
                     attributes = c("wbps_peptide_id", "external_gene_id")) %>%
  janitor::clean_names() %>%
  dplyr::select(transcript_id = 1, protein_id = 2) %>%
  mutate(transcript_id = str_remove(transcript_id, '.[0-9]$'))

detach('package:biomaRt', unload = TRUE)

unrooted_data <- left_join(unrooted_data, ce_proteins) %>%
  mutate(protein_id = str_remove(protein_id, ", isoform.*")) %>%
  mutate(label2 = case_when(
    is.na(protein_id) ~ transcript_id, 
    !is.na(protein_id) ~ protein_id
  )) 

t2 <- t1 %<+% unrooted_data +
  geom_tiplab(aes(label = label2)) +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

save_plot(here("t1_unrooted.pdf"), t2, base_height = 20)

# rooted ------------------------------------------------------------------

# reroot on the mgl family, node found by examining the unrooted tree
rooted <- phytools::reroot(iqtree, 566)

t3 <- ggtree(rooted, size = 1.3, layout = "circular", branch.length = "none")

rooted_data <- t3$data %>%
  separate(label, into = c("species", "transcript_id"), sep = "-", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(species = case_when(
    isTip == FALSE ~ as.character(NA),
    isTip == TRUE ~ species,
  )) %>%
  left_join(., ce_proteins) %>%
  mutate(label2 = case_when(
    is.na(protein_id) ~ transcript_id, 
    !is.na(protein_id) ~ protein_id
  ))  

t3 <- t3 %<+% rooted_data +
  geom_tiplab(aes(label = label2)) +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

save_plot(here("t2_rooted.pdf"), t3, base_height = 20)

# Collapse isoforms -------------------------------------------------------

phylo <- as.phylo(t3)

# get the list of tips for every internal node (will filter to only include nodes that have only a single species in tip offspring)
get_offspring_species <- function(node, phylo) {
  data <- phylo$data
  offspring <- offspring(data, node, tiponly = TRUE)
  unique(offspring$species)
}

# count the number of tips from that node
count_offspring_species <- function(node, phylo) {
  data <- phylo$data
  offspring <- offspring(data, node, tiponly = TRUE)
  length(offspring$species)
}

# get list of all nodes, check to see if the offspring() are the same species, and count the number
offspring <- rooted_data %>%
  filter(isTip == FALSE) %>%
  dplyr::select(node) %>%
  mutate(offspring_species = map(node, ~ get_offspring_species(.x, t3))) %>%
  filter(!str_detect(string = offspring_species, pattern = 'c\\(')) %>%
  mutate(offspring_species = unlist(offspring_species)) %>%
  mutate(n_tips_of_species = map(node, ~ count_offspring_species(.x, t3))) %>%
  mutate(n_tips_of_species = unlist(n_tips_of_species)) 

# node selection
collapse <- tibble(collapse = c('octr-1', 'gar-3', 'gar-1', 'gar-2', 'tyra-3', 'ser-6', 
                                'tyra-2', 'ser-2', 'dop-1', 'ser-7', 'ser-1', 'dop-5',
                                'dop-3', 'dop-2'),
                   node = c(790, 770, 731, 743, 662, 686, 
                            637, 614, 584, 569, 519, 504,
                            479, 466))

t4 <- ggtree(rooted, size = 0.25, layout = "circular", branch.length = "none") %<+% collapse

# collapse node with loop
for (i in 1:length(collapse$node)) {
  t4 <- ggtree::collapse(t4, collapse$node[i])
}

# get the tree data from the newly collapsed tree, label the collapsed nodes 
collapsed_data <- t4$data %>%
  separate(label, into = c("species", "transcript_id"), sep = "-", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(species = case_when(
    isTip == FALSE ~ as.character(NA),
    isTip == TRUE ~ species,
  )) %>%
  mutate(bootstrap = case_when(
    !is.na(transcript_id) ~ as.numeric(NA),
    species == "Root" ~ as.numeric(NA),
    is.na(transcript_id) ~ as.numeric(label)
  )) %>%
  left_join(., ce_proteins) %>%
  mutate(label = case_when(
    is.na(protein_id) ~ transcript_id, 
    !is.na(protein_id) ~ protein_id
  )) %>%
  left_join(., collapse)

t5 <- t4 %<+% collapsed_data +
  geom_tiplab(aes(label = str_c(species, transcript_id, str_replace_na(protein_id, ''), sep = '-'))) +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

save_plot(here("t3_collapsed.pdf"), t5, base_height = 20)

# family annotations ------------------------------------------------------

t6 <- t4 %<+% collapsed_data +
  geom_tiplab(aes(label = str_c(species, transcript_id, str_replace_na(protein_id, ''), sep = '-'))) +
  geom_label2(aes(subset = !isTip, label = bootstrap, fill = bootstrap)) +
  geom_text2(aes(subset = !isTip, label = node), nudge_x = 1, nudge_y = 1) +
  scale_fill_viridis_c() +
  NULL

# save_plot(here("t4_bootstrap.pdf"), t6, base_height = 20)

# choose nodes to have displayed bootstrap levels
bootstrap_label <- tribble(
  ~bootstrap_family, ~node,
  "Class C", 797,
  "octr-1", 784,
  "ador−1", 774,
  "gar-3", 752,
  "gar-1", 719,
  "gar-2", 733,
  "gar", 717,
  "tyra-3", 645,
  "ser-6", 679,
  "ser-3", 668,
  "Orphan", 713,
  "pcdr-1", 689,
  "tyra-2", 625,
  "ser-2", 602,
  "ser-4", 587,
  "dop-1", 571,
  "ser-7", 550,
  'ser-5', 534,
  "dop−4", 520,
  "ser-1", 505,
  "dop-5", 492,
  "dop-6", 483,
  "dop-3", 474,
  "dop-2", 446
)

# assign colors for species tip labels
species <- tibble(species = c('avite', 'bmala', 'btimo', 'bpaha', 'ooche', 'ovolv', 'asuum', 'tcani', 'bxylo', 'predi', 'scarp', 'aceyl', 'hcont', 'hpoly', 'celeg', 'ptric', 'sratt', 'tmuri', 'tspir', 'puniv', 'ppaci', 'dimmi', 'oflex', 'dmedi', 'namer', 'lsigm', 'lloa', 'eelap'),
                  clade = c('IIIc', 'IIIc', 'IIIc', 'IIIc', 'IIIc', 'IIIc', 'IIIb', 'IIIb', 'IV', 'IV', 'IV', 'V', 'V', 'V', 'V', 'IV', 'IV', 'I', 'I', 'IIIb', 'V', 'IIIc', 'IIIc', 'IIIb', 'V', 'IIIc', 'IIIc', 'IIIc'),
                  full = c('*A. viteae*', '*B. malayi*', '*B. timori*', '*B. pahangi*', '*O. ochengi*', '*O. volvulus*', '*A. suum*', '*T. canis*', '*B. xylophilus*', '*P. redivivus*', '*S. carpocapsae*', '*A ceylanicum*', '*H. contortus*', '*H. polygyrus*', '*C. elegans*', '*P. trichosuri*', '*S. ratti*', '*T. muris*', '*T. spiralis*', '*P. univalens*', '*P. pacificus*', '*D. immitis*', '*O. flexuosa*', '*D. medinensis*', '*N. americanus*', '*L. sigmodontis*', '*L. loa*', '*E. elaphi*'),
)

bootstrap_data <- collapsed_data %>%
  left_join(., bootstrap_label) %>%
  mutate(branch_type = case_when(
    !is.na(transcript_id) ~ species,
    is.na(transcript_id) ~ "inner"
  )) %>%
  left_join(., species) %>%
  mutate(isTip = case_when(
    !is.na(family) ~ TRUE,
    TRUE ~ isTip
  ))

t7 <- t4 %<+% bootstrap_data + 
  geom_tiplab(aes(color = clade), label = "-", size = 10) +
  geom_tippoint(aes(subset = !is.na(collapse)), color = "grey", size = 3) +
  geom_label2(aes(subset = !is.na(bootstrap_family), label = bootstrap), fill = 'grey30', color = "white", size = 2) +
  scale_color_manual(values = c(I = "#ABB065", IIIb = "#E495A5", IIIc = "dark blue", IV = "#39BEB1", V = "#ACA4E2"), name = "Clade") +
  NULL

save_plot(here("t5_family.pdf"), t7, base_height = 8, limitsize = FALSE)

t8 <- t7 +
  geom_cladelabel(797, offset.text = 2, hjust = 0, offset = 2, barsize = 2, fontsize = 4, label = "Class C (Outgroup)") +
  geom_cladelabel(784, offset.text = 2, hjust = 0.5, offset = 2, barsize = 2, fontsize = 4, label = "octr-1") +
  geom_cladelabel(774, offset.text = 2, hjust = 0.5, offset = 2, barsize = 2, fontsize = 4, label = "ador-1") +
  geom_cladelabel(752, offset.text = 2, hjust = 0.5, offset = 2, barsize = 2, fontsize = 4, label = "gar-3") +
  geom_cladelabel(719, offset.text = 2, hjust = 0.5, offset = 2, barsize = 2, fontsize = 4, label = "gar-1") +
  geom_cladelabel(733, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "gar-2") +
  geom_cladelabel(645, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "tyra-3") +
  geom_cladelabel(679, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "ser-6") +
  geom_cladelabel(668, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "ser-3") +
  geom_cladelabel(689, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "pcdr-1") +
  geom_cladelabel(625, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "tyra-2") +
  geom_cladelabel(602, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "ser-2") +
  geom_cladelabel(587, offset.text = 2, hjust = 1, offset = 2, barsize = 2, fontsize = 4, label = "ser-4") +
  geom_cladelabel(571, offset.text = 2, hjust = 0.75, offset = 2, barsize = 2, fontsize = 4, label = "dop-1") +
  geom_cladelabel(550, offset.text = 2, hjust = 0.75, offset = 2, barsize = 2, fontsize = 4, label = "ser-7") +
  geom_cladelabel(534, offset.text = 2, hjust = 0.5, offset = 2, barsize = 2, fontsize = 4, label = "ser-5") +
  geom_cladelabel(520, offset.text = 2, hjust = 0.5, offset = 2, barsize = 2, fontsize = 4, label = "dop-4") +
  geom_cladelabel(505, offset.text = 2, hjust = 0.25, offset = 2, barsize = 2, fontsize = 4, label = "ser-1") +
  geom_cladelabel(492, offset.text = 2, hjust = 0, offset = 2, barsize = 2, fontsize = 4, label = "dop-5") +
  geom_cladelabel(483, offset.text = 2, hjust = 0, offset = 2, barsize = 2, fontsize = 4, label = "dop-6") +
  geom_cladelabel(474, offset.text = 2, hjust = 0, offset = 2, barsize = 2, fontsize = 4, label = "dop-3") +
  geom_cladelabel(446, offset.text = 2, hjust = 0, offset = 2, barsize = 2, fontsize = 4, label = "dop-2") +
  NULL

save_plot(here("t6_clade.pdf"), t8, base_height = 8)
  
# new clades --------------------------------------------------------------

new_family_label <- tribble(
  ~new_family, ~node,
  "meloidogyne_enriched", 3066,
  "bursaphelenchus_enriched", 2974,
  "bursaphelenchus_enriched", 3453,
  "bursaphelenchus_enriched", 3502,
  "bursaphelenchus_enriched", 3454,
  "ppn_enriched", 4282
)

t9 <- t8 %<+% new_family_label +
  geom_label2(aes(subset = collapse_family %in% c("sri", "srsx", "srj"), label = collapse_family), fill = "grey10", color = "white", size = 5, nudge_x = -5, nudge_y = -7) +
  geom_label2(aes(subset = !is.na(new_family), label = bootstrap, fill = bootstrap), color = "white", size = 5) +
  geom_cladelabel(3066, offset.text = 2, offset = 2, barsize = 5, fontsize = 8, label = "Meloidogyne+", color = "#c96d44") +
  geom_cladelabel(3452, offset.text = 2, offset = 2, barsize = 5, fontsize = 8, label = "Bursaphelenchus+", color = '#c65999') +
  geom_cladelabel(4282, offset.text = 2, offset = 2, barsize = 5, fontsize = 8, label = "PPN+") +
  # geom_cladelabel(4148, offset.text = 2, offset = 2, barsize = 5, fontsize = 8, label = "", color = 'blue') +
  labs(color = "Genus", fill = "Bootstrap Support") +
  theme(legend.title = element_text(face = "bold", size = 16),
        legend.position = c(0.97, 0.77),
        legend.text = ggtext::element_markdown(size = 14),
        legend.key.size = unit(2, "line")
        ) +
  guides(colour = guide_legend(override.aes = list(size = 10, stroke = 10))) +
  NULL
# t9

save_plot(here("ggtree/t7_newfamily.pdf"), t9, base_height = 15, limitsize = FALSE)
