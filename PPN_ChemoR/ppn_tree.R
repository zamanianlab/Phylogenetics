library(tidyverse)
library(treeio)
library(tidytree)
library(ggtree)
library(cowplot)
library(conflicted)
library(here)


conflicted::conflict_prefer("filter", "dplyr")

# import tree -------------------------------------------------------------

iqtree <- read.newick(here("ppn_ce_3.aln.contree"))

t1 <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none")

# unrooted ----------------------------------------------------------------

celegans_liftover <- read_csv(here("ggtree/celegans_chemor.csv"), col_names = c("superfamily", "gene_id", "protein_id", "transcript_id", "sequence")) %>%
  slice(2:n()) %>%
  select(-sequence)

unrooted_data <- t1$data %>%
  separate(label, into = c("species", "transcript_id"), sep = "-", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(species = case_when(
    isTip == FALSE ~ as.character(NA),
    isTip == TRUE ~ species,
  )) %>%
  left_join(., select(celegans_liftover, transcript_id, protein_id)) %>%
  mutate(protein_id = str_remove(protein_id, ", isoform.*")) %>%
  mutate(label = case_when(
    is.na(protein_id) ~ transcript_id, 
    !is.na(protein_id) ~ protein_id
  )) 

t2 <- t1 %<+% unrooted_data +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

# save_plot(here("ggtree/t1_unrooted.pdf"), t2, base_height = 20)

# rooted ------------------------------------------------------------------

# reroot on the srw family, node found by examining the unrooted tree
rooted <- phytools::reroot(iqtree, 4520)

t3 <- ggtree(rooted, size = 1.3, layout = "circular", branch.length = "none")

rooted_data <- t3$data %>%
  separate(label, into = c("species", "transcript_id"), sep = "-", remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(species = case_when(
    isTip == FALSE ~ as.character(NA),
    isTip == TRUE ~ species,
  )) %>%
  left_join(., select(celegans_liftover, transcript_id, protein_id)) %>%
  mutate(protein_id = str_remove(protein_id, ", isoform.*")) %>%
  mutate(label = case_when(
    is.na(protein_id) ~ transcript_id, 
    !is.na(protein_id) ~ protein_id
  ))  

t3 <- t3 %<+% rooted_data +
  geom_tiplab() +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

# save_plot(here("ggtree/t2_rooted.pdf"), t3, base_height = 20)

# collapse Ce clades  --------------------------------------------------------------

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

# get list of all nodes, check to see if the offspring() is only Ce, and count the number
offspring <- rooted_data %>%
  filter(isTip == FALSE) %>%
  select(node) %>%
  mutate(offspring_species = map(node, ~ get_offspring_species(.x, t3))) %>%
  filter(!str_detect(string = offspring_species, pattern = 'c\\(')) %>%
  mutate(offspring_species = unlist(offspring_species)) %>%
  mutate(n_tips_of_species = map(node, ~ count_offspring_species(.x, t3))) %>%
  mutate(n_tips_of_species = unlist(n_tips_of_species)) 

# manual node selection
mrca <- tibble(collapse_family = c("sra", "srab", "srb", "srbc", "srbc", 
                                   "srd", "srd", "sre", "srg", "srg",
                                   "srh", "sri", 'srj', 'srm', 'srn',
                                   'sro', 'srr', 'srsx', 'srt', 'srt',
                                   'sru', 'srv', 'srv', 'srw', 'srx', 
                                   'srx', 'srx', 'srx', 'srxa', 'srz', 'str'),
               node = c(2487, 2610, as.numeric(NA), 4220, 4193,
                        3180, 3217, 2706, 3943, 4026, 
                        2780, 3003, 3419, 3498, as.numeric(NA), 
                        as.numeric(NA), 3551, 4598, 3801, 3780, 
                        4052, 4113, 4074, 4766, 3610,
                        3631, 3622, 3660, as.numeric(NA), 3580, 3256))

t4 <- ggtree(rooted, size = 0.25, layout = "circular", branch.length = "none") %<+% mrca

# collapse node with loop
for (i in 1:length(mrca$node)) {
  t4 <- ggtree::collapse(t4, mrca$node[i])
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
  left_join(., select(celegans_liftover, transcript_id, protein_id)) %>%
  mutate(protein_id = str_remove(protein_id, ", isoform.*")) %>%
  mutate(label = case_when(
    is.na(protein_id) ~ transcript_id, 
    !is.na(protein_id) ~ protein_id
  )) %>%
  left_join(., mrca)

t5 <- t4 %<+% collapsed_data +
  geom_tiplab(aes(label = str_c(species, transcript_id, protein_id))) +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

# save_plot(here("ggtree/t3_collapsed.pdf"), t5, base_height = 20)

# family annotations ------------------------------------------------------

t6 <- t4 %<+% collapsed_data +
  geom_tiplab(aes(label = str_c(species, transcript_id, protein_id))) +
  geom_label2(aes(subset = !isTip, label = bootstrap, fill = bootstrap)) +
  geom_text2(aes(subset = !isTip, label = node), nudge_x = 1, nudge_y = 1) +
  scale_fill_viridis_c() +
  NULL

# save_plot(here("ggtree/t4_bootstrap.pdf"), t6, base_height = 20)

# choose nodes to have displayed bootstrap levels
bootstrap_label <- tribble(
  ~bootstrap_family, ~node,
  "srw", 4676,
  "srxa", 4129,
  "srbc", 4159,
  "srsx", 4598,
  "sro", NA,
  "srz", 3560,
  "srv", 4064,
  "sru", 4042,
  "srg", 3941,
  "srr", 3551,
  "srt", 3733,
  "srx", 3601,
  "srb", 2749,
  "sre", 2651,
  "sra/srab", 2485,
  "sra", 2522,
  "sra", 2486,
  "srab", 2567,
  "srn", NA,
  "sri", 3003,
  "srh", 2775,
  "srd", 3161,
  "srm", NA,
  "srj", 3417,
  "str", 3244,
  "bxylo_srh", 2974,
  "bxylo_str", 3502
)

# assign colors for species tip labels
species <- tibble(species = c("cele", "rculi", "bmala", "hcont", "bxylo", "gpall", "grost", "hglyc", "maren", "mflor", "mente", "mgram", "minco", "mhapl", "mjava"),
                  clade = c("V", "I", "III", "V", "IV", "IV", "IV", "IV", "IV", "IV", "IV", "IV", "IV", "IV", "IV"),
                  genus = str_to_sentence(c("*caenorhabditis*", "*romanomermis*", "*brugia*", "*haemonchus*", "*bursaphelenchus*", "*globodera*", "*globodera*", "*heterodera*", "*meloidogyne*", "*meloidogyne*", "*meloidogyne*", "*meloidogyne*", "*meloidogyne*", "*meloidogyne*", "*meloidogyne*")),
                  color = c("grey", "grey", 'grey', "grey", '#c65999', '#7aa456', '#7aa456', "#777acd", "#c96d44", "#c96d44", "#c96d44", "#c96d44", '#c96d44', "#c96d44",  "#c96d44")
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
  geom_tiplab(aes(color = genus), label = "-", size = 5.5) +
  geom_tippoint(aes(subset = !is.na(collapse_family)), color = "grey", size = 3) +
  geom_label2(aes(subset = !is.na(bootstrap_family), label = bootstrap, fill = bootstrap), color = "white", size = 5) +
  scale_fill_viridis_c(limits = c(0, 100), option = "C", direction = -1) +
  scale_color_manual(limits = bootstrap_data$genus,
                     values = bootstrap_data$color) +
  NULL

# save_plot(here("ggtree/t5_family.pdf"), t7, base_height = 15, limitsize = FALSE)

t8 <- t7 +
  geom_cladelabel(4676, offset.text = 2, offset = 2, barsize = 5, fontsize = 7, label = "srw") +
  geom_cladelabel(4129, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srxa") +
  geom_cladelabel(4159, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srbc") +
  # geom_cladelabel(4598, offset.text = 12, offset = 2, barsize = 5, fontsize = 7, label = "srsx") +
  geom_cladelabel(4064, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srv") +
  geom_cladelabel(4042, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "sru") +
  geom_cladelabel(3941, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srg") +
  # geom_cladelabel(3551, offset.text = 9, hjust = -0.5, offset = 2, barsize = 5, fontsize = 7, label = "srr") +
  geom_cladelabel(3733, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srt") +
  geom_cladelabel(3601, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srx") +
  geom_cladelabel(2749, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "srb") +
  geom_cladelabel(2651, offset.text = 2, hjust = 1, offset = 2, barsize = 5, fontsize = 7, label = "sre") +
  geom_cladelabel(2522, offset.text = 2, hjust = 0.3, offset = 2, barsize = 5, fontsize = 7, label = "sra") +
  geom_cladelabel(2486, offset.text = 2, hjust = 0.3, offset = 2, barsize = 5, fontsize = 7, label = "sra") +
  geom_cladelabel(2567, offset.text = 2, hjust = 0.3, offset = 2, barsize = 5, fontsize = 7, label = "srab") +
  geom_cladelabel(3560, offset.text = 2, hjust = 0.3, offset = 2, barsize = 5, fontsize = 7, label = "srz") +
  # geom_cladelabel(3003, offset.text = 24, hjust = 2, offset = 2, barsize = 5, fontsize = 7, label = "sri") +
  geom_cladelabel(2778, offset.text = 2, offset = 2, barsize = 5, fontsize = 7, label = "srh") +
  geom_cladelabel(3161, offset.text = 2, offset = 2, barsize = 5, fontsize = 7, label = "srd") +
  # geom_cladelabel(3417, offset.text = 1 ,hjust = -0.5, offset = 2, barsize = 5, fontsize = 7, label = "srj") +
  geom_cladelabel(3244, offset.text = 4, hjust = 0.75, offset = 2, barsize = 5, fontsize = 7, label = "str") +
  NULL
# t8  

# save_plot(here("ggtree/t6_clade.pdf"), t8, base_height = 15, limitsize = FALSE)

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
