library(tidyverse)
library(treeio)
library(tidytree)
library(ggtree)
library(cowplot)
library(conflicted)
library(phangorn)
library(ZamanianLabThemes)
library(ggtext)


conflicted::conflict_prefer("filter", "dplyr")

# import tree -------------------------------------------------------------

setwd('/Users/njwheeler/Library/CloudStorage/Box-Box/ZamanianLab/Data/Genomics/PPN_ChemoR/phylo/5/iqtree/')

iqtree <- read.newick("ppn_ce_3.aln.contree")

t1 <- ggtree(iqtree, size = 1.3, layout = "circular", branch.length = "none")

# unrooted ----------------------------------------------------------------

celegans_liftover <- read_csv("ggtree/celegans_chemor.csv", col_names = c("superfamily", "gene_id", "protein_id", "transcript_id", "sequence")) %>%
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

# save_plot("ggtree/t2_rooted.pdf", t3, base_height = 50, limitsize = FALSE)

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
                        2780, 3003, 3419, 3499, as.numeric(NA), 
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
  geom_tiplab(aes(label = str_c(species, transcript_id, protein_id, sep = '|'))) +
  geom_text2(aes(subset = !isTip, label = node)) +
  NULL

# save_plot("ggtree/t3_collapsed.pdf", t5, base_height = 20)

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
  "srd", 3165,
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
    geom_tiplab(aes(color = genus, label = "-"), size = 1.5) +
    geom_tippoint(aes(subset = !is.na(collapse_family)), size = 2, color = 'grey') +
    geom_point2(aes(subset = !is.na(bootstrap_family), fill = bootstrap),
                color = 'white', shape = 21, size = 2) +
    scale_fill_gradientn(colors = c('#d8f3dc', '#b7e4c7', '#95d5b2', '#74c69d', '#52b788', '#40916c', '#2d6a4f', '#1b4332'),
                         limits = c(0, 100)) +
    scale_color_manual(limits = bootstrap_data$genus,
                       values = bootstrap_data$color) +
    NULL

# save_plot("ggtree/t5_family.pdf", t7, base_height = 15, limitsize = FALSE)

t8 <- t7 +
  geom_cladelabel(4676, offset.text = 2, offset = 4, barsize = 2, fontsize = 4, label = "srw") +
  geom_cladelabel(4129, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srxa") +
  geom_cladelabel(4159, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srbc") +
  # geom_cladelabel(4598, offset.text = 12, offset = 4, barsize = 2, fontsize = 4, label = "srsx") +
  geom_cladelabel(4064, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srv") +
  geom_cladelabel(4042, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "sru") +
  geom_cladelabel(3941, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srg") +
  # geom_cladelabel(3551, offset.text = 9, hjust = -0.5, offset = 4, barsize = 2, fontsize = 4, label = "srr") +
  geom_cladelabel(3733, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srt") +
  geom_cladelabel(3601, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srx") +
  geom_cladelabel(2749, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "srb") +
  geom_cladelabel(2651, offset.text = 2, hjust = 1, offset = 4, barsize = 2, fontsize = 4, label = "sre") +
  geom_cladelabel(2522, offset.text = 2, hjust = 0.3, offset = 4, barsize = 2, fontsize = 4, label = "sra") +
  geom_cladelabel(2486, offset.text = 2, hjust = 0.3, offset = 4, barsize = 2, fontsize = 4, label = "sra") +
  geom_cladelabel(2567, offset.text = 2, hjust = 0.3, offset = 4, barsize = 2, fontsize = 4, label = "srab") +
  geom_cladelabel(3560, offset.text = 2, hjust = 0.3, offset = 4, barsize = 2, fontsize = 4, label = "srz") +
  # geom_cladelabel(3003, offset.text = 24, hjust = 2, offset = 4, barsize = 2, fontsize = 4, label = "sri") +
  geom_cladelabel(2778, offset.text = 2, offset = 4, barsize = 2, fontsize = 4, label = "srh") +
  geom_cladelabel(3161, offset.text = 2, offset = 4, barsize = 2, fontsize = 4, label = "srd") +
  # geom_cladelabel(3417, offset.text = 1 ,hjust = -0.5, offset = 4, barsize = 2, fontsize = 4, label = "srj") +
  geom_cladelabel(3244, offset.text = 4, hjust = 0.75, offset = 4, barsize = 2, fontsize = 4, label = "str") +
  NULL
# t8

# save_plot("ggtree/t6_clade.pdf", t8, base_height = 15, limitsize = FALSE)

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
    geom_label2(aes(subset = collapse_family %in% c("sri", "srsx", "srj"), label = collapse_family), fill = "grey10", color = "white", size = 3, nudge_x = -5, nudge_y = -7) +
    geom_point2(aes(subset = !is.na(new_family), fill = bootstrap),
                color = 'white', shape = 21, size = 3) +
    geom_cladelabel(3066, offset.text = 2, offset = 4, barsize = 2, fontsize = 5, label = "Meloidogyne+\n(srd)", color = "#c96d44") +
    geom_cladelabel(3452, offset.text = 2, offset = 4, barsize = 2, fontsize = 5, label = "Bursaphelenchus+\n(srm)", color = '#c65999') +
    geom_cladelabel(4282, offset.text = 2, offset = 4, barsize = 2, fontsize = 5, label = "PPN+") +
    # geom_cladelabel(4148, offset.text = 2, offset = 2, barsize = 5, fontsize = 8, label = "", color = 'blue') +
    labs(color = "Genus", fill = "Bootstrap Support") +
    theme(
      plot.margin = margin(-10, -10, -10, -10),
      # legend.title = element_text(face = "bold", size = 16),
      legend.position = 'bottom',
      legend.text = ggtext::element_markdown(size = 8),
      # legend.key.size = unit(2, "line")
    ) +
    # guides(colour = guide_legend(override.aes = list(size = 10, stroke = 10))) +
    NULL

# save_plot("ggtree/t7_newfamily.pdf", t9, base_height = 15, limitsize = FALSE)


# get annotations ---------------------------------------------------------

family_data <- bootstrap_data %>% 
  mutate(family = case_when(
    node %in% Descendants(as.phylo(t9), node = 4674)[[1]] ~ 'srw',
    node %in% Descendants(as.phylo(t9), node = 4129)[[1]] ~ 'srxa',
    node %in% Descendants(as.phylo(t9), node = 3601)[[1]] ~ 'srx',
    node %in% Descendants(as.phylo(t9), node = 4159)[[1]] ~ 'srbc',
    node %in% Descendants(as.phylo(t9), node = 4598)[[1]] ~ 'srsx',
    node %in% Descendants(as.phylo(t9), node = 4320)[[1]] ~ 'srsx',
    node %in% Descendants(as.phylo(t9), node = 4596)[[1]] ~ 'srsx',
    node %in% Descendants(as.phylo(t9), node = 4282)[[1]] ~ 'ppn+',
    node %in% Descendants(as.phylo(t9), node = 4064)[[1]] ~ 'srv',
    node %in% Descendants(as.phylo(t9), node = 4042)[[1]] ~ 'sru',
    node %in% Descendants(as.phylo(t9), node = 3941)[[1]] ~ 'srg',
    node %in% Descendants(as.phylo(t9), node = 3733)[[1]] ~ 'srt',
    node %in% Descendants(as.phylo(t9), node = 2749)[[1]] ~ 'srb',
    node %in% Descendants(as.phylo(t9), node = 2651)[[1]] ~ 'sre',
    node %in% Descendants(as.phylo(t9), node = 2522)[[1]] ~ 'sra',
    node %in% Descendants(as.phylo(t9), node = 2486)[[1]] ~ 'sra',
    node %in% Descendants(as.phylo(t9), node = 2567)[[1]] ~ 'srab',
    node %in% Descendants(as.phylo(t9), node = 3599)[[1]] ~ 'sro',
    node %in% Descendants(as.phylo(t9), node = 3560)[[1]] ~ 'srz',
    node %in% Descendants(as.phylo(t9), node = 3551)[[1]] ~ 'srr',
    node %in% Descendants(as.phylo(t9), node = 3452)[[1]] ~ 'srm',
    node %in% Descendants(as.phylo(t9), node = 3417)[[1]] ~ 'srj',
    node %in% Descendants(as.phylo(t9), node = 3245)[[1]] ~ 'srj',
    node %in% Descendants(as.phylo(t9), node = 3256)[[1]] ~ 'str',
    node %in% Descendants(as.phylo(t9), node = 3065)[[1]] ~ 'srd',
    node %in% Descendants(as.phylo(t9), node = 3060)[[1]] ~ 'srn',
    node %in% Descendants(as.phylo(t9), node = 3003)[[1]] ~ 'sri',
    node %in% Descendants(as.phylo(t9), node = 2778)[[1]] ~ 'srh',
    isTip == TRUE & branch_type != 'inner' ~ 'Unplaced'
  )) %>% 
  filter(!is.na(family))

write_rds(family_data, 'ggtree/family_data.rds')

chemo_families <- tibble(family = c("srh", "str", "sri", "srd", "srj", "srm", "srn",
                                        "sre", "sra", "srab", "srb",
                                        "srx", "srt", "srg", "sru", "srv", "srxa",
                                        "srw", "srz", "srbc", "srsx", "srr", "sro", 'ppn+'),
                             superfamily = c("Str", "Str", "Str", "Str", "Str", "Str", "Str",
                                             "Sra", "Sra", "Sra", "Sra",
                                             "Srg", "Srg", "Srg", "Srg", "Srg", "Srg",
                                             "Solo", "Solo", "Solo", "Solo", "Solo", "Solo", "Solo"))


full_species <- tribble(~species, ~full, ~n50,
                        'rculi', '*R. culicivorax*', 17583,
                        'bmala', '*B. malayi*', 14214749,
                        'mente', '*M. enterolobii*', 10529,
                        'mflor', '*M. floridensis*', 13230,
                        'mjava', '*M. javanica*', 14103,
                        'mgram', '*M. graminicola*', 20427,
                        'mhapl', '*M. hapla*', 37501,
                        'minco', '*M. incognita*', 38543,
                        'grost', '*G. rostochiensis*', 88495,
                        'gpall', '*G. pallida*', 120163,
                        'maren', '*M. arenaria*', 204341,
                        'hglyc', '*H. glycines*', 304127,
                        'bxylo', '*B. xylophilus*', 949830,
                        'cele', '*C. elegans*', 17493829,
                        'hcont', '*H. contortus*', 47382676)

family_plot <- family_data %>% 
    mutate(species = factor(species),
           family = factor(family)) %>% 
    group_by(clade, species, family) %>% 
    tally() %>% 
    left_join(chemo_families) %>% 
    left_join(full_species) %>% 
    filter(family != 'Unplaced') %>% 
    ggplot(aes(x = full, y = family)) +
    geom_tile(aes(fill = n)) +
    scale_fill_gradientn(colors = c('#FFFFFF', 
                                     '#40916c', '#2d6a4f',
                                    '#1b4332', '#081c15'),
                         # trans = 'log2'
                         ) +
    facet_grid(rows = vars(superfamily), 
               cols = vars(clade),
               scales = 'free', space = 'free') +
    labs(x = 'Species', y = 'Family', fill = 'Gene count') +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_markdown(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_markdown(size = 8),
      axis.ticks = element_line(size = 0.25, color = 'grey'),
      panel.border = element_rect(colour = "grey", size = 0.25, fill = NA),
      panel.grid = element_line(color = 'white'),
      legend.title = element_markdown(),
      legend.position = 'bottom'
    ) +
    NULL

# save_plot('ggtree/t10_heatmap.pdf',
#           plot_grid(t9, family_plot,
#                     rel_widths = c(1, 0.6)),
#           base_height = 5,
#           base_width = 10)


(species_bar <- family_data %>% 
    group_by(species, genus) %>% 
    summarise(n = n()) %>% 
    left_join(full_species) %>% 
    filter(!species %in% c('bmala', 'cele', 'hcont', 'rculi')) %>% 
    ggplot(aes(fct_reorder(full, -n), n)) +
    geom_col(aes(fill = genus)) +
    geom_text(aes(y = n + 5, 
                  color = genus,
                  label = prettyNum(n50, big.mark = ',')),
              fontface = 'italic', hjust = 0) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_color_manual(limits = family_data$genus,
                      values = family_data$color) +
    scale_fill_manual(limits = family_data$genus,
                      values = family_data$color) +
    labs(x = 'Species', y = 'Gene count', color = 'Genus', fill = 'Genus') +
    coord_flip() +
    theme_minimal() +
    theme(
      # axis.title.x = element_blank(),
      axis.text.y = element_markdown(size = 8),
      axis.ticks = element_line(size = 0.25, color = 'grey'),
      panel.border = element_rect(colour = "grey", size = 0.25, fill = NA),
      panel.grid = element_line(color = 'white'),
      legend.text = element_markdown(),
      legend.position = 'empty'
    ) +
    NULL)

# save_plot('ggtree/n50_bar.pdf', species_bar, base_width = 6)
