library(tidyverse)
library(here)
library(conflicted)
library(biomaRt)

conflict_prefer("filter", "dplyr")


# Get orthologes of C. elegans AmineRs ------------------------------------
seeds <- read_csv(here('AmineR', 'aux', 'celegans_aminer.txt'), col_names = 'seed')

# establish the BioMart that will be used
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# set parameters
filters = c('wbps_gene_id')
value = list(seeds$seed)
attributes = c('wbps_peptide_id',
               'acviteprjeb1697_homolog_ensembl_peptide',
               'anceylprjna72583_homolog_ensembl_peptide',
               'assuumprjna62057_homolog_ensembl_peptide',
               'brmalaprjna10729_homolog_ensembl_peptide',
               'brpahaprjeb497_homolog_ensembl_peptide',
               'brtimoprjeb4663_homolog_ensembl_peptide',
               'buxyloprjea64437_homolog_ensembl_peptide',
               'caelegprjna13758_homolog_ensembl_peptide',
               'diimmiprjeb1797_homolog_ensembl_peptide',
               'drmediprjeb500_homolog_ensembl_peptide',
               'elelapprjeb502_homolog_ensembl_peptide',
               'hacontprjeb506_homolog_ensembl_peptide',
               'hepolyprjeb15396_homolog_ensembl_peptide',
               'lisigmprjeb3075_homolog_ensembl_peptide',
               'loloaprjna246086_homolog_ensembl_peptide',
               'neamerprjna72135_homolog_ensembl_peptide',
               'onflexprjna230512_homolog_ensembl_peptide',
               'onocheprjeb1204_homolog_ensembl_peptide',
               'onvolvprjeb513_homolog_ensembl_peptide',
               'parediprjna186477_homolog_ensembl_peptide',
               'paunivprjna386823_homolog_ensembl_peptide',
               'patricprjeb515_homolog_ensembl_peptide',
               'prpaciprjna12644_homolog_ensembl_peptide',
               'stcarpprjna202318_homolog_ensembl_peptide',
               'strattprjeb125_homolog_ensembl_peptide',
               'tocaniprjna248777_homolog_ensembl_peptide',
               'trspirprjna12603_homolog_ensembl_peptide',
               'trmuriprjeb126_homolog_ensembl_peptide',
               'wubancprjna275548_homolog_ensembl_peptide'
               )

# pull the desired data
orthologs <- getBM(mart = mart, 
                  filters = filters,
                  value = value,
                  attributes = attributes) %>%
  janitor::clean_names() %>%
  pivot_longer(cols = 1:30) %>%
  dplyr::select(name, transcript_id = value) %>%
  filter(transcript_id != '') %>%
  distinct() %>%
  arrange(transcript_id)

# use the orthologs of Ce AmineRs as seeds for a new query
orthologs2 <- getBM(mart = mart, 
                    filters = 'wbps_peptide_id',
                    value = orthologs$transcript_id,
                    attributes = attributes) %>%
  janitor::clean_names() %>%
  pivot_longer(cols = 1:30) %>%
  dplyr::select(name, transcript_id = value) %>%
  filter(transcript_id != '') %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  arrange(transcript_id)

# get the paralogs 
paralogs <- getBM(mart = mart, 
                  filters = 'wbps_peptide_id',
                  value = orthologs2$transcript_id,
                  attributes = c('species_id_key', 'wbps_paralog_paralog_peptide_id')) %>%
  janitor::clean_names() %>%
  dplyr::select(name = species_id_key, transcript_id = 2) %>%
  filter(transcript_id != '') %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  arrange(transcript_id)

final_list <- bind_rows(orthologs, orthologs2, paralogs) %>%
  arrange(transcript_id) %>%
  distinct(transcript_id, .keep_all = TRUE) %>%
  mutate(species = case_when(
    name == 'wbps_peptide_id' ~ 'caelegprjna13758_homolog_ensembl_peptide',
    TRUE ~ name)) %>%
  mutate(species = str_remove(species, "_homolog_ensembl_peptide")) %>%
  mutate(
    species = str_sub(species, 1, 6),
    species = case_when(
      species == 'loloap' ~ 'loloa',
      TRUE ~ species
    )) %>%
  mutate(transcript_id = case_when(
    species %in% c('brmala', 'caeleg', 'onvolv', 'prpaci', 'trmuri') ~ str_remove(transcript_id, '\\.[0-9]$'),
    TRUE ~ transcript_id
  )) %>%
  left_join(., read_delim(here('AmineR', 'aux', 'species_liftover.txt'), delim = '\t', col_names = c('species', 'full', 'BioProject'))) %>%
  distinct(transcript_id, .keep_all = TRUE) 


out <- dplyr::select(final_list, transcript_id, full) %>% 
  group_by(full) %>%
  group_nest()

# write out transcript_ids by species
out %>% pwalk(~write_csv(x = .y, path = str_glue("~/Box/ZamanianLab/Data/Genomics/AmineR/identification_pipeline/3/", .x, "_1_ids.txt")))

write_csv(dplyr::select(final_list, transcript_id, path),
          "~/Box/ZamanianLab/Data/Genomics/AmineR/identification_pipeline/3/all_amineR.txt",
          col_names = FALSE)
