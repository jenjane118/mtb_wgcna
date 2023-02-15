library(shinyjs)
library(tidyverse)
library(RColorBrewer)
library(DT)
library(JBrowseR)

###############----------------**data**---------------##################

gene_info <- readRDS(here("R_data/complete_geneInfo.RData"))
gene_info <- gene_info %>% mutate(gene_ID = sub("gene:", "", gene_ID))
cds_mods <- read_tsv(here("Output/cds_funcats.tsv"), show_col_types = F)
cds_mods <- cds_mods %>% mutate(gene_ID = sub("gene:", "", gene_ID))
counts_condition <- readRDS(here("R_data/expression_cond_df.RData"))
counts_condition <- counts_condition %>% mutate(gene_ID = sub("gene:", "", gene_ID))
cond_order <- read_csv(here("Output/dend_labels.csv"), show_col_types = F)
cond_labels <- cond_order$condition
module_info <- read_tsv(here("Output/module_summary.tsv"), show_col_types = F)
mod_colors <- module_info$moduleColor
utr_df <- readRDS(here("R_data/utr_mods_df.Rdata"))
srna_df <- readRDS(here("R_data/srna_mods_df.Rdata"))
cds_df <- readRDS(here("R_data/cds_mods_df.RData"))

# source helper functions
source(here("shinyapp/helpers.R"))


# JBrowseR server for hosting local data files (have to have this running before run app)
data_server <- serve_data(here("shinyapp/sequences"))

