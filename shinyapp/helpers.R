
#####jBrowseR functions #######################

# create the necessary JB2 assembly configuration
assembly <- assembly(
  "http://127.0.0.1:5000/AL123456.3.fasta.gz",
  bgzip = TRUE
)

# create configuration for a JB2 GFF FeatureTrack
annotations_track <- track_feature(
  "http://127.0.0.1:5000/comb_filtered_ov_30_11.gff3.gz",
  assembly
)

# create the tracks array to pass to browser
tracks <- tracks(annotations_track)

# set up the default session for the browser
default_session <- default_session(
  assembly,
  c(annotations_track),
  display_assembly = F
)


#### functions ########

# find predicted name for putative ncRNA
ncrna_name <- function(ncrna){
  p_name <- gene_info %>%
    filter(gene_ID == ncrna) %>%
    select(pred_name) %>% pull()
  return(p_name)
}

#gene attribute/descr table for cds
gene_attrib <- function(cds){
  #need cds_funcats.tsv file (cds_mods)
  g_desc <- cds_mods %>%
    filter(tolower(gene_ID) == tolower(cds)) %>%
    select(Name, Functional_Category)
  return(g_desc)
}

#test gene locus input
test_input <- function(user_input){
  user_input <- tolower(user_input)
  if (user_input %in% tolower(gene_info$gene_ID)){
    if (grepl("rv", user_input)){
      type <- 'CDS'
    }else if (grepl("putative", user_input)){
      type <- 'ncrna'
    }else{
      type <- 'annotated'
    }
  }else{
    type <- 'invalid'
  }
  return(type)
}

# display selectable datatable of hubs
make_hubtable <- function(module){
  #mm_name <- paste("MM", module, sep="")
  cds_table <- cds_df %>% filter(moduleColor==module &
                                   MM > 0.70) %>%
    select(gene_ID, Name, start, end, MM) %>%
    mutate_if(is.numeric, ~round(., 2))
  srna_table <- srna_df %>% filter(mod_col==module &
                                     MM > 0.70) %>%
    select(pred_srna, srna_name, start, stop, MM) %>%
    mutate_if(is.numeric, ~round(., 2)) %>%
    dplyr::rename(Name=srna_name, gene_ID=pred_srna, end=stop)
  utr_table <- utr_df %>% filter(mod_col==module &
                                   MM > 0.70) %>%
    select(pred_utr, utr, start, stop, MM) %>%
    mutate_if(is.numeric, ~round(., 2)) %>%
    dplyr::rename(Name=utr, gene_ID=pred_utr, end=stop)
  hub_table <- rbind(cds_table, srna_table, utr_table)
  hub_table <- hub_table %>% arrange(desc(MM))
  return(hub_table)
}


#find UTRs adjacent to desired cds
find_utrs <- function(transcript){
  #requires utr_mods_df.RData
  utr_dt <- utr_df %>% 
    filter(tolower(nearest)==tolower(transcript)) %>%
    select(pred_utr, tss, utr, mod_col, MM) %>%
    mutate_if(is.numeric, ~round(., 2)) %>%
    dplyr::rename(pred_name=utr)
  return(utr_dt)
}


#find antisense opposite desired cds
find_as <- function(transcript){
  as_dt <- srna_df %>% 
    filter(tolower(ov_orf)==tolower(transcript)) %>%
    select(pred_srna, tss, srna_name, mod_col, MM) %>%
    mutate_if(is.numeric, ~round(., 2)) %>%
    dplyr::rename(pred_name=srna_name)
  return(as_dt)
}


#plot function for expression (needs 'counts_condition' dataframe)
# theme is masked by JBrowseR
txt_boxplot <- function(transcript){
  txt_df   <- counts_condition %>%
    filter(tolower(gene_ID)==tolower(transcript))
  txt_df$condition <- factor(txt_df$condition, 
                             levels = unique(cond_labels))
  ggplot(txt_df, aes(x=condition, y=counts)) +
    geom_boxplot(aes(fill=condition)) +
    coord_flip() +
    theme_bw() +
    # the x tick labels wouldn't show up in rshiny if at bottom
    ggplot2::theme(#axis.ticks.x = element_blank(),
      axis.text.y = element_text(face="bold"),
      panel.border = element_blank(),
      legend.position = "none") +
    ylab("Normalised counts")
}

# Return module color (eventually have text in appropriate color)
modColor <- function(transcript){
  #find moduleColor
  module_color <- gene_info %>% 
    filter(tolower(gene_ID)==tolower(transcript)) %>%
    select(moduleColor) %>% pull()
  return(module_color)
}

MM_value <- function(txt){
  #find module membership value
  modcol <- modColor(txt)
  if (modcol != ""){
    mm_name <- paste("MM", modcol, sep="")
    mmem <- round(gene_info %>%
                    filter(tolower(gene_ID)==tolower(txt)) %>%
                    select(all_of(mm_name)) %>% pull(), 2)
    return(mmem)
  }else{
    NULL
  }
}

distr_plot <- function(module){
  distr_df <- module_info %>% 
    filter(moduleColor==module) %>%
    select(n_utr, n_srna, n_ncrna, n_cds) %>% 
    pivot_longer(cols = c(n_utr, n_srna, n_ncrna, n_cds), names_to="type", values_to="number")
  ggplot(distr_df) +
    geom_col(aes(x=type, y=number, fill=type), show.legend = F) +
    scale_fill_brewer(palette = "YlGnBu") +
    theme_bw() +
    scale_x_discrete(labels = c("CDS", "annot-sRNA", "pred-sRNA", "pred-UTR"))
}


