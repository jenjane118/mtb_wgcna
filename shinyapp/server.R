
#define server logic to retrieve transcript info, plot expression

server <- function(input, output) {
  
  library(shiny)
  library(shinyjs)
  library(dplyr)
  library(ggplot2)
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
  
  #### select module #############################
  ### input$module_color
  
  # caption of type distr plot
  plotCaption1 <- reactive({
    paste("Distribution of", input$module_color, "transcript types")
  })
  #render type distr plot caption
  output$plotCaption1 <- renderText({
    plotCaption1()
  })
  # create reactive plot
  make_distPlot <- reactive({
    req(input$module_color)
    distr_plot(input$module_color)
  })
  #render distribution type plot
  output$distrPlot <- renderPlot({
    make_distPlot()
  })
  
  # caption for hub table
  tableCaption <- reactive({
    paste("Hub members (MM>0.70) of ", input$module_color)
  })
  # render hubtable caption
  output$tableCaption <- renderText({
    tableCaption()
  })
  make_hublist <- reactive({
    req(input$module_color)
    my_dt <- make_hubtable(input$module_color)
    my_dt %>% select(gene_ID, Name, MM)
  })
  
  # render selectable hub table
  output$hubList <- renderDataTable({
    make_hublist()
    }, selection = "single", 
    options = list(columnDefs = list(list(className = "dt-center", 
                           targets = "_all"))
    )
  )
  
  # set up event listener on table selection
  table_selection <- reactive({
    hubTable <- make_hubtable(input$module_color)
    hubTable[input$hubList_rows_selected, ]
  })
  
  #uses data table input
  location_module <- reactive({
    print(table_selection())
    if (is.na(table_selection()$start[1])) {
      "AL123456.3:1..500"
    } else {
      str_glue(
        "AL123456.3:",
        "{table_selection()$start}",
        "..",
        "{table_selection()$end}"
      )
    }
  })
  
  # link the UI with the browser widget
  output$browserOutput_module <- renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly,
      tracks = tracks,
      location = location_module(),
      defaultSession = default_session
    )
  )
  
  ###### select transcript ########################
  #input$transcript_locus
  
  # reactive caption for expression plot
  plotCaption2 <- eventReactive( input$go, {
    paste("Expression over conditions:", input$transcript_locus)
  })
  #render caption for expression plot
  output$expr_caption <- renderText({
    plotCaption2()
  })
  
  make_exprPlot <- eventReactive(input$go, {
    req(input$transcript_locus)
    txt_boxplot(input$transcript_locus)
  })
  #render expression plot
  output$exprPlot <- renderPlot({
    make_exprPlot()
  })
  
  #reactive caption for module membership value
  mm_caption <- eventReactive( input$go, {
    my_col <- modColor(input$transcript_locus)
    paste("Module Membership in", my_col, "Module:")
  })
  #render caption for mod mem value
  output$caption2 <- renderText({
    mm_caption()
  })
  #reactive output for mm value
  find_mmValue <- eventReactive( input$go, {
    req(input$transcript_locus)
    validate(
      need(tolower(input$transcript_locus) %in% tolower(gene_info$gene_ID), "Invalid transcript")
    )
    MM_value(input$transcript_locus)
  })
  #render text of module membership value
  output$modMembership <- renderText({
    find_mmValue()
  })
  
  #caption for module assignment
  mod_assign_text <- eventReactive( input$go, {
    paste("Module assignment:", input$transcript_locus)
  })
  output$caption1 <- renderText({
    mod_assign_text()
  }) 
  
  #reactive module assignment
  find_module <- eventReactive( input$go, {
    req(input$transcript_locus)
    validate(
      need(test_input(input$transcript_locus) != "invalid", 
           "Transcript is not found in data. Please use H37Rv locus number or valid putative ncRNA ID.")
    )
    my_col <- modColor(input$transcript_locus)
    paste('<span style=\"background-color:', my_col,
          '\">', my_col, '<span>', sep="")
  })
  # render text of module assignment
  output$moduleAssign <- renderText({
    find_module()
  })
  
  #create gene descr
  make_geneDescr <- eventReactive(input$go, {
    req(input$transcript_locus)
    type <- test_input(input$transcript_locus)
    if (type=='CDS') {
      gene_attrib(input$transcript_locus)
    }else{
      NULL
    }
  })
  
  #render gene description
  output$geneDescr <- renderTable({
    make_geneDescr()
  })
  
  # caption for pred name 
  re_ncrnaCaption <- eventReactive( input$go, {
    req(input$transcript_locus)
    type <- test_input(input$transcript_locus)
    if (type=='ncrna') {
      print("Predicted locus ID:")
    }else{
      NULL
    }
  })
  #render caption for predicted name 
  output$ncrnaCaption <- renderText({
    re_ncrnaCaption()
  })
  
  # reactive text for predname
  find_predname <- eventReactive( input$go, {
    req(input$transcript_locus)
    type <- test_input(input$transcript_locus)
    if (type=='ncrna') {
      ncrna_name(input$transcript_locus)
    }else{
      NULL
    }
  })
  # render predicted name text
  output$predName <- renderText({
    find_predname()
  })
  
  #reactive caption
  utr_caption <- eventReactive( input$go, {
    req(input$transcript_locus)
    type <- test_input(input$transcript_locus)
    if (type=='CDS') {
      print("Adjacent UTRs:")
    }else{
      NULL
    }
  })
  #render caption for adjacent utrs 
  output$caption_utrs <- renderText({
    utr_caption()
  })
  
  # #reactive table of adj utrs
  list_utrs <- eventReactive( input$go, {
    req(input$transcript_locus)
    type <- test_input(input$transcript_locus)
    if (type=='CDS') {
      find_utrs(input$transcript_locus)
    }else{
      NULL
    }
  })
  # # render table of adjacent utrs
  output$findUTRs <- renderTable({
    list_utrs()
  })
  
  # make antisense table caption
  caption_as <- eventReactive( input$go, {
    type <- test_input(input$transcript_locus)
    if (type=='CDS') {
      print("Overlapping antisense RNAs:")
    }else{
      NULL
    }
  })
  #render caption antisense
  output$caption_antisense <- renderText({
    caption_as()
  })

  #reactive table of adj antisense
  list_as <- eventReactive( input$go, {
    #requires srna_mods_df.RData
    req(input$transcript_locus)
    type <- test_input(input$transcript_locus)
    if (type=='CDS') {
      find_as(input$transcript_locus)
    }else{
      NULL
    }
  })
  # render table of antisense
  output$findAntisense <- renderTable({
    list_as()
  })

  
########## select coordinates ####################  
  ## input$start, input$end

  
  output$predicted_transcript_caption <- renderText({
    req(input$start, input$end)
      print("Predicted ncRNA transcripts within coordinates:")
  })
  
  #function to return all transcripts within coordinates
  get_srna <- reactive({
    req(input$start, input$end)
    validate(
      need(input$start < input$end, "Start coordinate must be less than end coordinate."),
      )
    srna_df %>% 
      filter(strand == input$strand & start >= input$start & stop <= input$end) %>%
      select(pred_srna, srna_name, mod_col, MM, ov_orf )
  })
  
  get_utrs <- reactive({
    req(input$start, input$end)
    validate(
      need(input$start < input$end, "Start coordinate must be less than end coordinate.")
    )
    utr_df %>%
      filter(strand == input$strand & start >= input$start & stop <= input$end) %>%
      select(pred_utr, utr, mod_col, MM, tss)
  })
  
  output$srna_table <- renderTable({
      get_srna()
  })
  
  output$utr_table <- renderTable({
      get_utrs()
  })
  
  #JBrowseR for coordinates page
  
  location <- reactive({
    if (input$start == 1) {
      "AL123456.3:1..500"
    } else {
      str_glue(
        "AL123456.3:",
        "{input$start}",
        "..",
        "{input$end}"
      )
    }
  })
  
  # link the UI with the browser widget (coordinates page)
  output$browserOutput <- renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly,
      tracks = tracks,
      location = location(),
      #location = "AL123456.3:1..2500" ,
      defaultSession = default_session
    )
  )
  
  
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
  
  
  
} #server  
  
  
  
