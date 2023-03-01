# server logic for shinyapp 'Mtb WGCNA module explorer'
# author: Jennifer J Stiens
# email: j.j.stiens@gmail.com

server <- function(input, output) {
  
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
    paste("Members of ", input$module_color, "with MM > ", input$mmSlider)
  })
  # render hubtable caption
  output$tableCaption <- renderText({
    tableCaption()
  })
  make_hublist <- reactive({
    req(input$module_color, input$mmSlider)
    my_dt <- make_hubtable(input$module_color, input$mmSlider)
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
  
  # set up event listener on table selection (not currently used)
  table_selection <- reactive({
    hubTable <- make_hubtable(input$module_color)
    hubTable[input$hubList_rows_selected, ]
  })
  
  #download all 'hubs'
  output$downloadHubs <- downloadHandler(
    filename = function(){
      paste(input$module_color, "_hubs", '.csv', sep='')
    },
    content = function(file) {
      write.csv(make_hublist(), file)
    }
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
    ncrna <- c("sRNA", "UTR")
    if (type %in% ncrna) {
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
    ncrna <- c("sRNA", "UTR")
    if (type %in% ncrna) {
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
  
  #download button for utr table
  observeEvent(input$transcript_locus, {
    toggleState("downloadUTRs",
                condition=input$transcript_locus
    )
  })
  #download table of adjacent UTRs
  output$downloadUTRs <- downloadHandler(
    filename = function(){
      paste('UTRS_', input$transcript_locus, '.csv', sep='')
    },
    content = function(file) {
      write.csv(list_utrs(), file)
    }
  ) 
  
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

  #download button for antisense
  observeEvent(input$transcript_locus, {
    toggleState("downloadAntisense",
                condition=input$transcript_locus
    )
  })
  #download table of antisense
  output$downloadAntisense <- downloadHandler(
    filename = function(){
      paste('antisense_', input$transcript_locus, '.csv', sep='')
    },
    content = function(file) {
      write.csv(list_as(), file)
    }
  )
  
  # JBrowseR window for transcript input
  location_transcript <- reactive({
    req(input$transcript_locus)
    type = test_input(input$transcript_locus)
    if (type=='CDS'){
      buf_start <- cds_df %>% filter(gene_ID==input$transcript_locus) %>%
        select(start) %>% pull() - 200
      buf_stop <- cds_df %>% filter(gene_ID==input$transcript_locus) %>%
        select(end) %>% pull() + 200
    }else if (type=='sRNA'){
      buf_start <- srna_df %>% filter(pred_srna==input$transcript_locus) %>%
        select(start) %>% pull() - 200
      buf_stop <- srna_df %>% filter(pred_srna==input$transcript_locus) %>%
        select(stop) %>% pull() + 200
    }else if (type=='UTR'){
      buf_start <- utr_df %>% filter(pred_utr==input$transcript_locus) %>%
        select(start) %>% pull() - 200
      buf_stop <- utr_df %>% filter(pred_utr==input$transcript_locus) %>%
        select(stop) %>% pull() + 200
    }
    str_glue(
      "AL123456.3:",
      "{buf_start}",
      "..",
      "{buf_stop}"
    )
  })
  # activate download button only if location function not null
  observeEvent(list_as(), {
    toggleState("browserOutput_transcript",
                condition= !is.null(location_transcript())
                
    )
  })
  # render JBrowseR widget
  output$browserOutput_transcript <- renderJBrowseR(
    JBrowseR(
      "View",
      assembly = assembly,
      tracks = tracks,
      location = location_transcript(),
      defaultSession = default_session
    )
  )
  
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
      filter(start >= input$start & stop <= input$end) %>%
      select(pred_srna, srna_name, mod_col, MM, ov_orf )
  })
  
  get_utrs <- reactive({
    req(input$start, input$end)
    validate(
      need(input$start < input$end, "Start coordinate must be less than end coordinate.")
    )
    utr_df %>%
      filter(start >= input$start & stop <= input$end) %>%
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
      buf_start <- input$start - 200
      buf_stop <- input$end + 200
      str_glue(
        "AL123456.3:",
        "{buf_start}",
        "..",
        "{buf_stop}"
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

  #if there are sRNAs...
  observeEvent(get_srna(), {
    toggleState("downloadSRNAs",
                condition=nrow(get_srna()) > 0
    )
  })
  
  #download table of sRNAs
  output$downloadSRNAs <- downloadHandler(
    filename = function(){
      paste('sRNAs_', input$start, "_", input$end, '.csv', sep='')
    },
    content = function(file) {
      write.csv(get_srna(), file)
    }
  )
  
  #if there are UTRs...
  observeEvent(get_utrs(), {
    toggleState("downloadUTRs_coord",
                condition=nrow(get_utrs()) > 0
    )
  })
  
  #download table of UTRs
  output$downloadUTRs_coord <- downloadHandler(
    filename = function(){
      paste('UTRs_', input$start, "_", input$end, '.csv', sep='')
    },
    content = function(file) {
      write.csv(get_utrs(), file)
    }
  )
  
} #server  
  
  
  
