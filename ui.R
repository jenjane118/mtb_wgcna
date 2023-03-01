#define UI for MTB modules app

ui <- fluidPage(
  useShinyjs(),
  # title
  titlePanel("MTB WGCNA Modules Explorer"),
  
  #sidebar layout with input and output definitions
  sidebarLayout(
    #sidebar panel
    sidebarPanel(
      # input: select to explore by module or transcript name
      selectInput("explore", "Explore by:",
                  c("module", "transcript", "coordinates"), selected = "module"),
      
      #only if module is chosen see these choices
      conditionalPanel(
        condition = "input.explore == 'module' ",
        selectInput(
          "module_color", "Module",
          choices = mod_colors, selected = "green"),
        # distribution of module components
        h4(textOutput("plotCaption1")),
        #plot distribution of module transcripts
        plotOutput("distrPlot")
        ),  #conditionalPanel explore=module
      
      conditionalPanel(
        condition = "input.explore == 'transcript' ",
          #input: enter text of transcript name
          textInput("transcript_locus", "Transcript Name"),
          #action button to delay until button entered
          actionButton("go", "GO"),
          h4(textOutput("caption1")),
          h4(htmlOutput("moduleAssign")),
          br(),
          h4(textOutput("caption2")),
          h4(textOutput("modMembership")),
          br(),
          h5(tableOutput("geneDescr")),
          br(),
          h5(textOutput("ncrnaCaption")),
          h4(textOutput("predName"))
      ), #conditionalPanel transcript
    
      conditionalPanel(
        condition = "input.explore == 'coordinates' ",
        numericInput("start", "Start (rel to + strand)", 
                     value=1, min=1, max=4411532),
        numericInput("end", "Stop", value=10000, min=1, max=4411533)
      ) #condPanel coordinates
      
      
    ), #sidebarPanel
    
    
    #main panel to display output
    mainPanel(
      
      #output: if exploring by module
      conditionalPanel(
        condition = "input.explore == 'module' ",
        #filter by MM
        sliderInput("mmSlider", "MM filter", 
                    min=0, max=0.9, value=0.7, step=0.1, round=-1),
        #module output: list of hubs
        h4(textOutput("tableCaption")),
        dataTableOutput("hubList"),
        br(),
        div(id="dwnbutton", 
            downloadButton("downloadHubs", "Download")
        )
      ),  #conditionalPanelmodule
      
      #output: if exploring by transcript
      conditionalPanel(
        condition = "input.explore == 'transcript' ",
        h4(textOutput("expr_caption")),
        plotOutput("exprPlot"),
        JBrowseROutput("browserOutput_transcript"),
        h4(textOutput("caption_utrs")),
        tableOutput("findUTRs"),
        div(id="dwnbutton", 
            downloadButton("downloadUTRs", "Download UTRs")
        ),
        h4(textOutput("caption_antisense")),
        tableOutput("findAntisense"),
        div(id="dwnbutton", 
            downloadButton("downloadAntisense", "Download Antisense")
        )
      ),  #conditionalPanel transcript
      
      #output if exploring by coordinates
      conditionalPanel(
        condition = "input.explore == 'coordinates' ",
        JBrowseROutput("browserOutput"),
        h3(textOutput("predicted_transcript_caption")),
        tableOutput("srna_table"),
        div(id="dwnbutton", 
            downloadButton("downloadSRNAs", "Download sRNAs")
        ),
        tableOutput("utr_table"),
        div(id="dwnbutton", 
            downloadButton("downloadUTRs_coord", "Download UTRs")
        )
      )#condPanel coordinates
      
    )#mainPanel
    
  )#sidebarLayout
  
) #fluidPage

