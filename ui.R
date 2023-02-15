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
          #radioButtons("visuBtn", NULL, 
          #          choices = c("Expression Plots" = "plot",
          #                      "Associated non-coding RNA" = "table")),
          h5(textOutput("ncrnaCaption")),
          h4(textOutput("predName"))
      ), #conditionalPanel transcript
    
      conditionalPanel(
        condition = "input.explore == 'coordinates' ",
        selectInput("strand", "Strand", 
                    choices = c("+","-"),
                    selected = "+"),
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
        #module output: list of hubs
        h4(textOutput("tableCaption")),
        dataTableOutput("hubList"),
        JBrowseROutput("browserOutput_module")
       
      ),  #conditionalPanelmodule
      #output: if exploring by transcript
      conditionalPanel(
        condition = "input.explore == 'transcript' ",
          h4(textOutput("expr_caption")),
          plotOutput("exprPlot"),
          h4(textOutput("caption_utrs")),
          tableOutput("findUTRs"),
          h4(textOutput("caption_antisense")),
          tableOutput("findAntisense")
      ),  #conditionalPanel transcript
      
      #output if exploring by coordinates
      conditionalPanel(
        condition = "input.explore == 'coordinates' ",
        h3(textOutput("predicted_transcript_caption")),
        tableOutput("srna_table"),
        tableOutput("utr_table"),
        JBrowseROutput("browserOutput")
      )#condPanel coordinates
      
    )#mainPanel
    
  )#sidebarLayout
  
) #fluidPage

