# RNA Editing Shiny App for plotting and subsetting

source("src/VariantPlottingFunctions.R")
source("src/MAFPrepFunctions.R")

# Define UI for data upload app ----
ui <- fluidPage(
  title = "Variant Analysis Plots",
  sidebarLayout(
    sidebarPanel(
      # Input: Select a file ----
      fileInput("mafFile", "Upload MAF File",
                multiple = FALSE,
                accept = c(".maf")),
      hr(),
      # conditionalPanel(
      #   'input.dataset === "MAF"',
      #   checkboxGroupInput("cond", "Condition", levels(maf$Condition), selected=levels(maf$Condition))),
      selectInput("plotType", "Plot Type",
                  c("Violin" = "VAFViolin",
                    "Box Plot" = "VAFBoxPlot",
                    "Top Genes" = "TopGenes",
                    "Variants Per Gene" = "VariantsPerGene",
                    "Mutated Genes Per Sample - Boxplot" = "MutatedGenesPerSampleBox")),
      hr(),
      actionButton(inputId = "plot_now", label = "Plot"),
      downloadButton('downloadPlot',"Save Plot"),
      hr()
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", id="plottab", plotOutput("Plot"))
        ),
      hr(),
      tabsetPanel(
        id = 'dataset',
        tabPanel("MAF", id="maftab", DT::dataTableOutput("maf"))
      )
    )
  )
)



# Define server logic to read selected file ----
server <- function(input, output) {
  
  uploadMaf <- eventReactive(input$mafFile$datapath, {
    # input$mafFile will be NULL initially. After the user selects and uploads a file
    req(input$mafFile)
    tryCatch({
      maf <- fread(input$mafFile$datapath,
                   stringsAsFactors = FALSE,
                   sep = "\t")},
      error = function(e) {stop(safeError(e))})
      maf$Condition <- factor(maf$Tumor_Sample_Barcode, levels=sort(unique(maf$Tumor_Sample_Barcode)))
      maf
  })
  
  
  output$maf <- DT::renderDataTable({
    maf <- uploadMaf()  
    DT::datatable(maf)
    # DT::datatable(subset(maf, Condition %in% input$cond))
  })
  
    
  myPlot <- function() {
    maf <- uploadMaf()
    
    if(input$plotType == "VAFViolin"){
      ggViolin(maf)
    }
    else if(input$plotType == "VAFBoxPlot") {
      ggBox(maf)
    }
    else if(input$plotType == "TopGenes") {
      topGenes(maf)
    }
    else if(input$plotType == "VariantsPerGene") {
      VariantsPerGene(maf)
    }
    else if(input$plotType == "MutatedGenesPerSampleBox") {
      MutatedGenesPerSampleBox(maf)
    }
  }
  output$Plot <-renderPlot({
    myPlot()
  })
  
  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste0(input$dataset, "_", input$plotType, "_", gsub("-", "", Sys.Date()), sep=".pdf")
    },
    content = function (file) {
      pdf(file)
      print(myPlot())
      dev.off()
    }
  )
}

# Create Shiny app ----
shinyApp(ui, server)