

# setwd("~/Desktop/BenToolsTab/dockerTest")
# setwd("~/BenTools gh/BenTools_shiny")

# runApp("BenTools_shiny")
source("helper.R")

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 20MB.
options(shiny.maxRequestSize = 20*1024^2)

# server ----
server <- function(input, output) {
  # loads file
  first_file <- reactive({
    req(input$file$datapath)
    # add warnings for total size of LIST_DATA
    LIST_DATA <<-
      LoadTableFile(input$file$datapath, input$file$name, LIST_DATA)
    names(LIST_DATA$table_file)
  })
  
  # renders text
  output$txtout <- renderText({
    paste(input$txt, input$slider, format(input$date), sep = ", ")
  })
  
  # renders table
  output$table <- renderTable({
    first_file()
    head(LIST_DATA[[1]][[1]], 4)
    
  })
  
  # records check box on/off for common list
  observe({
    input$checkGroupCommon
    LIST_DATA$gene_info <<-
      CheckBoxOnOff("common", input$checkGroupCommon, LIST_DATA$gene_info)
  })
  
  # renders check box
  output$fileNamesCommon <- renderUI({
    req(input$file$name)
    check_on <-
      c(sapply((LIST_DATA$gene_info), function(g)
        sapply(g, "[[", 4)), "!NULL")
    checkboxGroupInput(
      "checkGroupCommon",
      label = h3("Plot on/off"),
      choices = first_file(),
      selected = check_on
    ) # build and use control list TODO
  })
  
  observe({
    input$plotBinRange
    kplotBinRange <<- c(kplotBinRange[1:2], input$plotBinRange)
  })
  
  # Specification of range within an interval to plot
  output$rangeBin <- renderUI({ 
    req(input$file$name) # way of triggering only once ? TODO
    sliderInput("plotBinRange", label = h3("Plot Range:"),
                min = kplotBinRange[1], max = kplotBinRange[2], value = kplotBinRange[3:4])
  })
  
}

# UI ----
ui <- tagList(
  navbarPage(
    theme = shinythemes::shinytheme("cerulean"),
    "BenTools",
    tabPanel(
      "Navbar 1",
      sidebarPanel(tabsetPanel(tabPanel(
        "Load file(s)",
        fileInput(
          "file",
          "File input:",
          accept = c('.table'),
          multiple = TRUE
        )
        
      ), tabPanel(
        "Plot Options",
        radioButtons("myMath", 
                     "Set math function", 
                     choices = c("mean", "sum", "median", "var"), 
                     selected = "mean")
        
      )
      ),
      tabsetPanel(tabPanel(
        "Common Gene list",
        uiOutput("fileNamesCommon"),
        uiOutput("rangeBin")
        
      ))
      ),
      mainPanel(tabsetPanel(
        tabPanel(
          "Tab 1",
          h4("Table"),
          tableOutput("table"),
          h4("Verbatim text output"),
          verbatimTextOutput("txtout"),
          h1("Header 1"),
          h2("Header 2"),
          h3("Header 3"),
          h4("Header 4"),
          h5("Header 5")
        ),
        tabPanel("Tab 2"),
        tabPanel("Tab 3")
      ))
    ),
    tabPanel("Navbar 2"),
    tabPanel("Navbar 3")
  )
)

# exicute ----
shinyApp(ui = ui, server = server)