
# setwd("~/Desktop/BenToolsTab/dockerTest")

# runApp("BenTools_shiny")
source("helper.R")

# run load needed pakages using my_pakages(x) ----
suppressPackageStartupMessages(my_packages(
  c(
    "shiny",
    "shinythemes",
    "ggplot2",
    "tcltk2",
    "dplyr",
    "tidyr",
    "readr",
    "fastcluster",
    "RColorBrewer"
  )
))

# server ----
server <- function(input, output) {
  first_file <- reactive({
    req(input$file$datapath)
    read.table(input$file$datapath, header = TRUE, stringsAsFactors= FALSE, comment.char = "")
  })
  
  output$txtout <- renderText({
    paste(input$txt, input$slider, format(input$date), sep = ", ")
  })
  output$table <- renderTable({
    if(is.null(first_file())){
      head(cars, 4)
    } else{
      head(first_file(),4)
    }
    
    
  })
  observe({
    if(is.null(first_file())){
      print(head(first_file()))
    }
  })
  
  output$fileNames <- renderUI({
    req(input$file$name)
    
      checkboxGroupInput("checkGroup", label = h3("plot on/off"), 
                         choices = input$file$name,selected=input$file$name)
   
  })
  
}

# UI ----
ui <- tagList(
  navbarPage(
    theme = shinythemes::shinytheme("cerulean"),
    "BenTools",
    tabPanel("Navbar 1",
             sidebarPanel(
               tabsetPanel(
                 tabPanel("load file",
               fileInput("file", "File input:", accept = c('.table')),
               radioButtons('fileType', 'File type',
                            c('Table' = '.table',
                              'Gene list' = '.txt',
                              'color file' = '.color.txt'),
                            '.table')
             )),
             tabsetPanel(
               tabPanel("Common Gene list",
                        uiOutput("fileNames")

               ))
             ),
             mainPanel(
               tabsetPanel(
                 tabPanel("Tab 1",
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
               )
             )
    ),
    tabPanel("Navbar 2"),
    tabPanel("Navbar 3")
  )
)

# exicute ----
shinyApp(ui = ui, server = server)