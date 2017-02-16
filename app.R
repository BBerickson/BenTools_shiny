
# setwd("~/Desktop/BenToolsTab/dockerTest")
# setwd("~/BenTools gh/BenTools_shiny")

# runApp("BenTools_shiny")
source("helper.R")

# server ----
server <- function(input, output) {
  
  # loads file
  first_file <- reactive({
    req(input$file$datapath)
    # add warnings for total size of LIST_DATA
    LIST_DATA <<- LoadTableFile(input$file$datapath, input$file$name, LIST_DATA)
    names(LIST_DATA$table_file)
  })
  
  # renders text
  output$txtout <- renderText({
    paste(input$txt, input$slider, format(input$date), sep = ", ")
  })
  
  # renders table
  output$table <- renderTable({
    first_file()
    head(LIST_DATA[[1]][[1]],4)
   
  })
  
  # sets fileType based on radio button
  observe({
    #print(input$checkGroup)
  })
  
  # renders check box
  output$fileNames <- renderUI({
    req(input$file$name)
    check_on <- c(sapply((LIST_DATA$gene_info), function(g) sapply(g, "[[",4)))
    checkboxGroupInput("checkGroup", label = h3("plot on/off"), 
                         choices = first_file(), selected = check_on) # build and use control list TODO
   
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
               fileInput("file", "File input:", accept = c('.table'), multiple = TRUE)
               
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