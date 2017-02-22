


# setwd("~/Desktop/BenToolsTab/dockerTest")
# setwd("~/BenTools gh/BenTools_shiny")

# runApp("BenTools_shiny")
source("helper.R")

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 20MB.
options(shiny.maxRequestSize = 20 * 1024 ^ 2)

# server ----
server <- function(input, output) {
  # load file, save data and set up info ----
  LIST_DATA <- reactiveValues(
    table_file = list(),
    # [[]] gene X1 X2 ...
    gene_file = list(),
    # holds $common genes from files and $gene file(s)
    gene_info = list(),
    # for holding gene file info in a list of lists, a set for $common and each $gene file(s) [c("dot", "line", "color", plot?, NickName,nrom)]
    clust = list()
    # Cluster holder
  )
  # loads file(s) 
  first_file <- reactive({
    req(input$file$datapath)
    # add warnings for total size of LIST_DATA
    LIST_DATA <-
      LoadTableFile(input$file$datapath, input$file$name, LIST_DATA)
    names(LIST_DATA$table_file)
  })
  
  # records check box on/off for common list and builds data frame and info for plot
  observe({
    req(first_file())
    LIST_DATA$gene_info <-
      CheckBoxOnOff("common", input$checkGroupCommon, LIST_DATA$gene_info)
    if (!is.null(LIST_DATA$table_file)) {
      Make_Data_Frame <<- MakeDataFrame(LIST_DATA) # set out
      if (!is.null(Make_Data_Frame[[1]])) {
        Apply_Math <<- ApplyMath(Make_Data_Frame, input$myMath)
      }
    }
  })
  
  

  observe({
    input$myMath
    if (!is.null(Make_Data_Frame[[1]])) {
      Apply_Math <<- ApplyMath(Make_Data_Frame, input$myMath)
      print(GGplotLineDot(Apply_Math)) # renderplot
    }
    
  })
  # set up and control check box actions common ----
  
  # renders check box when file is loaded
  output$fileNamesCommon <- renderUI({
    print("hi")
    checkboxGroupInput(
      "checkGroupCommon",
      label = h3("Plot on/off"),
      choices = first_file(),
      selected = unique(c(
        sapply(isolate(LIST_DATA$gene_info$common), "[[", 4),
        last(first_file())))
    )
  })
  
  # plots when acction button is pressed
  observe({
    req(input$myplot)
    if (!is.null(Apply_Math)) {
      Y_Axis_Lable <<- YAxisLable()
      Lines_Lables_List <<- LinesLablesList()
      print(GGplotLineDot(Apply_Math)) # renderplot
    }
  })
  
  # set up and control bin range action ----
  
  # Specification of range within an interval to plot
  output$rangeBin <- renderUI({
    req(input$file$name)
    sliderInput( #  try and only render 1 time? how to update min max and values?
      "plotBinRange",
      label = h3("Plot Range:"),
      min = kplotBinRange[1],
      max = kplotBinRange[2],
      value = kplotBinRange[3:4]
    )
  })
  # plots when range bin slider is triggered
  observe({
    req(input$plotBinRange)
    kplotBinRange <<- c(kplotBinRange[1:2], input$plotBinRange)
    if (!is.null(Make_Data_Frame[[1]])) {
      print(GGplotLineDot(Apply_Math)) # renderplot
    }
  })
  
}

# UI -----
ui <- tagList(
  navbarPage(
    theme = shinythemes::shinytheme("cerulean"),
    "BenTools",
    tabPanel(
      "Navbar 1",
      sidebarPanel(
        tabsetPanel(
          
        # load file(s) tab ----
        tabPanel(
          "Load file(s)",
          fileInput(
            "file",
            label = "File input:",
            accept = c('.table'),
            multiple = TRUE
          ),
          uiOutput("rangeBin")
        ),
        
        # plot options tab ----
        tabPanel(
          "Plot Options",
          radioButtons(
            "myMath",
            "Set math function",
            choices = c("mean", "sum", "median", "var"),
            selected = "mean"
          )
          
        )
      ),
      
      #common genes tab ----
      tabsetPanel(
        tabPanel(
          "Common Gene list",
          uiOutput("fileNamesCommon"),
          actionButton("myplot", "Update Plot")
          
        )
      )),
      # main panel for plot ----
      mainPanel(tabsetPanel(
        tabPanel("Tab 1",
                 h4("Table Plot"),
                 plotOutput("plot")),
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