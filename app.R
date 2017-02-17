

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
  
  # records check box on/off for common list and builds data frame and info for plot
  observe({
    input$checkGroupCommon
    LIST_DATA$gene_info <<-
      CheckBoxOnOff("common", input$checkGroupCommon, LIST_DATA$gene_info)
    
    print("hi")
    
  })
  
  make_data_frame <- reactive({
    tt <- first_file()  #MakeDataFrame(LIST_DATA)
    print((tt))
    c(tt,tt)
  })
  
  
  # makes applied math data frame
  apply_math <- reactive({
    req(input$file$name)
    makedataframe <- make_data_frame()
    print(names(makedataframe))
    # ApplyMath(
    #   makedataframe$list_data_frame,
    #   makedataframe$use_col,
    #   makedataframe$use_dot,
    #   makedataframe$use_line,
    #   makedataframe$use_size,
    #   makedataframe$use_x_label,
    #   makedataframe$legend_space
    #   )
  })
  
  # renders plot
  observe({
    req(input$file$name)
    # applymath <- apply_math()
    # print(names(applymath))
    # output$plot <- renderPlot({
    #   GGplotF(
    #     applymath$list_long_data_frame,
    #     applymath$use_col,
    #     applymath$use_dot,
    #     applymath$use_line,
    #     applymath$use_size,
    #     applymath$use_y_label,
    #     applymath$use_x_label,
    #     applymath$use_plot_breaks,
    #     applymath$virtical_line_data_frame,
    #     applymath$use_plot_breaks_labels,
    #     applymath$use_plot_limits,
    #     applymath$use_y_limits,
    #     applymath$legend_space
    #   )
    # })
  })
  
  
  # renders check box
  output$fileNamesCommon <- renderUI({
    req(input$file$name)
    checkboxGroupInput(
      "checkGroupCommon",
      label = h3("Plot on/off"),
      choices = first_file(),
      selected = "!NULL"
    )
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
        #submitButton("Update View"),
        uiOutput("rangeBin")
        
      ))
      ),
      mainPanel(tabsetPanel(
        tabPanel(
          "Tab 1",
          h4("Table Plot"),
          plotOutput("plot")
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