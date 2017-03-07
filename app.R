# Created by Benjamin Erickson BBErickson@gmail.com

# setwd("~/Desktop/BenToolsTab/dockerTest")
# setwd("~/BenTools gh/BenTools_shiny")
# runApp("BenTools_shiny")

source("helper.R")

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 20MB.
options(shiny.maxRequestSize = 20 * 1024 ^ 2)

# server ----
server <- function(input, output, session) {
  
  LIST_DATA <- list(
    table_file = list(),
    # [[]] gene X1 X2 ...
    gene_file = list(),
    # holds $common genes from files and $gene file(s)
    gene_info = list(),
    # for holding gene file info in a list of lists, a set for $common and each $gene file(s) [c("dot", "line", "color", plot?, NickName,nrom)]
    clust = list(), # Cluster holder
    x_plot_range = c(0, 0, 0, 0),
    STATE = 0
  )      
  Y_Axis_Lable <- NULL
  Lines_Lables_List <- NULL
  
  
  # load file, save data and set up info ----
  reactive_values <- reactiveValues(Make_Data_Frame = list(NULL),
                                    Apply_Math = NULL,
                                    Plot_Options = NULL)
  
  # loads file(s) 
  first_file <- reactive({
    req(input$tablefile$datapath)
    print("load file")
    # add warnings for total size of LIST_DATA
    LIST_DATA <<-
      LoadTableFile(input$tablefile$datapath, input$tablefile$name, LIST_DATA)
    names(LIST_DATA$table_file)
  })
  
  # records check box on/off for common list and builds data frame and info for plot
  observe({
    input$checkGroupCommon
    if (length(LIST_DATA$table_file) > 0) {
    print("checkbox on/off")
    LIST_DATA$gene_info <<-
      CheckBoxOnOff("common", input$checkGroupCommon, LIST_DATA$gene_info)
    reactive_values$Make_Data_Frame <- MakeDataFrame(LIST_DATA)
    reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      if (!is.null(isolate(reactive_values$Make_Data_Frame[[1]]))) {
        reactive_values$Apply_Math <- ApplyMath(isolate(reactive_values$Make_Data_Frame), isolate(input$myMath),isolate(reactive_values$Plot_Options))
      }
  }
  })
  
  observe({
    input$myMath
    if (length(LIST_DATA$table_file) > 0) {
    print("math")
    if (!is.null(isolate(reactive_values$Make_Data_Frame[[1]]))) {
      reactive_values$Apply_Math <- ApplyMath(isolate(reactive_values$Make_Data_Frame), input$myMath,isolate(reactive_values$Plot_Options))
      #print(GGplotLineDot(isolate(reactive_values$Apply_Math))) # renderplot
      output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
    }
    }
  })
  
  # set up and control check box actions common ----
  
  # renders check box when file is loaded
  output$fileOnOffCommon <- renderUI({
    req(first_file())
    print("render checkbox")
    checkboxGroupInput(
      "checkGroupCommon",
      label = h3("Plot on/off"),
      choices = first_file(),
      selected = unique(c(
        sapply(LIST_DATA$gene_info$common, "[[", 5),
        last(first_file())))
    )
  })
  
  # plots when acction button is pressed
  observe({
    req(input$myplot)
    print("plot button")
    if (!is.null(isolate(reactive_values$Apply_Math))) {
      Y_Axis_Lable <<- YAxisLable()
      Lines_Lables_List <<- LinesLablesList(use_pos_plot_ticks = c(LIST_DATA$x_plot_range[1:2]),
                                            use_label_plot_ticks = c(LIST_DATA$x_plot_range[1:2]))
      #print(GGplotLineDot(isolate(reactive_values$Apply_Math))) # renderplot
      output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
    }
  })
  
  # set up and control bin range action ----
  
  # Specification of range within an interval to plot
  output$rangeBin <- renderUI({
    req(input$tablefile$name)
    print("render slider")
    sliderInput( #  try and only render 1 time? how to update min max and values?
      "plotBinRange",
      label = h3("Plot Range:"),
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = LIST_DATA$x_plot_range[3:4]
    )
  })
  
  # plots when range bin slider is triggered
  observe({
    input$plotBinRange
    if (length(LIST_DATA$table_file) > 0) {
      if(LIST_DATA$STATE == 1) {
    print("slider")
    LIST_DATA$x_plot_range <<- c(LIST_DATA$x_plot_range[1:2], input$plotBinRange)
    if (!is.null(isolate(reactive_values$Make_Data_Frame[[1]]))) {
      #print(GGplotLineDot(isolate(reactive_values$Apply_Math))) # renderplot
      output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
    }
      } else {
        LIST_DATA$STATE <<- 1
      }
    }
  })
  
  # common file plot options ----
  output$fileOptionsCommon <- renderUI({
    req(input$tablefile$name)
    print("common options")
    tagList(
      radioButtons("commonOptionsSelect",
                   "choose",
                   choices = first_file())
    )
    
  })
  
  # quick color set change
  observe({
    col <- input$kbrewer
    if (!is.null(isolate(reactive_values$Apply_Math))) {
    print("kbrewer")
    kListColorSet <<- brewer.pal(8, col)
    lapply(names(LIST_DATA$gene_info), function(i) {
      lapply(seq_along(LIST_DATA$gene_info[[i]]), function(j) {
        color_safe <- j %% length(kListColorSet)
        if (color_safe == 0) {
          color_safe <- 1
        }
        LIST_DATA$gene_info[[i]][[j]][4] <<- kListColorSet[color_safe]
      })
    })
    reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
    
    output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
    }
  })
  
  # update color shown
  observe({
    input$kbrewer
    my_sel <- input$commonOptionsSelect
    if(is.null(my_sel)){
      my_sel <- names(LIST_DATA$gene_info$common)[1]
    }
    
    if (!is.null(isolate(reactive_values$Apply_Math))) {
    print("color update")
    updateColourInput(session, "col", value = paste(LIST_DATA$gene_info$common[[my_sel]][4]))
    updateTextInput(session, "rgbtohex",value = RgbToHex(my_hex = input$col))  
    }
  })
  
  observe({
    req(input$myrgb)
    print("color rgb")
    if (!is.null(isolate(reactive_values$Apply_Math))) {
      updateColourInput(session, "col", value = RgbToHex(my_rgb = isolate(input$rgbtohex)))
    }
  })
  
  # update color selected
  observe({
    input$col
    if (!is.null(isolate(reactive_values$Apply_Math))) {
      if(input$col != LIST_DATA$gene_info$common[[isolate(input$commonOptionsSelect)]][4]){
        print("color new")
        LIST_DATA$gene_info$common[[isolate(input$commonOptionsSelect)]][4] <<- input$col
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
      }
    }
  })
}

# UI -----
ui <- dashboardPage(
  dashboardHeader(title = "Ben Tools"),
  dashboardSidebar(sidebarMenu(
    id = "tabs",
    menuItem("Load Data", tabName = "loaddata", icon = icon("file")),
    
    menuItem("Plot", tabName = "mainplot", icon = icon("area-chart"))
  )),
  dashboardBody(tabItems(
    # load data tab
    tabItem(tabName = "loaddata",
            fluidRow(
              box(
                width = 4,
                fileInput(
                  "tablefile",
                  label = "Load table file",
                  accept = c('.table'),
                  multiple = TRUE
                ),
                fileInput("genefile",
                          label = "Load gene list",
                          accept = c('.txt')),
                fileInput(
                  "colorfile",
                  label = "Load color list",
                  accept = c('.color.txt')
                )
              ),
              box(
                colourInput("col", "Select color HEX"),
                textInput("rgbtohex", "RGB"),
                actionButton("myrgb", "Update HEX color")
                
              )
            )),
    
    # First tab content
    tabItem(tabName = "mainplot",
            fluidRow(
              box(width = 12, plotOutput("plot"))
              ),
            fluidRow(
              box(width = 4,
                uiOutput("rangeBin")
              ),
              box(width = 4,
                uiOutput("fileOptionsCommon"),
                uiOutput("fileOnOffCommon"),
                actionButton("myplot", "Update Plot")
              ),
              box(
                title = "math", collapsed = TRUE,
                collapsible = T,
                width = 3,
                radioButtons(
                  "myMath",
                  "Set math function",
                  choices = c("mean", "sum", "median", "var"),
                  selected = "mean"
                )
              ),
              box(width= 3, selectInput("kbrewer", "quick color set", choices = kBrewerList, selected = kBrewerList[3]))
            )
            )
  ))
)

# exicute ----
shinyApp(ui = ui, server = server)