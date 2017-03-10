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
    # for holding gene file info in a list of lists, a set for $common and each $gene file(s) [c("set", "dot", "line", "color", plot?, nrom)]
    clust = list(), # Cluster holder
    x_plot_range = c(0, 0, 0, 0),
    STATE = c(0, "common")
  )      
  Y_Axis_Lable <- NULL
  Lines_Lables_List <- NULL
  
  
  # load file, save data and set up info ----
  reactive_values <- reactiveValues(Make_Data_Frame = list(NULL),
                                    Apply_Math = NULL,
                                    Plot_Options = NULL)
  
  # loads data file(s) 
  first_file <- reactive({
    req(input$filetable$datapath)
    print("load file")
    # add warnings for total size of LIST_DATA
    LIST_DATA <<-
      LoadTableFile(input$filetable$datapath, input$filetable$name, LIST_DATA)
    names(LIST_DATA$table_file)
  })
  
  # loads data file(s) 
  gene_file <- reactive({
    req(input$filegene$datapath)
    print("load gene file")
    # add warnings for total size of LIST_DATA
    names(LIST_DATA$gene_info)
  })
  
  # update when data file is loaded
  observe({ 
    req(first_file())
    print("update load file")
    updateRadioButtons(session, "radiodataoption", choices = first_file(), selected = last(first_file()))
    updateCheckboxGroupInput(session,"checkboxonoff", choices = first_file(), 
                             selected = LIST_DATA$gene_info[[input$selectgenelistonoff]][[first_file()]]["onoff"])
    if(LIST_DATA$STATE[1] == 0){
      print("set slider")
      updateSliderInput(session, "sliderplotBinRange",
        min = LIST_DATA$x_plot_range[1],
        max = LIST_DATA$x_plot_range[2],
        value = LIST_DATA$x_plot_range[3:4]
      )
      LIST_DATA$STATE[1] <<- 1
    } else{
      print("slider already set")
    }
  })
  
  # update when gene list is loaded
  observe({ 
    req(gene_file())
    print("select gene list options")
    updateSelectInput(session, "selectgenelistoptions", choices = names(LIST_DATA$gene_info))
    updateSelectInput(session, "selectgenelistonoff", choices = names(LIST_DATA$gene_info))
  })
  
  
  # update desplay selected item info
  observe({
    my_sel <- input$radiodataoption
    req(isolate(first_file()))
    my_list <- isolate(input$selectgenelistoptions)
    print("options update")
    updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mycol"]))
    updateTextInput(session, "textrgbtohex",value = RgbToHex(my_hex = isolate(input$colourhex)))
    updateTextInput(session, "textnickname", value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["set"]))
    updateSelectInput(session, "selectdot", selected = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mydot"]))
    updateSelectInput(session, "selectline", selected = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["myline"]))
  })
  
  observe({
    req(isolate(first_file()))
    input$selectgenelistoptions
    print("update options on gene list change")
    updateRadioButtons(session, "radiodataoption", selected = isolate(input$radiodataoption))
  })
  
  # save new plot options
  observe({
    req(input$actionoptions)
    isolate(if (!is.null(names(LIST_DATA$gene_info))) {
      print("new options")
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"] <<- input$textnickname
      LIST_DATA$table_file[[input$radiodataoption]]["set"] <<- paste(input$textnickname)
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mydot"] <<- input$selectdot
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["myline"] <<- input$selectline
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
    })
  })
  
  # update color based on rgb input
  observe({
    req(input$actionmyrgb)
    print("color rgb")
    isolate(
    if (!is.null(reactive_values$Apply_Math)) {
      updateColourInput(session, "colourhex", value = RgbToHex(my_rgb = input$textrgbtohe))
    })
  })
  
  # update and save color selected
  observe({
    input$colourhex
    if (!is.null(isolate(reactive_values$Apply_Math))) {
      if(input$colourhex != LIST_DATA$gene_info[[input$selectgenelistoptions]][[isolate(input$radiodataoption)]]["mycol"]){
        print("color new")
        LIST_DATA$gene_info[[input$selectgenelistoptions]][[isolate(input$radiodataoption)]]["mycol"] <<- input$colourhex
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
      }
    }
  })
  
  
  # update check box on off with selecting gene list
  observe({
    input$selectgenelistonoff
    req(isolate(first_file()))
    print("update on off check box on gene list change")
    updateCheckboxGroupInput(session,"checkboxonoff", 
                             selected = LIST_DATA$gene_info[[input$selectgenelistonoff]][[first_file()]]["onoff"])
  })
  
  #records check box on/off for common list and builds data frame and info for plot
  observe({
    input$checkboxonoff
    req(isolate(first_file()))
    if(LIST_DATA$STATE[2] == isolate(input$selectgenelistonoff)){
      print("checkbox on/off")
      LIST_DATA$gene_info <<-
        CheckBoxOnOff(isolate(input$selectgenelistonoff), input$checkboxonoff, LIST_DATA$gene_info)
      reactive_values$Make_Data_Frame <- MakeDataFrame(LIST_DATA)
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      if (!is.null(isolate(reactive_values$Make_Data_Frame[[1]]))) {
        reactive_values$Apply_Math <- ApplyMath(isolate(reactive_values$Make_Data_Frame), isolate(input$myMath),isolate(reactive_values$Plot_Options))
      print(reactive_values$Plot_Options)
      print(GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))) # renderplot
      }
      print(isolate(reactive_values$Make_Data_Frame[[1]]))
      print(reactive_values$Apply_Math)
    } else {
      print("just updating what is seen")
      LIST_DATA$STATE[2] <<- isolate(input$selectgenelistonoff)
    }
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
  
  observe({
    input$myMath
    if (length(LIST_DATA$table_file) > 0) {
      print("math")
      if (!is.null(isolate(reactive_values$Make_Data_Frame[[1]]))) {
        reactive_values$Apply_Math <- ApplyMath(isolate(reactive_values$Make_Data_Frame), input$myMath,isolate(reactive_values$Plot_Options))
        output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
      }
    }
  })
  

  # plots when range bin slider is triggered
  # observe({
  #   input$sliderplotBinRange
  #   req(isolate(first_file()))
  #   print("slider")
  #   LIST_DATA$x_plot_range <<- c(LIST_DATA$x_plot_range[1:2], input$plotBinRange)
  #   if (!is.null(isolate(reactive_values$Make_Data_Frame[[1]]))) {
  #     #print(GGplotLineDot(isolate(reactive_values$Apply_Math))) # renderplot
  #     output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
  #   }
  # })
  
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
    updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[isolate(input$selectgenelistoptions)]][[isolate(input$radiodataoption)]]["mycol"]))
    updateTextInput(session, "textrgbtohex",value = RgbToHex(my_hex = input$colourhex))
    reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
    
    output$plot <- renderPlot({GGplotLineDot(isolate(reactive_values$Apply_Math), LIST_DATA$x_plot_range[3:4], isolate(reactive_values$Plot_Options))})
    }
  })
  
  # data and list plot options ----
  
  
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
                  "filetable",
                  label = "Load table file",
                  accept = c('.table'),
                  multiple = TRUE
                ),
                fileInput("filegene",
                          label = "Load gene list",
                          accept = c('.txt')),
                fileInput(
                  "filecolor",
                  label = "Load color list",
                  accept = c('.color.txt')
                )
              ),
              box(
                selectInput("selectgenelistoptions", "Select Gene list", choices = "common"),
                radioButtons("radiodataoption","select data", choices = "Load data file"),
                selectInput("selectdot", "Select dot type", choices = kDotOptions),
                selectInput("selectline", "Select line type", choices = kLineOptions),
                textInput("textnickname", "Nick Name"),
                actionButton("actionoptions", "Update nick name")
              ),
              box(
                colourInput("colourhex", "Select color HEX"),
                textInput("textrgbtohex", "RGB"),
                actionButton("actionmyrgb", "Update HEX color")
              )
            )),
    
    # main plot tab
    tabItem(tabName = "mainplot",
            fluidRow(
              box(width = 12, plotOutput("plot"))
              ),
            fluidRow(
              box(width = 4,
                  sliderInput("sliderplotBinRange",
                    label = h3("Plot Range:"),
                    min = 0,
                    max = 80,
                    value = c(0,80)
                  )
              ),
              box(width = 4,
                selectInput("selectgenelistonoff", "Select Gene list", choices = "common"),
                checkboxGroupInput("checkboxonoff", h3("Plot on/off"), choices = "Load data file"),
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