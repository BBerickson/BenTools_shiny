# Created by Benjamin Erickson BBErickson@gmail.com

# setwd("~/Desktop/BenToolsTab/dockerTest")
# setwd("~/BenTools gh/BenTools_shiny")
# runApp("BenTools_shiny")

source("helper.R")

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 20MB.
options(shiny.maxRequestSize = 30 * 1024 ^ 2)



# server ----
server <- function(input, output, session) {
  # Globals and reacive values ----
  reactive_values <- reactiveValues(
    Y_Axis_Lable = NULL,
    Lines_Lables_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL
  )
  
  # show hide checkbox and action button ----
  observeEvent(input$tabs, {
    req(first_file())
    toggle("checkboxonoff",
           condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0))
    toggle(
      "selectgenelistonoff",
      condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0)
    )
    toggle(
      "selectlineslables",
      condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0)
    )
    toggle("actionmyplot",
           condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] == 2))
  })
  
  # loads data file(s) ----
  first_file <- reactive({
    req(input$filetable$datapath)
    isolate({
      print("load file")
      # add warnings for total size of LIST_DATA TODO
      LIST_DATA <<-
        LoadTableFile(input$filetable$datapath,
                      input$filetable$name,
                      LIST_DATA)
      if (LIST_DATA$STATE[1] == 0) {
        show("hidemainplot")
        show("startoff")
        print("1st slider and plot lines Ylable")
        reactive_values$Y_Axis_Lable <- YAxisLable()
        reactive_values$Lines_Lables_List <- LinesLablesList()
        updateSliderInput(
          session,
          "sliderplotBinRange",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
        LIST_DATA$STATE[1] <<- 1
      }
      names(LIST_DATA$table_file)
    })
  })
  
  
  # loads data file(s) ----
  gene_file <- reactive({
    req(input$filegene$datapath)
    print("load gene file")
    # add warnings for total size of LIST_DATA
    names(LIST_DATA$gene_info)
  })
  
  # update when data file is loaded ----
  observeEvent(first_file(), {
    print("update load file")
    # prevents over flow of text ... needs some added padding #TODO
    addClass("radiodataoption",class = "ofhidden")
    addClass("checkboxonoff",class = "ofhidden")
    updateRadioButtons(
      session,
      "radiodataoption",
      choices = first_file(),
      selected = last(first_file())
    )
    updateCheckboxGroupInput(session,
                             "checkboxonoff",
                             choices = first_file(),
                             selected = c(sapply(LIST_DATA$gene_info[[input$selectgenelistonoff]], "[[", 5)))
  })
  
  # update when gene list is loaded ----
  observeEvent(gene_file(), {
    print("select gene list options")
    updateSelectInput(session,
                      "selectgenelistoptions",
                      choices = names(LIST_DATA$gene_info))
    updateSelectInput(session,
                      "selectgenelistonoff",
                      choices = names(LIST_DATA$gene_info))
  })
  
  # update desplay selected item info ----
  observeEvent(input$radiodataoption, {
    req(first_file())
    my_sel <- input$radiodataoption
    my_list <- input$selectgenelistoptions
    print("options update")
    updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mycol"]))
    updateTextInput(session,
                    "textnickname",
                    value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["set"]))
    updateSelectInput(session,
                      "selectdot",
                      selected = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mydot"]))
    updateSelectInput(session,
                      "selectline",
                      selected = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["myline"]))
  })
  
  # triggers update on changing gene list ----
  observeEvent(input$selectgenelistoptions, {
    req(first_file())
    print("update options on gene list change")
    updateRadioButtons(session, "radiodataoption", selected = input$radiodataoption)
  })
  
  # record new nickname ----
  observeEvent(input$actionoptions, {
    req(first_file())
    print("new nickname")
    LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"] <<-
      input$textnickname
    LIST_DATA$table_file[[input$radiodataoption]]["set"] <<-
      paste(input$textnickname)
    if (LIST_DATA$STATE[1] == 1) {
      reactive_values$Apply_Math <-
        ApplyMath(LIST_DATA,
                  input$myMath,
                  r_checkbox_gene_relative_frequency = 0)
      if (!is.null(reactive_values$Apply_Math)) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
             reactive_values$Plot_Options, input$sliderplotYRange, 
            reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
          )
      }
    }
  })
  
  # records new dot options ----
  observeEvent(input$selectdot, {
    if (!is.null(names(LIST_DATA$gene_info))) {
      if (input$selectdot != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mydot"]) {
        print("new dot")
        LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mydot"] <<-
          input$selectdot
        if (LIST_DATA$STATE[1] == 1 &
            !is.null(reactive_values$Apply_Math)) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
               reactive_values$Plot_Options, input$sliderplotYRange,
              reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
            )
        }
      }
    }
  })
  
  # records new line options ----
  observeEvent(input$selectline, {
    if (!is.null(names(LIST_DATA$gene_info))) {
      if (input$selectline != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["myline"]) {
        print("new line")
        LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["myline"] <<-
          input$selectline
        if (LIST_DATA$STATE[1] == 1 &
            !is.null(reactive_values$Apply_Math)) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
               reactive_values$Plot_Options, input$sliderplotYRange, 
              reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
            )
        }
      }
    }
  })
  
  # update color based on rgb text input ----
  observeEvent(input$actionmyrgb, {
    print("color rgb")
    updateColourInput(session, "colourhex", value = RgbToHex(my_rgb = input$textrgbtohex))
  })
  
  # update and save color selected ----
  observeEvent(input$colourhex, {
    print("update text color")
    updateTextInput(session,
                    "textrgbtohex",
                    value = RgbToHex(my_hex = input$colourhex))
    if (!is.null(names(LIST_DATA$gene_info))) {
      if (input$colourhex != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"]) {
        print("color new")
        LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"] <<-
          input$colourhex
        if (LIST_DATA$STATE[1] == 1 &
            !is.null(reactive_values$Apply_Math)) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
               reactive_values$Plot_Options, input$sliderplotYRange, 
              reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
            )
        }
      }
    }
  })
  
  # update check box on off with selecting gene list ----
  observeEvent(input$selectgenelistonoff, {
    req(first_file())
    print("update on off check box on gene list change")
    updateCheckboxGroupInput(session, "checkboxonoff",
                             selected = LIST_DATA$gene_info[[input$selectgenelistonoff]][[first_file()]]["onoff"])
  })
  
  #records check box on/off ----
  observe({
    input$checkboxonoff
    isolate({
      req(first_file())
      if (LIST_DATA$STATE[2] == input$selectgenelistonoff) {
        print("checkbox on/off")
        LIST_DATA$gene_info <<-
          CheckBoxOnOff(input$selectgenelistonoff,
                        input$checkboxonoff,
                        LIST_DATA$gene_info)
        LIST_DATA$STATE[1] <<- 2
        toggle("actionmyplot", condition = (input$tabs == "mainplot"))
        disable("selectlineslables")
        disable("hidemainplot")
      } else {
        print("just updating what is seen")
        LIST_DATA$STATE[2] <<- input$selectgenelistonoff
      }
    })
  })
  
  # plots when acction button is pressed ----
  observeEvent(input$actionmyplot, {
    req(first_file())
    print("plot button")
    reactive_values$Apply_Math <-
      ApplyMath(LIST_DATA,
                input$myMath,
                r_checkbox_gene_relative_frequency = 0)
    if (!is.null(reactive_values$Apply_Math)) {
      myYBinRange <- MyYSetValues(reactive_values$Apply_Math, input$sliderplotBinRange)
      updateSliderInput(session,
                        "sliderplotYRange",
                        min = myYBinRange[3],
                        max = myYBinRange[4],
                        value = myYBinRange[1:2],
                        step = ((myYBinRange[2]-myYBinRange[1])/20))
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      enable("hidemainplot")
      enable("selectlineslables")
    } else{
      disable("hidemainplot")
      disable("selectlineslables")
    }
    hide("actionmyplot")
    LIST_DATA$STATE[1] <<- 1
  })
  
  # renders plot ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  
  #update plot with math selection ----
  observeEvent(input$myMath, {
    req(first_file())
    print("math")
    reactive_values$Apply_Math <-
      ApplyMath(LIST_DATA,
                input$myMath,
                r_checkbox_gene_relative_frequency = 0)
    if (!is.null(reactive_values$Apply_Math)) {
      myYBinRange <- MyYSetValues(reactive_values$Apply_Math, input$sliderplotBinRange)
      updateSliderInput(session,
                        "sliderplotYRange",
                        min = myYBinRange[3],
                        max = myYBinRange[4],
                        value = myYBinRange[1:2],
                        step = ((myYBinRange[2]-myYBinRange[1])/20))
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
    }
  })
  
  
  #plots when range bin slider is triggered ----
  observeEvent(input$sliderplotBinRange, {
    req(first_file())
    if (!is.null(reactive_values$Apply_Math)) {
      print("bin slider")
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
           reactive_values$Plot_Options, input$sliderplotYRange, 
          reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
        )
    }
  })
  
  #plots when Y range slider is triggered ----
  observeEvent(input$sliderplotYRange, {
    req(first_file())
    if (!is.null(reactive_values$Apply_Math)) {
      print("Y slider")
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
           reactive_values$Plot_Options, input$sliderplotYRange, 
          reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
        )
    }
  })

  # quick lines and lables preset change #TODO finish update ----
  observeEvent(input$selectlineslables, {
    req(first_file())
    print("update lines and lables")
    myset <- LinesLablesPreSet(input$selectlineslables)
      updateNumericInput(session,"numericbody1", value = myset[1])
      updateNumericInput(session,"numericbody2", value = myset[2])
      updateNumericInput(session,"numerictss", value = myset[3])
      updateNumericInput(session,"numerictes", value = myset[4])
      updateNumericInput(session,"numericbinsize", value = myset[5])
      updateNumericInput(session,"numericlabelspaceing", value = myset[6])
   
    
    reactive_values$Lines_Lables_List <- 
      LinesLablesList(myset[1],
                      myset[2],
                      myset[3],
                      myset[4],
                      myset[5],
                    LIST_DATA$x_plot_range[2],
                    myset[6])
   
  })
  
  # action button update lines and lables ----
  observeEvent(input$actionlineslabels,{
    req(first_file())
    print("action lines and lables")
    myset <- c(input$numericbody1,
    input$numericbody2,
    input$numerictss,
    input$numerictes,
    input$numericbinsize,
    input$numericlabelspaceing)
    print(myset)
    myset[is.na(myset)] <- 0
    
    for(i in seq_along(myset)){
        if(myset[i] < 1){
          myset[i] <- 0
        } else if(i %in% c(1:4,6) & myset[i] > LIST_DATA$x_plot_range[2]){
          myset[i] <- LIST_DATA$x_plot_range[2]
        }
    }
    updateNumericInput(session,"numericbody1", value = myset[1])
    updateNumericInput(session,"numericbody2", value = myset[2])
    updateNumericInput(session,"numerictss", value = myset[3])
    updateNumericInput(session,"numerictes", value = myset[4])
    updateNumericInput(session,"numericbinsize", value = myset[5])
    updateNumericInput(session,"numericlabelspaceing", value = myset[6])
    
    reactive_values$Lines_Lables_List <- 
      LinesLablesList(myset[1],
                      myset[2],
                      myset[3],
                      myset[4],
                      myset[5],
                      LIST_DATA$x_plot_range[2],
                      myset[6])
    
  })
  
  # replot with smooth update ----
  observeEvent(input$checkboxsmooth, {
    req(first_file())
    print(input$checkboxsmooth)
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options, 
          input$sliderplotYRange, 
          reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
        )
  })
  
  # observe lines and labes update and update plot ----
  observeEvent(reactive_values$Lines_Lables_List,{
    req(input$actionmyplot)
    reactive_values$Plot_controler <-
      GGplotLineDot(
        reactive_values$Apply_Math,
        input$sliderplotBinRange,
        reactive_values$Plot_Options, input$sliderplotYRange, 
        reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
      )
  })
  
  # quick color set change ----
  observeEvent(input$kbrewer, {
    req(first_file())
    print("kbrewer")
    kListColorSet <<- brewer.pal(8, input$kbrewer)
    if (!is.null(LIST_DATA$gene_info[[1]])) {
      print("kbrewer update")
      lapply(names(LIST_DATA$gene_info), function(i) {
        lapply(seq_along(LIST_DATA$gene_info[[i]]), function(j) {
          color_safe <- j %% length(kListColorSet)
          if (color_safe == 0) {
            color_safe <- 1
          }
          LIST_DATA$gene_info[[i]][[j]][4] <<-
            kListColorSet[color_safe]
        })
      })
      updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"]))
      if (!is.null(reactive_values$Apply_Math)) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
             reactive_values$Plot_Options, input$sliderplotYRange, 
            reactive_values$Lines_Lables_List, use_smooth = input$checkboxsmooth
          )
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
    
    menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
    hidden(
      checkboxGroupInput("checkboxonoff", h3("Plot on/off"), choices = "Load data file"),
      selectInput("selectgenelistonoff", "Select Gene list", choices = "common"),
      actionButton("actionmyplot", "Update Plot"),
      selectInput("selectlineslables", 
                  label = "quick set lines and lables", 
                  choices = c("543 bins 20,20,40","543 bins 10,10,10", "5' 1k 1k 80bins" ,"3'","4")
      )
    )
  )),
  dashboardBody(useShinyjs(),
                inlineCSS(list(.ofhidden = "overflow: auto")),
                tabItems(
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
                              hidden(fileInput(
                                "filegene",
                                label = "Load gene list",
                                accept = c('.txt')
                              )),
                              hidden(fileInput(
                                "filecolor",
                                label = "Load color list",
                                accept = c('.color.txt')
                              ))
                            )
                            ,
                            hidden(div(
                              id = "startoff",
                              box(
                                selectInput("selectgenelistoptions", "Select Gene list", choices = "common"),
                                radioButtons("radiodataoption", "select data", choices = "Load data file"),
                                selectInput("selectdot", "Select dot type", choices = kDotOptions),
                                selectInput("selectline", "Select line type", choices = kLineOptions),
                                textInput("textnickname", "Nick Name"),
                                actionButton("actionoptions", "Update nick name")
                              ),
                              box(width = 3,
                                colourInput("colourhex", "Select color HEX"),
                                textInput("textrgbtohex", "RGB"),
                                actionButton("actionmyrgb", "Update HEX color")
                              )
                            ))
                          )),
                  
                  # main plot tab
                  tabItem(tabName = "mainplot",
                          fluidRow(box(
                            width = 12, plotOutput("plot")
                          )),
                          hidden(div(
                            id = "hidemainplot",  fluidRow(
                              
                              box(title = "Normalization", width = 3, collapsible = TRUE,
                                  checkboxInput("checkboxrf", label = "relative frequency"),
                                  checkboxInput("checkboxrgf", label = "relative gene frequency"),
                                  checkboxInput("checkboxsmooth", label = "smooth"),
                                  numericInput("numericnormbin", "Norm to bin", value = 0)),
                              
                              box(title = "Lines and Labels", width = 3, collapsible = TRUE, collapsed = TRUE,
                                  numericInput("numericbody1", "5|4 bin",value = 20),
                                  numericInput("numericbody2", "4|3 bin",value = 40),
                                  numericInput("numerictss", "TSS bin",value = 15),
                                  numericInput("numerictes", "TES bin",value = 45),
                                  numericInput("numericbinsize", "bp/bin",value = 100, min = 20, max = 1000, step = 5),
                                  numericInput("numericlabelspaceing", "every bin",value = 5),
                                  actionButton("actionlineslabels", "Update Lines and Lables")
                                
                              ),
                              box(title = "Sliders",
                                width = 6, collapsible = TRUE,
                                sliderInput(
                                  "sliderplotBinRange",
                                  label = "Plot Bin Range:",
                                  min = 0,
                                  max = 80,
                                  value = c(0, 80)
                                ),
                                sliderInput(
                                  "sliderplotYRange",
                                  label = "Plot Y hight:",
                                  min = 0,
                                  max = 1,
                                  value = c(0, 1)
                                )
                              ),
                              box(
                                title = "math",
                                collapsed = TRUE,
                                collapsible = T,
                                width = 3,
                                radioButtons(
                                  "myMath",
                                  "Set math function",
                                  choices = c("mean", "sum", "median", "var"),
                                  selected = "mean"
                                )
                              ),
                              box(
                                width = 3,
                                selectInput(
                                  "kbrewer",
                                  "quick color set",
                                  choices = kBrewerList,
                                  selected = kBrewerList[3]
                                )
                              )
                            )
                          )))
                ))
)

# exicute ----
shinyApp(ui = ui, server = server)