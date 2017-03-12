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
  
  # Globals and reacive values ----
  reactive_values <- reactiveValues(
    Y_Axis_Lable = NULL,
    Lines_Lables_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL
  )
  
  # show hide checkbox and action button ----
  observe({
    req(first_file())
    toggle("checkboxonoff",condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0))
    toggle("selectgenelistonoff",condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0))
    toggle("actionmyplot",condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] == 2 ))
    })
  
  # loads data file(s) ----
  first_file <- reactive({
    req(input$filetable$datapath)
    isolate({
      print("load file")
      # add warnings for total size of LIST_DATA
      LIST_DATA <<-
        LoadTableFile(input$filetable$datapath,
                      input$filetable$name,
                      LIST_DATA)
      if (LIST_DATA$STATE[1] == 0) {
        print("1st slider and plot lines Ylable")
        reactive_values$Y_Axis_Lable <- YAxisLable()
        reactive_values$Lines_Lables_List <-
          LinesLablesList(
            use_pos_plot_ticks = c(LIST_DATA$x_plot_range[1:2]),
            use_label_plot_ticks = c(LIST_DATA$x_plot_range[1:2])
          )
        updateSliderInput(
          session,
          "sliderplotBinRange",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
        LIST_DATA$STATE[1] <<- 1
      }
      show("hidemainplot")
      show("startoff")
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
  observe({
    req(first_file())
    isolate({
      print("update load file")
      updateRadioButtons(
        session,
        "radiodataoption",
        choices = first_file(),
        selected = last(first_file())
      )
      updateCheckboxGroupInput(
        session,
        "checkboxonoff",
        choices = first_file(),
        selected = c(sapply(LIST_DATA$gene_info[[input$selectgenelistonoff]], "[[", 5))
      )
    })
  })
  
  # update when gene list is loaded ----
  observe({
    req(gene_file())
    isolate({
      print("select gene list options")
      updateSelectInput(session,
                        "selectgenelistoptions",
                        choices = names(LIST_DATA$gene_info))
      updateSelectInput(session,
                        "selectgenelistonoff",
                        choices = names(LIST_DATA$gene_info))
    })
  })
  
  # update desplay selected item info ----
  observe({
    my_sel <- input$radiodataoption
    isolate({
      req(first_file())
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
  })
  
  # triggers update on changing gene list ----
  observe({
    input$selectgenelistoptions
    isolate({
      req(first_file())
      print("update options on gene list change")
      updateRadioButtons(session, "radiodataoption", selected = input$radiodataoption)
    })
  })
  
  # record new nickname ----
  observe({
    req(input$actionoptions)
    isolate({
      req(first_file())
      print("new nickname")
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"] <<-
        input$textnickname
      LIST_DATA$table_file[[input$radiodataoption]]["set"] <<-
        paste(input$textnickname)
      if(LIST_DATA$STATE[1] == 1){
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
            reactive_values$Plot_Options
          )
        }
      }
    })
  })
  
  # records new dot options ----
  observe({
    input$selectdot
    isolate({
      if (!is.null(names(LIST_DATA$gene_info))) {
        if (input$selectdot != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mydot"]) {
          print("new dot")
          LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mydot"] <<-
            input$selectdot
          if(LIST_DATA$STATE[1] == 1 & !is.null(reactive_values$Apply_Math)) {
              reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
              reactive_values$Plot_controler <-
                GGplotLineDot(
                  reactive_values$Apply_Math,
                  input$sliderplotBinRange,
                  reactive_values$Plot_Options
                )
            }
          }
      }
    })
  })
  
  # records new line options ----
  observe({
    input$selectline
    isolate({
      if (!is.null(names(LIST_DATA$gene_info))) {
        if (input$selectline != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["myline"]) {
          print("new line")
          LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["myline"] <<-
            input$selectline
          if(LIST_DATA$STATE[1] == 1 & !is.null(reactive_values$Apply_Math)) {
            reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
            reactive_values$Plot_controler <-
              GGplotLineDot(
                reactive_values$Apply_Math,
                input$sliderplotBinRange,
                reactive_values$Plot_Options
              )
          }
        }
      }
    })
  })
  
  # update color based on rgb text input ----
  observe({
    req(input$actionmyrgb)
    isolate({
      print("color rgb")
      updateColourInput(session, "colourhex", value = RgbToHex(my_rgb = input$textrgbtohe))
    })
  })
  
  # update and save color selected ----
  observe({
    input$colourhex
    isolate({
      print("update text color")
      updateTextInput(session,
                      "textrgbtohex",
                      value = RgbToHex(my_hex = input$colourhex))
      if (!is.null(names(LIST_DATA$gene_info))) {
        if (input$colourhex != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"]) {
          print("color new")
          LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"] <<-
            input$colourhex
          if(LIST_DATA$STATE[1] == 1 & !is.null(reactive_values$Apply_Math)) {
            reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
            reactive_values$Plot_controler <-
              GGplotLineDot(
                reactive_values$Apply_Math,
                input$sliderplotBinRange,
                reactive_values$Plot_Options
              )
          }
        }
      }
    })
  })
  
  # update check box on off with selecting gene list ----
  observe({
    input$selectgenelistonoff
    isolate({
      req(first_file())
      print("update on off check box on gene list change")
      updateCheckboxGroupInput(session, "checkboxonoff",
                               selected = LIST_DATA$gene_info[[input$selectgenelistonoff]][[first_file()]]["onoff"])
    })
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
          toggle("actionmyplot",condition = (input$tabs == "mainplot"))
          disable("hidemainplot")
        } else {
          print("just updating what is seen")
          LIST_DATA$STATE[2] <<- input$selectgenelistonoff
        }
    })
  })
  
  # plots when acction button is pressed ----
  observe({
    req(input$actionmyplot)
    isolate({
      req(first_file())
      print("plot button")
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
            reactive_values$Plot_Options
          )
        enable("hidemainplot")
      } else{
        disable("hidemainplot")
      }
      hide("actionmyplot")
      LIST_DATA$STATE[1] <<- 1
    })
  })
  
  # renders plot ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  
  #update plot with math selection ----
  observe({
    input$myMath
    isolate({
      req(first_file())
      print("math")
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
              reactive_values$Plot_Options
            )
        }
    })
  })
  
  
  #plots when range bin slider is triggered ----
  observe({
    input$sliderplotBinRange
    isolate({
      req(first_file())
      print("slider")
      if (!is.null(reactive_values$Apply_Math)) {
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options
          )
      }
    })
  })
  
  # quick color set change ----
  observe({
    col <- input$kbrewer
    isolate({
      req(first_file())
      print("kbrewer")
      kListColorSet <<- brewer.pal(8, col)
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
        if(!is.null(reactive_values$Apply_Math)) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
              reactive_values$Plot_Options
            )
        }
      }
    })
  })
  
  
}

# UI -----
ui <- dashboardPage(
  dashboardHeader(title = "Ben Tools"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Load Data", tabName = "loaddata", icon = icon("file")),
      
     menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
     hidden( 
      checkboxGroupInput("checkboxonoff", h3("Plot on/off"), choices = "Load data file"),
      selectInput("selectgenelistonoff", "Select Gene list", choices = "common"),
      actionButton("actionmyplot", "Update Plot")
     )
    )
  ),
  dashboardBody(useShinyjs(),
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
                              fileInput("filegene",
                                        label = "Load gene list",
                                        accept = c('.txt')),
                              fileInput(
                                "filecolor",
                                label = "Load color list",
                                accept = c('.color.txt')
                              )
                            ),
                            hidden(div(id = "startoff", box(
                              selectInput("selectgenelistoptions", "Select Gene list", choices = "common"),
                              radioButtons("radiodataoption", "select data", choices = "Load data file"),
                              selectInput("selectdot", "Select dot type", choices = kDotOptions),
                              selectInput("selectline", "Select line type", choices = kLineOptions),
                              textInput("textnickname", "Nick Name"),
                              actionButton("actionoptions", "Update nick name")
                            ),
                            box(
                              colourInput("colourhex", "Select color HEX"),
                              textInput("textrgbtohex", "RGB"),
                              actionButton("actionmyrgb", "Update HEX color")
                            )))
                          )),
                  
                  # main plot tab
                  tabItem(tabName = "mainplot",
                          fluidRow(box(
                            width = 12, plotOutput("plot")
                          )),
                          hidden( div( id = "hidemainplot",  fluidRow(
                            box(
                              width = 4,
                              sliderInput(
                                "sliderplotBinRange",
                                label = h3("Plot Range:"),
                                min = 0,
                                max = 80,
                                value = c(0, 80)
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
                          ))))
                ))
)

# exicute ----
shinyApp(ui = ui, server = server)