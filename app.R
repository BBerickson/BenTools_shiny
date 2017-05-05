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
  # remove on non-local deployment
  session$onSessionEnded(stopApp)
  
  # Globals and reacive values ----
  
  # LIST_DATA <- list(
  #   table_file = list(),
  #   # [[]] gene X1 X2 ...
  #   gene_file = list(),
  #   # holds $common genes from files and $gene file(s)
  #   gene_info = list(),
  #   # for holding gene file info in a list of lists, a set for $common and each $gene file(s) [c("set", "dot", "line", "color", plot?, nrom)]
  #   clust = list(), # Cluster holder
  #   x_plot_range = c(0, 0),
  #   STATE = c(0, "common", 0) # flow control, gene list flow control
  # )      
  
  reactive_values <- reactiveValues(
    Y_Axis_Lable = NULL,
    Y_Axis_numbers = NULL,
    Lines_Lables_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL,
    Y_Axis_plot = NULL,
    mynorm = "none"
  )
  
  # show hide checkbox and action button ----
  observeEvent(input$tabs, {
    req(first_file())
    toggle(
      "showpicker",
      condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0)
      )
    toggle(
      "showpickersort",
      condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0 & LIST_DATA$STATE[3] != 0)
    )
    toggle(
      "selectlineslablesshow",
      condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0)
    )
    toggle("actionmyplotshow",
           condition = (input$tabs == "mainplot" & LIST_DATA$STATE[1] == 2))
  })
  
  # loads data file(s) ----
  first_file <- reactive({
    req(input$filetable$datapath)
    isolate({
      print("load file")
      # add warnings for total size of LIST_DATA TODO
      disable("startoff")
      disable("hidemainplot")
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...', value = 0, style = "old", {
      LIST_DATA <<-
        LoadTableFile(input$filetable$datapath,
                      input$filetable$name,
                      LIST_DATA)
                   })
      updateSelectInput(session,
                        "selectgenelistoptions",
                        choices = names(LIST_DATA$gene_info), selected = LIST_DATA$STATE[2])
      if (LIST_DATA$STATE[1] == 0) {
        show("filegene")
        show("checkboxconvert")
        show("downloadGeneList")
        show("filecolor")
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
        shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
        LIST_DATA$STATE[1] <<- 1
      }
      enable("startoff")
      enable("hidemainplot")
      reset("filetable")
      names(LIST_DATA$table_file)
    })
  })
  
  # loads gene list file ----
  observeEvent(input$filegene,{
    print("load gene file")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
    # load info, update select boxes, switching works and chaning info and ploting
    LIST_DATA <<- LoadGeneFile(input$filegene$datapath,
                               input$filegene$name,
                               LIST_DATA, input$checkboxconvert)
                 })
    updateSelectInput(session,
                      "selectgenelistoptions",
                      choices = names(LIST_DATA$gene_info), selected = LIST_DATA$STATE[2])
    onoff <- NULL
    for(i in sapply(strsplit(names(LIST_DATA$gene_file), "\nn ="), "[[",1)){
      onoff <- c(onoff,paste(i, names(LIST_DATA$table_file),sep = "-"))
    }
    
    TF <- c(sapply(names(LIST_DATA$gene_info), 
                   function(i) sapply(LIST_DATA$gene_info[[i]], "[[",5) != 0))
    mycolors <- paste("color", c(sapply(names(LIST_DATA$gene_info), function(i) sapply(LIST_DATA$gene_info[[i]], "[[",4))), sep = ":")
    updatePickerInput(session, "pickeronoff",
                      choices = onoff, 
                      selected = onoff[TF],
                      choicesOpt = list(style = mycolors)
    )
    reset("filegene")
    # add warnings for total size of LIST_DATA

  })
  
  # loads color file ----
  observeEvent(input$filecolor,{
    req(first_file())
    my_sel <- input$radiodataoption
    my_list <- input$selectgenelistoptions
    print("load color file")
    # load info, update select boxes, switching works and chaning info and ploting
    LIST_DATA <<- LoadColorFile(input$filecolor$datapath,
                               LIST_DATA)
    updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mycol"]))
    if (LIST_DATA$STATE[1] == 1 &
        !is.null(reactive_values$Apply_Math)) {
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options, 
          input$sliderplotYRange, 
          reactive_values$Lines_Lables_List, 
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
    }
    reset("filecolor")
    })
  
  # update when data file is loaded ----
  observeEvent(first_file(), {
    print("update load file")
    updateAwesomeRadio(
      session,
      "radiodataoption",
      choices = first_file(),
      selected = last(first_file())
    )
    onoff <- NULL
    for(i in sapply(strsplit(names(LIST_DATA$gene_file), "\nn ="), "[[",1)){
      onoff <- c(onoff,paste(i, names(LIST_DATA$table_file),sep = "-"))
    }
    
    TF <- c(sapply(names(LIST_DATA$gene_info), 
                   function(i) sapply(LIST_DATA$gene_info[[i]], "[[",5) != 0))
    mycolors <- paste("color", c(sapply(names(LIST_DATA$gene_info), function(i) sapply(LIST_DATA$gene_info[[i]], "[[",4))), sep = ":")
    updatePickerInput(session, "pickeronoff",
                      choices = onoff, 
                      selected = onoff[TF],
                      choicesOpt = list(style = mycolors)
                      )
  })
  
  # update desplay selected item info ----
  observeEvent(c(input$radiodataoption, input$selectgenelistoptions), {
    req(first_file())
    my_sel <- input$radiodataoption
    my_list <- input$selectgenelistoptions
    print("options update")
    updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mycol"]))
    updateTextInput(session,
                    "textnickname",
                    value = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["set"]))
    updateNumericInput(session, 
                       "normfactor", 
                       value = as.numeric(LIST_DATA$gene_info[[my_list]][[my_sel]]["rnorm"]))
    updateSelectInput(session,
                      "selectdot",
                      selected = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["mydot"]))
    updateSelectInput(session,
                      "selectline",
                      selected = paste(LIST_DATA$gene_info[[my_list]][[my_sel]]["myline"]))
    if(my_list == names(LIST_DATA$gene_info)[1]){
      enable("normfactor")
      disable("downloadGeneList")
    } else {
      disable("normfactor")
      enable("downloadGeneList")
    }
  })
  
  # saves gene list ----
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      paste(strsplit(input$selectgenelistoptions, " ")[[1]][1], ".txt", sep = "")
    },
    content = function(file) {
      new_comments <- paste("#", Sys.Date())
      tt <- c(new_comments, LIST_DATA$gene_file[[input$selectgenelistoptions]]$full)
      write_lines(tt, file)
    }
  )
  
  # record new nickname and norm factor ----
  observeEvent(input$actionoptions, {
    req(first_file())
    if (!is.na(input$normfactor) & !input$normfactor %in% c(0,1) & input$normfactor != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"]) {
      print("norm")
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"] <<- 
        input$normfactor
      LIST_DATA$table_file[[input$radiodataoption]] <<-
        mutate(LIST_DATA$table_file[[input$radiodataoption]],
               score = score / as.numeric(input$normfactor))
      if (LIST_DATA$STATE[1] == 1) {
        reactive_values$Apply_Math <-
          ApplyMath(LIST_DATA,
                    input$myMath,
                    input$checkboxrgf,
                    input$checkboxrf,
                    input$numericnormbin)
      }
    } else if(!is.na(input$normfactor)){
      updateNumericInput(session, 
                         "normfactor", 
                         value = as.numeric(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"]))
    }

    if(nchar(input$textnickname)>0){
      print("new nickname")
      oldnickname <- paste(gsub("(.{17})", "\\1\n", input$selectgenelistoptions), 
                           gsub("(.{17})", "\\1\n", 
                                LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"]), sep = '\n')
      
    LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"] <<-
      input$textnickname
    LIST_DATA$table_file[[input$radiodataoption]]["set"] <<-
      input$textnickname
    } else {
      updateTextInput(session,
                      "textnickname",
                      value = paste(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"]))
    }
    
    if (LIST_DATA$STATE[1] == 1) {
      if (!is.null(reactive_values$Apply_Math)) {
      newnickname <- paste(gsub("(.{17})", "\\1\n", input$selectgenelistoptions), 
                           gsub("(.{17})", "\\1\n", input$textnickname), sep = '\n')
      reactive_values$Apply_Math <- reactive_values$Apply_Math %>% 
        mutate(set = replace(set, set == oldnickname, newnickname))
      
        reactive_values$Plot_Options <- reactive_values$Plot_Options %>% 
          mutate(set = replace(set, set == oldnickname, newnickname))
        
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options, 
            input$sliderplotYRange, 
            reactive_values$Lines_Lables_List, 
            input$checkboxsmooth,
            reactive_values$Y_Axis_Lable
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
              reactive_values$Plot_Options, 
              input$sliderplotYRange, 
              reactive_values$Lines_Lables_List, 
              input$checkboxsmooth,
              reactive_values$Y_Axis_Lable
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
              reactive_values$Plot_Options, 
              input$sliderplotYRange, 
              reactive_values$Lines_Lables_List, 
              input$checkboxsmooth,
              reactive_values$Y_Axis_Lable
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
    req(first_file())
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
              reactive_values$Plot_Options, 
              input$sliderplotYRange, 
              reactive_values$Lines_Lables_List, 
              input$checkboxsmooth,
              reactive_values$Y_Axis_Lable
            )
        }
      }
    }
  })
    
  #records check box on/off ----
  observeEvent(input$pickeronoff, ignoreNULL = FALSE,{
      req(first_file())
        print("checkbox on/off")
        checkboxonoff <- list()
        for(i in input$pickeronoff){
          tt <- strsplit(sub("-"," ", i)," ")[[1]]
          selectgenelistonoff <- grep(tt[1], names(LIST_DATA$gene_file),value = T)
          checkboxonoff[[selectgenelistonoff]] <- c(checkboxonoff[[selectgenelistonoff]], tt[2])
        }
        LIST_DATA$gene_info <<-
          CheckBoxOnOff(checkboxonoff,
                        LIST_DATA$gene_info)
        LIST_DATA$STATE[1] <<- 2
        toggle("actionmyplotshow", condition = (input$tabs == "mainplot"))
        disable("selectlineslablesshow")
        disable("hidemainplot")
  })
  
  # plots when action button is pressed ----
  observeEvent(input$actionmyplot, {
    req(first_file())
    print("plot button")
    reactive_values$Apply_Math <-
      ApplyMath(LIST_DATA,
                input$myMath,
                input$checkboxrgf,
                input$checkboxrf,
                input$numericnormbin)
    if (!is.null(reactive_values$Apply_Math)) {
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      enable("hidemainplot")
      enable("selectlineslablesshow")
    } else{
      disable("hidemainplot")
      disable("selectlineslablesshow")
    }
    hide("actionmyplotshow")
    LIST_DATA$STATE[1] <<- 1
  })
  
  # updates y axis limits
  observeEvent(reactive_values$Apply_Math,{
    print("upate y axix on new math")
  test <- reactive_values$Y_Axis_numbers
    reactive_values$Y_Axis_numbers <- MyXSetValues(reactive_values$Apply_Math, 
                                                   input$sliderplotBinRange)
    
    updateSliderInput(session,
                      "sliderplotYRange",
                      min = reactive_values$Y_Axis_numbers[3],
                      max = reactive_values$Y_Axis_numbers[4],
                      value = reactive_values$Y_Axis_numbers[1:2],
                      step = ((reactive_values$Y_Axis_numbers[4] - 
                                 reactive_values$Y_Axis_numbers[3])/20))
    # Forces update if y Asis values stay the same
    if(sum(test) == sum(reactive_values$Y_Axis_numbers)){
      reactive_values$Y_Axis_plot <- input$actionmyplot[1]
    }
    
  })
  
  # renders plot ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  
  # updates applymath ----
  observeEvent(c(input$myMath, input$checkboxrgf, input$numericnormbin, input$checkboxrf),{
    req(first_file())
    if(is.na(input$numericnormbin)){
      updateNumericInput(session, "numericnormbin", value = 0)
    } else if(input$numericnormbin < 0){
      updateNumericInput(session, "numericnormbin", value = 0)
    } else if(input$numericnormbin > LIST_DATA$x_plot_range[2]){
      updateNumericInput(session, "numericnormbin", value = LIST_DATA$x_plot_range[2])
    } else {
      updateNumericInput(session, "numericnormbin", value = round(input$numericnormbin))
    }
    
    if(input$checkboxrgf & reactive_values$mynorm == "checkboxrf"){               
      updateCheckboxInput(session, "checkboxrf", value = FALSE)
      reactive_values$mynorm <- "checkboxrgf"
    } else if(input$checkboxrf & reactive_values$mynorm == "checkboxrgf"){
      updateCheckboxInput(session, "checkboxrgf", value = FALSE)
      updateNumericInput(session, "numericnormbin", value = 0)
      reactive_values$mynorm <- "checkboxrf"
    } else if(input$checkboxrf & reactive_values$mynorm == "numericnormbin"){
      updateNumericInput(session, "numericnormbin", value = 0)
      reactive_values$mynorm <- "checkboxrf"
    } else if(input$numericnormbin > 0 & reactive_values$mynorm == "checkboxrf"){
      updateCheckboxInput(session, "checkboxrf", value = FALSE)
      reactive_values$mynorm <- "numericnormbin"
    } else if(input$checkboxrgf){
      reactive_values$mynorm <- "checkboxrgf"
    } else if(input$checkboxrf){
      reactive_values$mynorm <- "checkboxrf"
    } else if(input$numericnormbin > 0){
      reactive_values$mynorm <- "numericnormbin"
    } else {
      reactive_values$mynorm <- "none"
    }
    reactive_values$Y_Axis_Lable <- YAxisLable(input$myMath, input$checkboxrf, input$checkboxrgf, input$numericnormbin, input$checkboxsmooth)
    if(LIST_DATA$STATE[1]==1){
    print("apply math")
    reactive_values$Apply_Math <-
      ApplyMath(LIST_DATA,
                input$myMath,
                input$checkboxrgf,
                input$checkboxrf,
                input$numericnormbin)
    }
  })
 
  #plots when bin slider or y slider is triggered ----
  observeEvent(c(reactive_values$Lines_Lables_List, input$sliderplotBinRange, reactive_values$Y_Axis_plot, input$sliderplotYRange), {
    req(first_file())
    if (!is.null(reactive_values$Apply_Math)) {
      print("bin slider or L&L making ggplot")
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options, 
          input$sliderplotYRange, 
          reactive_values$Lines_Lables_List, 
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
    }
  })

  # quick lines and lables preset change #TODO finish update ----
  observeEvent(input$selectlineslables, {
    req(first_file())
    print("quick Lines & Lables")
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
    reactive_values$Y_Axis_Lable <- YAxisLable(input$myMath, input$checkboxrf, input$checkboxrgf, input$numericnormbin, input$checkboxsmooth)
    reactive_values$Plot_controler <-
      GGplotLineDot(
        reactive_values$Apply_Math,
        input$sliderplotBinRange,
        reactive_values$Plot_Options, 
        input$sliderplotYRange, 
        reactive_values$Lines_Lables_List, 
        input$checkboxsmooth,
        reactive_values$Y_Axis_Lable
      )
  })
  
  # quick color set change ----
  observeEvent(input$kbrewer, {
    req(first_file())
    print("kbrewer")
    kListColorSet <<- brewer.pal(8, input$kbrewer)
    common_name <- names(LIST_DATA$gene_info)[1]
    if (!is.null(LIST_DATA$gene_info[[1]])) {
      print("kbrewer update")
      lapply(names(LIST_DATA$gene_info), function(i) {
        lapply(seq_along(LIST_DATA$gene_info[[i]]), function(j) {
          color_safe <- j %% length(kListColorSet)
          if (color_safe == 0) {
            color_safe <- 1
          }
          if(i == common_name){
            LIST_DATA$gene_info[[i]][[j]][4] <<-
              kListColorSet[color_safe]
          } else {
            LIST_DATA$gene_info[[i]][[j]][4] <<-
              RgbToHex(my_hex = kListColorSet[color_safe], tint = T)
          }
        })
      })
      updateColourInput(session, "colourhex", value = paste(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"]))
      if (!is.null(reactive_values$Apply_Math)) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options, 
            input$sliderplotYRange, 
            reactive_values$Lines_Lables_List, 
            input$checkboxsmooth,
            reactive_values$Y_Axis_Lable
          )
      }
    }
  })
  
  # sort tool action ----
  observeEvent(input$actionsorttool,{
    req(first_file())
    print("sort tool")
    dtsort <- SortTop(LIST_DATA, input$pickeronoff, 
                      input$slidersortbinrange[1], 
                      input$slidersortbinrange[2], input$slidersortpercent, input$selectsorttop)
    print(dtsort)
    output$sorttable <- renderDataTable(dtsort,
                                        options = list(
                                          pageLength = 5,
                                          searching = FALSE
                                          #initComplete = I("function(settings, json) {alert('Done.');}")
                                          )
    )
    
    show("showpickersort")
  })
  
  # sort tool action ----
  observeEvent(input$actionsortquick,{
    req(first_file())
    print("quick sort")
  })
  
  
  
  shinyjs::addClass(selector = "body", class = "sidebar-collapse")
}

# UI -----
ui <- dashboardPage(
  dashboardHeader(title = "Ben Tools"),
  dashboardSidebar(sidebarMenu(
    id = "tabs",
    menuItem("Load Data", tabName = "loaddata", icon = icon("file")),
    menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
    hidden(div(style = "padding-left: 15%",
               id = "showpicker",
               pickerInput(inputId = "pickeronoff", width = "75%",
                           label = h4("Select what to plot"), 
                           choices = "Load data file",
                           multiple = T,
                           options = list(`actions-box` = TRUE,`selected-text-format` = "count > 0")
               ))),
    hidden(div(style = "padding-left: 15%",
               id = "showpickersort",
               pickerInput(inputId = "pickersorttop", width = "75%",
                           label = h4("Sort gene list"), 
                           choices = "Load data file",
                           multiple = T,
                           options = list(`actions-box` = TRUE,`selected-text-format` = "count > 1")
               ))),
    menuItem("Sort Tool", tabName = "sorttool", icon = icon("gears")),
    
    
    hr(style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
    hidden(div(id = "actionmyplotshow", style = "padding-left: 20%",
      actionButton("actionmyplot", "Update Plot", icon = icon("area-chart"),
                   style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
    hidden(div(id = "selectlineslablesshow", style = "padding-left: 10%",  
      selectInput("selectlineslables", width = "85%",
                  label = "quick set lines and lables", 
                  choices = c("543 bins 20,20,40","543 bins 10,10,10", "5' 1k 1k 80bins" ,"3'","4")
      )
    ))
  )),
  dashboardBody(
    useShinyjs(),
                tabItems(
                  # load data tab
                  tabItem(tabName = "loaddata", 
                          fluidRow(tags$style(".nav-tabs-custom .nav-tabs li.active a { background-color: transparent; border-color: #2e6da4; } "),
                            tabBox(title = "Load files", id = "loadfiles",
                              width = 4, 
                              tabPanel("Table", 
                                       div(style = "height: 200px",
                                           fluidRow(div(style = "padding: 5px 15px",
                                                        fileInput(
                                "filetable", 
                                label = "Load table file",
                                accept = c('.table'),
                                multiple = TRUE)
                              )),
                              helpText("load table file(s) or bedGraph file(s)")
                              )),
                              tabPanel(title = "Gene",  
                                       div(style = "height: 200px",
                                           fluidRow(div(style = "padding: 5px 15px",
                                                        hidden(fileInput(
                                "filegene",
                                label = "Load gene list",
                                accept = c('.txt')
                              ),
                              checkboxInput("checkboxconvert",
                                            "gene list partial matching,      !!!can be slow!!!", value = FALSE),
                              downloadButton("downloadGeneList", "Save Gene List")
                                       )
                                           )))
                              ),
                              
                              tabPanel(title = "Color",
                                       div(style = "height: 200px",
                                           fluidRow(div(style = "padding: 5px 15px",
                                       hidden(fileInput(
                                "filecolor",
                                label = "Load color list",
                                accept = c('.color.txt')
                              )))))
                              )
                            ),
                    
                            hidden(div(
                              id = "startoff",
                              box(title =  "Select Gene list", width = 8,status = "primary", solidHeader = T,
                                selectInput("selectgenelistoptions", "", choices = "common"),
                                div(style = "height: 170px; overflow: auto",fluidRow(div(style = "padding: 0px 30px",
                                awesomeRadio("radiodataoption", "select data", choices = "Load data file")
                                )))
                                ),
                              box(title =  "Set Plot Options", width = 4, status = "primary", solidHeader = T,
                                selectInput("selectdot", "Select dot type", choices = kDotOptions),
                                selectInput("selectline", "Select line type", choices = kLineOptions),
                                
                                fluidRow(
                                  # column(1,tags$br()),
                                box(width = 12, status = "info", background = "light-blue",
                                    colourInput("colourhex", "Select color HEX"),
                                    tags$hr(),
                                    textInput("textrgbtohex", "RGB"),
                                    actionButton("actionmyrgb", "Update HEX color")
                                    ))
                              ),
                              box(title =  "Set Plot Options", width = 4, status = "primary", solidHeader = T,
                                textInput("textnickname", "Update Nickname"),
                                numericInput("normfactor", "Set norm factor, score/rpm", value = 1),
                                actionButton("actionoptions", "Update"),
                                helpText("Need to update Nickname and/or nrom factor")
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
                              
                              box(title = "Sliders", status = "primary", solidHeader = T,
                                  width = 7, collapsible = TRUE,
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
                              
                              box(style='padding:2px;', 
                                title = "math",
                                collapsed = F,
                                collapsible = T,
                                width = 2, status = "primary", solidHeader = T,
                                awesomeRadio(
                                  "myMath",label = 
                                    " ",
                                  choices = c("mean", "sum", "median", "var"),
                                  selected = "mean"
                                )
                              ),
                              box(title = "Normalization", status = "primary", solidHeader = T, 
                                  width = 3, collapsible = TRUE,
                                  checkboxInput("checkboxrf", label = "relative frequency"),
                                  checkboxInput("checkboxrgf", label = "relative gene frequency"),
                                  checkboxInput("checkboxsmooth", label = "smooth"),
                                  numericInput("numericnormbin", "Norm to bin", value = 0)
                              ),
                              box(title = "Lines and Labels", width = 9, status = "primary", solidHeader = T,
                                  collapsible = TRUE, collapsed = TRUE,
                                  div(style="padding:2px; display:inline-block", numericInput("numerictss", "TSS bin",value = 15, min = 0, max = 100)),
                                  div(style="padding:2px; display:inline-block", numericInput("numerictes", "TES bin",value = 45, min = 0, max = 100)),
                                  div(style="padding:2px; display:inline-block", numericInput("numericbinsize", 
                                                                                 "bp/bin",value = 100, min = 20, max = 1000, step = 5)),
                                  
                                  div(style="padding:2px; display:inline-block", numericInput("numericbody1", "5|4 bin",value = 20, min = 0, max = 100)),
                                  div(style="padding:2px; display:inline-block",numericInput("numericbody2", "4|3 bin",value = 40, min = 0, max = 100)),
                                  div(style="padding:2px; display:inline-block",numericInput("numericlabelspaceing", "every bin",value = 5, min = 0, max = 100)),
                                  div(style= "padding-left:33%", actionButton("actionlineslabels", "Update Lines and Lables"))
                              ),
                              box(title = "Quick Color Change",
                                width = 3, status = "primary", solidHeader = T,
                                selectInput(
                                  "kbrewer",
                                  "color brewer",
                                  choices = kBrewerList,
                                  selected = kBrewerList[6]
                                )
                              )
                            )
                          ))),
                  # main plot tab
                  tabItem(tabName = "sorttool",
                          
                          box(title = "Sort tool", status = "primary", solidHeader = T,
                            width = 8,
                            selectInput(
                              "selectsorttop",
                              "Sort Options",
                              choices = c("Top%", "Bottom%"),
                              selected = "Top%"
                            ),
                            sliderInput(
                              "slidersortpercent",
                              label = "% select:",
                              min = 1,
                              max = 100,
                              value = 80
                            ),
                            sliderInput(
                              "slidersortbinrange",
                              label = "Select Bin Range:",
                              min = 0,
                              max = 80,
                              value = c(0, 80)
                            ),
                            actionButton("actionsorttool", "Sort selected"),
                            actionButton("actionsortquick", "Quick% Sort selected"),
                            dataTableOutput('sorttable')
                          ))
                ))
)

# exicute ----
shinyApp(ui = ui, server = server)