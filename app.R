# Created by Benjamin Erickson BBErickson@gmail.com

# setwd("~/Desktop/BenToolsTab/dockerTest")
# setwd("~/BenTools gh/BenTools_shiny")
# runApp("BenTools_shiny")
# runApp("app.R")

source("helper.R")

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 30MB.
options(shiny.maxRequestSize = 30 * 1024 ^ 2)



# server ----
server <- function(input, output, session) {
  # remove on non-local deployment
  session$onSessionEnded(stopApp)
  
  # Globals and reacive values ----
  
  reactive_values <- reactiveValues(
    Y_Axis_Lable = NULL,
    Y_Axis_numbers = NULL,
    Lines_Lables_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL,
    Picker_controler = NULL,
    Y_Axis_plot = NULL,
    mynorm = "none"
  )
  
  # change tab controls ----
  observeEvent(input$tabs, ignoreInit = TRUE, {
    print("switching tabs")
    toggle("showpicker",
           condition = (input$tabs == "mainplot" &
                          LIST_DATA$STATE[1] != 0))
    toggle("showsorttoolpicker",
           condition = (input$tabs == "sorttool" &
                          LIST_DATA$STATE[1] != 0))
    toggle(
      "showratiotoolpicker",
      condition = (input$tabs == "ratiotool" &
                     LIST_DATA$STATE[1] != 0)
    )
    toggle(
      "showclustertoolpicker",
      condition = (input$tabs == "clustertool" &
                     LIST_DATA$STATE[1] != 0)
    )
    
    if (input$tabs == "sorttool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectsortfile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectsortfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      if (sum(grepl("Sort\nn =", names(LIST_DATA$gene_file))) == 0) {
        updateSliderInput(
          session,
          "slidersortbinrange",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
      }
    }
    
    if (input$tabs == "ratiotool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectratiofile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(session,
                        "selectratiofile",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol)
      if (sum(grepl("Ratio_", names(LIST_DATA$gene_file))) == 0) {
        updateSliderInput(
          session,
          "sliderbinratio1",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
        updateSliderInput(
          session,
          "sliderbinratio2",
          min = 0,
          max = LIST_DATA$x_plot_range[2],
          value = c(0,0)
        )
      }
    }
    
    if (input$tabs == "clustertool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectclusterfile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(session,
                        "selectclusterfile",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol)
      if (sum(grepl("Cluster_1\nn =", names(LIST_DATA$gene_file))) == 0) {
        updateSliderInput(
          session,
          "sliderbincluster",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
      }
    }
    
    toggle(
      "selectlineslablesshow",
      condition = (input$tabs == "mainplot" &
                     LIST_DATA$STATE[1] != 0)
    )
    # first time switch tab auto plot
    if (input$tabs == "mainplot" & LIST_DATA$STATE[4] == 0) {
      reactive_values$Apply_Math <-
        ApplyMath(
          LIST_DATA,
          input$myMath,
          input$checkboxrgf,
          input$checkboxrf,
          input$numericnormbin
        )
      if (!is.null(reactive_values$Apply_Math)) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        enable("showmainplot")
        enable("selectlineslablesshow")
        LIST_DATA$STATE[c(1, 4)] <<- 1
      } else{
        disable("showmainplot")
        disable("selectlineslablesshow")
      }
    } else {
      toggle("actionmyplotshow",
             condition = (input$tabs == "mainplot" &
                            LIST_DATA$STATE[1] == 2))
    }
    
  })
  
  # loads data file(s) ----
  observeEvent(input$filetable, {
    print("load file")
    # add warnings for total size of LIST_DATA TODO
    disable("startoff")
    disable("showmainplot")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     LoadTableFile(input$filetable$datapath,
                                   input$filetable$name,
                                   LIST_DATA)
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
    } else {
      return()
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    if (LIST_DATA$STATE[1] == 0) {
      show("filegene1")
      show("checkboxconvert")
      show("downloadGeneList")
      show("checkboxsavesplit")
      show("filecolor")
      show("showmainplot")
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
      num_bins <-
        collapse(summarise(LIST_DATA$table_file[[1]], max(bin)))[[1]]
      if (num_bins == 80) {
      } else if (num_bins <= 30) {
        updateSelectInput(session, "selectlineslables", selected = "543 bins 10,10,10")
      }
      
      
      LIST_DATA$STATE[1] <<- 1
    }
    if (LIST_DATA$STATE[4] != 0) {
      LIST_DATA$STATE[4] <<- 3
    }
    enable("startoff")
    enable("showmainplot")
    reset("filetable")
    ff <- names(LIST_DATA$table_file)
    updateAwesomeRadio(session,
                       "radiodataoption",
                       choices = ff,
                       selected = last(ff))
  })
  
  # loads gene list file ----
  observeEvent(input$filegene1, {
    print("load gene file")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   # load info, update select boxes, switching works and chaning info and ploting
                   LD <- LoadTableFile(
                     input$filegene1$datapath,
                     input$filegene1$name,
                     LIST_DATA,
                     TRUE,
                     input$checkboxconvert
                   )
                 })
    if (!is.null(LD)) {
      LIST_DATA <<- LD
    }
    if (LIST_DATA$STATE[4] != 0) {
      LIST_DATA$STATE[4] <<- 3
    }
    reset("filegene1")
    updateCheckboxInput(session, "checkboxconvert", value = FALSE)
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
  })
  
  # loads color file ----
  observeEvent(input$filecolor, {
    my_sel <- input$radiodataoption
    my_list <- input$selectgenelistoptions
    print("load color file")
    # load info, update select boxes, switching works and chaning info and ploting
    LIST_DATA <<- LoadColorFile(input$filecolor$datapath,
                                LIST_DATA, my_list)
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
  
  # update desplay selected item info ----
  observeEvent(c(input$radiodataoption, input$selectgenelistoptions),
               ignoreInit = TRUE,
               {
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
                 if (my_list == names(LIST_DATA$gene_info)[1]) {
                   show("normfactor")
                   hide("actionremovegene")
                 } else {
                   hide("normfactor")
                   show("actionremovegene")
                 }
               })
  
  # saves gene list ----
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      if (input$loadfiles == "Color") {
        paste(Sys.Date(), ".color.txt", sep = "")
      } else {
        paste(gsub("\nn = ", " n = ", input$selectgenelistoptions),
              Sys.Date(),
              ".txt",
              sep = "_")
      }
    },
    content = function(file) {
      if (input$loadfiles == "Color") {
        new_comments <- NULL
        for (i in names(LIST_DATA$gene_info[[input$selectgenelistoptions]])) {
          new_comments <-
            c(new_comments, paste(i, LIST_DATA$gene_info[[input$selectgenelistoptions]][[i]][4]))
        }
      } else {
        new_comments <- paste("#", Sys.Date(), "\n# File(s) used:")
        new_comments <-
          c(new_comments, paste("#", names(LIST_DATA$table_file)))
        new_comments <-
          c(new_comments,  paste("\n#", gsub(
            "\nn = ", " n = ",  input$selectgenelistoptions
          )))
        new_comments <-
          c(new_comments, paste("#", gsub(
            "\nn = ", " n = ",
            paste(LIST_DATA$gene_file[[input$selectgenelistoptions]]$info)
          )))
        if (input$checkboxsavesplit) {
          new_comments <-
            c(new_comments, gsub(
              ";|\\+;|\\-;|\\|",
              "\t",
              pull(LIST_DATA$gene_file[[input$selectgenelistoptions]]$use)
            ))
        } else {
          new_comments <-
            c(new_comments, pull(LIST_DATA$gene_file[[input$selectgenelistoptions]]$use))
        }
        
      }
      write_lines(new_comments, file)
    }
  )
  
  
  # record new nickname and norm factor ----
  observeEvent(input$actionoptions, ignoreInit = TRUE, {
    if (!is.na(input$normfactor) &
        !input$normfactor %in% c(0, 1) &
        input$normfactor != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"]) {
      print("norm")
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"] <<-
        input$normfactor
      LIST_DATA$table_file[[input$radiodataoption]] <<-
        mutate(LIST_DATA$table_file[[input$radiodataoption]],
               score = score / as.numeric(input$normfactor))
      if (LIST_DATA$STATE[1] == 1) {
        reactive_values$Apply_Math <-
          ApplyMath(
            LIST_DATA,
            input$myMath,
            input$checkboxrgf,
            input$checkboxrf,
            input$numericnormbin
          )
      }
    } else if (!is.na(input$normfactor)) {
      updateNumericInput(session,
                         "normfactor",
                         value = as.numeric(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"]))
    }
    
    if (nchar(input$textnickname) > 0) {
      print("new nickname")
      oldnickname <-
        paste(
          gsub("(.{17})", "\\1\n", input$selectgenelistoptions),
          gsub("(.{17})", "\\1\n",
               LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["set"]),
          sep = '\n'
        )
      
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
        newnickname <-
          paste(
            gsub("(.{17})", "\\1\n", input$selectgenelistoptions),
            gsub("(.{17})", "\\1\n", input$textnickname),
            sep = '\n'
          )
        reactive_values$Apply_Math <-
          reactive_values$Apply_Math %>%
          mutate(set = replace(set, set == oldnickname, newnickname))
        
        reactive_values$Plot_Options <-
          reactive_values$Plot_Options %>%
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
  observeEvent(input$colourhex, ignoreInit = TRUE, {
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
  
  # reactive picker watcher (watches everything) ----
  observeEvent(reactiveValuesToList(input)[gsub("\nn = ", "-bensspace-", names(LIST_DATA$gene_info))],
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 reactive_values$picker <-
                   reactiveValuesToList(input)[gsub("\nn = ", "-bensspace-", names(LIST_DATA$gene_info))]
                 
                 
               })
  
  #records check box on/off ----
  observeEvent(reactive_values$picker,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 if (LIST_DATA$STATE[4] != 0) {
                   if (LIST_DATA$STATE[4] == 3) {
                     LIST_DATA$STATE[4] <<- 1
                     return()
                   }
                   print("checkbox on/off")
                   
                   ttt <- reactive_values$picker
                   
                   checkboxonoff <- list()
                   for (i in names(ttt)) {
                     for (tt in ttt[i]) {
                       selectgenelistonoff <- gsub("-bensspace-", "\nn = ", i)
                       checkboxonoff[[selectgenelistonoff]] <-
                         c(checkboxonoff[[selectgenelistonoff]], tt)
                     }
                     
                   }
                   LIST_DATA$gene_info <<-
                     CheckBoxOnOff(checkboxonoff,
                                   LIST_DATA$gene_info)
                   
                   if (LIST_DATA$STATE[4] == 2) {
                     LIST_DATA$STATE[1] <<- 2
                     print("toggle on/off")
                     toggle("actionmyplotshow",
                            condition = (input$tabs == "mainplot"))
                     # reactive_values$Plot_controler <- plot(0,type='n',axes=FALSE,ann=FALSE)
                     disable("selectlineslablesshow")
                     disable("showmainplot")
                   } else {
                     LIST_DATA$STATE[4] <<- 2
                   }
                 }
                 
               })
  
  # plots when action button is pressed ----
  observeEvent(input$actionmyplot, ignoreInit = TRUE, {
    print("plot button")
    reactive_values$Apply_Math <-
      ApplyMath(
        LIST_DATA,
        input$myMath,
        input$checkboxrgf,
        input$checkboxrf,
        input$numericnormbin
      )
    if (!is.null(reactive_values$Apply_Math)) {
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      enable("showmainplot")
      enable("selectlineslablesshow")
    } else{
      disable("showmainplot")
      disable("selectlineslablesshow")
      text = paste("Nothing selected to plot.\n")
      reactive_values$Plot_controler <- ggplot() +
        annotate(
          "text",
          x = 4,
          y = 25,
          size = 8,
          label = text
        ) +
        theme_void()
      return()
    }
    hide("actionmyplotshow")
    LIST_DATA$STATE[1] <<- 1
  })
  
  # updates y axis limits
  observeEvent(reactive_values$Apply_Math, {
    print("upate y axix on new math")
    test <- reactive_values$Y_Axis_numbers
    reactive_values$Y_Axis_numbers <-
      MyXSetValues(reactive_values$Apply_Math,
                   input$sliderplotBinRange)
    
    updateSliderInput(
      session,
      "sliderplotYRange",
      min = reactive_values$Y_Axis_numbers[3],
      max = reactive_values$Y_Axis_numbers[4],
      value = reactive_values$Y_Axis_numbers[1:2],
      step = ((
        reactive_values$Y_Axis_numbers[4] -
          reactive_values$Y_Axis_numbers[3]
      ) / 20
      )
    )
    # Forces update if y Asis values stay the same
    if (sum(test) == sum(reactive_values$Y_Axis_numbers)) {
      reactive_values$Y_Axis_plot <- input$actionmyplot[1]
    }
    
  })
  
  # renders plot ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  
  # updates applymath ----
  observeEvent(
    c(
      input$myMath,
      input$checkboxrgf,
      input$numericnormbin,
      input$checkboxrf
    ),
    ignoreInit = TRUE,
    {
      if (is.na(input$numericnormbin)) {
        updateNumericInput(session, "numericnormbin", value = 0)
      } else if (input$numericnormbin < 0) {
        updateNumericInput(session, "numericnormbin", value = 0)
      } else if (input$numericnormbin > LIST_DATA$x_plot_range[2]) {
        updateNumericInput(session, "numericnormbin", value = LIST_DATA$x_plot_range[2])
      } else {
        updateNumericInput(session,
                           "numericnormbin",
                           value = round(input$numericnormbin))
      }
      
      if (input$checkboxrgf &
          reactive_values$mynorm == "checkboxrf") {
        updateCheckboxInput(session, "checkboxrf", value = FALSE)
        reactive_values$mynorm <- "checkboxrgf"
      } else if (input$checkboxrf &
                 reactive_values$mynorm == "checkboxrgf") {
        updateCheckboxInput(session, "checkboxrgf", value = FALSE)
        updateNumericInput(session, "numericnormbin", value = 0)
        reactive_values$mynorm <- "checkboxrf"
      } else if (input$checkboxrf &
                 reactive_values$mynorm == "numericnormbin") {
        updateNumericInput(session, "numericnormbin", value = 0)
        reactive_values$mynorm <- "checkboxrf"
      } else if (input$numericnormbin > 0 &
                 reactive_values$mynorm == "checkboxrf") {
        updateCheckboxInput(session, "checkboxrf", value = FALSE)
        reactive_values$mynorm <- "numericnormbin"
      } else if (input$checkboxrgf) {
        reactive_values$mynorm <- "checkboxrgf"
      } else if (input$checkboxrf) {
        reactive_values$mynorm <- "checkboxrf"
      } else if (input$numericnormbin > 0) {
        reactive_values$mynorm <- "numericnormbin"
      } else {
        reactive_values$mynorm <- "none"
      }
      reactive_values$Y_Axis_Lable <-
        YAxisLable(
          input$myMath,
          input$checkboxrf,
          input$checkboxrgf,
          input$numericnormbin,
          input$checkboxsmooth
        )
      if (LIST_DATA$STATE[1] == 1) {
        print("apply math")
        reactive_values$Apply_Math <-
          ApplyMath(
            LIST_DATA,
            input$myMath,
            input$checkboxrgf,
            input$checkboxrf,
            input$numericnormbin
          )
      }
    }
  )
  
  #plots when bin slider or y slider is triggered ----
  observeEvent(
    c(
      reactive_values$Lines_Lables_List,
      input$sliderplotBinRange,
      reactive_values$Y_Axis_plot,
      input$sliderplotYRange
    ),
    ignoreInit = TRUE,
    {
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
    }
  )
  
  # quick lines and lables preset change #TODO finish update ----
  observeEvent(input$selectlineslables, ignoreInit = TRUE, {
    print("quick Lines & Lables")
    myset <- LinesLablesPreSet(input$selectlineslables)
    updateNumericInput(session, "numericbody1", value = myset[1])
    updateNumericInput(session, "numericbody2", value = myset[2])
    updateNumericInput(session, "numerictss", value = myset[3])
    updateNumericInput(session, "numerictes", value = myset[4])
    updateNumericInput(session, "numericbinsize", value = myset[5])
    updateNumericInput(session, "numericlabelspaceing", value = myset[6])
    
    
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
  observeEvent(input$actionlineslabels, ignoreInit = TRUE, {
    print("action lines and lables")
    myset <- c(
      input$numericbody1,
      input$numericbody2,
      input$numerictss,
      input$numerictes,
      input$numericbinsize,
      input$numericlabelspaceing
    )
    myset[is.na(myset)] <- 0
    
    for (i in seq_along(myset)) {
      if (myset[i] < 1) {
        myset[i] <- 0
      } else if (i %in% c(1:4, 6) &
                 myset[i] > LIST_DATA$x_plot_range[2]) {
        myset[i] <- LIST_DATA$x_plot_range[2]
      }
    }
    updateNumericInput(session, "numericbody1", value = myset[1])
    updateNumericInput(session, "numericbody2", value = myset[2])
    updateNumericInput(session, "numerictss", value = myset[3])
    updateNumericInput(session, "numerictes", value = myset[4])
    updateNumericInput(session, "numericbinsize", value = myset[5])
    updateNumericInput(session, "numericlabelspaceing", value = myset[6])
    
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
  observeEvent(input$checkboxsmooth, ignoreInit = TRUE, {
    reactive_values$Y_Axis_Lable <-
      YAxisLable(
        input$myMath,
        input$checkboxrf,
        input$checkboxrgf,
        input$numericnormbin,
        input$checkboxsmooth
      )
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
  observeEvent(input$kbrewer, ignoreInit = TRUE, {
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
          if (i == common_name) {
            LIST_DATA$gene_info[[i]][[j]][4] <<-
              kListColorSet[color_safe]
          } else {
            nn <- which(names(LIST_DATA$gene_info) == i)
            LIST_DATA$gene_info[[i]][[j]][4] <<-
              RgbToHex(my_hex = kListColorSet[color_safe], tint = nn * .15)
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
  
  # sort tool picker control ----
  observeEvent(input$selectsortfile, ignoreInit = TRUE, {
    print("sort picker update")
    updatePickerInput(
      session,
      "pickersortfile",
      choices = names(LIST_DATA$table_file),
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectsortfile]], "[[", 4)
      ), sep = ":"))
    )
  })
  
  # update plot picker ----
  observeEvent(reactive_values$Picker_controler, ignoreNULL = FALSE, ignoreInit = TRUE,{
    print("plot pickers update")
    output$DynamicGenePicker <- renderUI({
      pickerlist <- list()
      for (i in names(LIST_DATA$gene_info)) {
        pickerlist[[i]] <-
          list(
            pickerInput(
              inputId = gsub("\nn = ", "-bensspace-", i),
              label = i,
              width = "99%",
              choices = names(LIST_DATA$table_file),
              selected = names(LIST_DATA$table_file)[c(sapply(LIST_DATA$gene_info[[i]], "[[", 5) != 0)],
              multiple = T,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 0"
              ),
              choicesOpt = list(style = paste("color", c(
                sapply(LIST_DATA$gene_info[[i]], "[[", 4)
              ), sep = ":"))
            )
          )
      }
      pickerlist
    })
  })
  
  # Remove data file ----
  observeEvent(input$actionremovefile, ignoreInit = TRUE,{
    print("remove file")
    LIST_DATA <<- RemoveFile(LIST_DATA, input$radiodataoption) 
    ff <- names(LIST_DATA$table_file)
    updateAwesomeRadio(session,
                       "radiodataoption",
                       choices = ff,
                       selected = last(ff))
    
    reactive_values$Picker_controler <- names(LIST_DATA$table_file)
  })
  
  # Remove gene list ----
  observeEvent(input$actionremovegene, ignoreInit = TRUE,{
    print("remove gene list")
    LIST_DATA <<- RemoveGeneList(LIST_DATA, input$selectgenelistoptions) 
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
  })
  
  # sort tool picker enable/disable ----
  observeEvent(input$pickersortfile, ignoreNULL = FALSE, ignoreInit = TRUE,{
    if (!is.null(input$pickersortfile)) {
      enable("enablemainsort")
    } else {
      disable("enablemainsort")
      disable('hidesorttable')
    }
  })
  
  # sort tool action ----
  observeEvent(input$actionsorttool, ignoreInit = TRUE, {
    print("sort tool")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- SortTop(
                     LIST_DATA,
                     input$selectsortfile,
                     input$pickersortfile,
                     input$slidersortbinrange[1],
                     input$slidersortbinrange[2],
                     input$slidersortpercent,
                     input$selectsorttop
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
    } else {
      return()
    }
    newnames <-
      gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full))
    if(any(grep("Sort", names(LIST_DATA$gene_info)) > 0)){
      dt <- datatable(
        LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full,
        rownames = FALSE,
        colnames = strtrim(newnames, 24),
        class = 'cell-border stripe compact',
        filter = 'top',
        caption = 'rank percent',
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = TRUE,
          columnDefs = list(
            list(className = 'dt-center ', targets = "_all"),
            list(
              targets = 0,
              render = JS(
                "function(data, type, row, meta) {",
                "return type === 'display' && data.length > 24 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 19) + '...</span>' : data;",
                "}"
              )
            )
          )
        )
      ) %>% formatPercentage(names(LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full)[-1])
    } else {
      dt <- datatable(LIST_DATA$gene_file[[1]]$empty, 
                      rownames = FALSE,
                      colnames = strtrim(newnames, 24),
                      options = list(searching = FALSE))
    }
    output$sorttable <- DT::renderDataTable(dt)
    enable('hidesorttable')
    
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectsortfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectsortfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
  })
  
  # sort quick tool action ----
  observeEvent(input$actionsortquick, ignoreInit = TRUE, {
    print("quick sort")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- SortTop(
                     LIST_DATA,
                     input$selectsortfile,
                     input$pickersortfile,
                     input$slidersortbinrange[1],
                     input$slidersortbinrange[2],
                     input$slidersortpercent,
                     "Quick%"
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
    } else {
      return()
    }
    
    newnames <-
      gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full))
    if(any(grep("Sort", names(LIST_DATA$gene_info)) > 0)){
    dt <- datatable(
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full,
      rownames = FALSE,
      colnames = strtrim(newnames, 24),
      class = 'cell-border stripe compact',
      filter = 'top',
      caption = 'rank percent',
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = TRUE,
        autoWidth = TRUE,
        columnDefs = list(
          list(className = 'dt-center ', targets = "_all"),
          list(
            targets = 0,
            render = JS(
              "function(data, type, row, meta) {",
              "return type === 'display' && data.length > 24 ?",
              "'<span title=\"' + data + '\">' + data.substr(0, 19) + '...</span>' : data;",
              "}"
            )
          )
        )
      )
    ) %>% formatPercentage(names(LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full)[-1])
    } else {
      dt <- datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames, 24),
                    options = list(searching = FALSE))
    }
    
    output$sorttable <- DT::renderDataTable(dt)
    enable('hidesorttable')
  })
  
  # sort tool gene list $use ----
  observeEvent(input$sorttable_rows_all, ignoreInit = TRUE, {
    newname <- strtrim(gsub("(.{30})", "\\1... ", paste0("Sort\nn = ", length(input$sorttable_rows_all), "-", gsub("\nn = ", ":",input$selectsortfile))),33)
    if(newname != LIST_DATA$STATE[2]){
    print("sort filter $use")
    oldname <- grep("Sort\nn =", names(LIST_DATA$gene_info))
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$sorttable_rows_all])
    
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectsortfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectsortfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  # ratio tool picker control ----
  observeEvent(input$selectratiofile, ignoreInit = TRUE, {
    print("ratio picker update")
    updatePickerInput(
      session,
      "pickerratio1file",
      choices = names(LIST_DATA$table_file),
      selected = "",
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
      ), sep = ":"))
    )
    updatePickerInput(
      session,
      "pickerratio2file",
      choices = names(LIST_DATA$table_file),
      selected = "",
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
      ), sep = ":"))
    )
  })
  
  # ratio tool picker enable/disable ----
  observeEvent(c(input$pickerratio1file, input$pickerratio2file), ignoreInit = TRUE,
               ignoreNULL = FALSE,
               {
                 if (input$pickerratio1file != "" | input$pickerratio2file != "") {
                   enable("enablemainratio")
                 } else {
                   disable("enablemainratio")
                   disable('hideratiotable')
                 }
               })
  
  # keep ratio number above 0 ----
  observeEvent(input$numericratio, ignoreInit = TRUE, {
    print("ratio num")
    if(is.numeric(input$numericratio)){
      if(input$numericratio < 0){
        updateNumericInput(session, "numericratio", value = 2)
      }
    } else {
      updateNumericInput(session, "numericratio", value = 2)
      }
    })
  
  # keep ration bins from overlaping ----
  observeEvent(c(input$sliderbinratio1, input$sliderbinratio2), ignoreInit = TRUE, {
    print("ration bin")
    if(sum(input$sliderbinratio2)>0){
      if(input$sliderbinratio1[2] > input$sliderbinratio2[1]){
        updateSliderInput(session, "sliderbinratio2", value = input$sliderbinratio1[2])
      }
    }
  })
  
  # ratio tool gene lists $use ----
  observeEvent(input$ratio1table_rows_all, ignoreInit = TRUE, {
    newname <- paste("Ratio_Up_file1\nn =", length(input$ratio1table_rows_all))
    oldname <- grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("ratio1 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio1table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectratiofile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectratiofile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  observeEvent(input$ratio2table_rows_all, ignoreInit = TRUE, {
    newname <- paste("Ratio_Up_file2\nn =", length(input$ratio1table_rows_all))
    oldname <- grep("Ratio_Up_file2\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("ratio2 filter $use")
    
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio2table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectratiofile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectratiofile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  observeEvent(input$ratio3table_rows_all, ignoreInit = TRUE, {
    newname <-  paste("Ratio_No_Diff\nn =", length(input$ratio3table_rows_all))
    oldname <- grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("no ratio filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio3table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectratiofile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectratiofile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  # Ratio tool action ----
  observeEvent(input$actionratiotool, ignoreInit = TRUE, {
    print("ratio tool action")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     CompareRatios(
                       LIST_DATA,
                       input$selectratiofile,
                       input$pickerratio1file,
                       input$pickerratio2file,
                       input$sliderbinratio1[1],
                       input$sliderbinratio1[2],
                       input$sliderbinratio2[1],
                       input$sliderbinratio2[2],
                       input$numericratio
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
    } else {
      return()
    }
    newnames1 <- gsub("(.{20})", "\\1... ", input$pickerratio1file)
    if(any(grep("Ratio_Up_file1\nn =",names(LIST_DATA$gene_info))>0)){
      mytab <- "Up Fold Change file 1"
    output$ratio1table <-
      DT::renderDataTable(
        datatable(
          LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn =",
                                    names(LIST_DATA$gene_info))]]$full,
          rownames = FALSE,
          colnames = strtrim(newnames1, 24),
          class = 'cell-border stripe compact',
          filter = 'top',
          caption = 'Ratio_Up_file1',
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            scrollY = TRUE,
            autoWidth = TRUE,
            columnDefs = list(
              list(className = 'dt-center ', targets = "_all"),
              list(
                targets = 0,
                render = JS(
                  "function(data, type, row, meta) {",
                  "return type === 'display' && data.length > 44 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                  "}"
                )
              )
            )
          )
        )
      )
    } else {
      output$ratio1table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames1, 24),
                    options = list(searching = FALSE)))
    }
    newnames2 <- gsub("(.{20})", "\\1... ", input$pickerratio2file)
    if(any(grep("Ratio_Up_file2\nn =",names(LIST_DATA$gene_info))>0)){
    output$ratio2table <-
      DT::renderDataTable(
        datatable(
          LIST_DATA$gene_file[[grep("Ratio_Up_file2\nn =",
                                    names(LIST_DATA$gene_info))]]$full,
          rownames = FALSE,
          colnames = strtrim(newnames2, 24),
          class = 'cell-border stripe compact',
          filter = 'top',
          caption = 'Ratio_Up_file2',
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            scrollY = TRUE,
            autoWidth = TRUE,
            columnDefs = list(
              list(className = 'dt-center ', targets = "_all"),
              list(
                targets = 0,
                render = JS(
                  "function(data, type, row, meta) {",
                  "return type === 'display' && data.length > 44 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                  "}"
                )
              )
            )
          )
        )
      )
    } else {
      output$ratio2table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames2, 24),
                    options = list(searching = FALSE)))
    }
    if(any(grep("Ratio_No_Diff\nn =",names(LIST_DATA$gene_info))>0)){
      newnames3 <- gsub("\n", " ", grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_info), value = T))
    output$ratio3table <-
      DT::renderDataTable(
        datatable(
          LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn =",
                                    names(LIST_DATA$gene_info))]]$full,
          rownames = FALSE,
          colnames = newnames3,
          class = 'cell-border stripe compact',
          filter = 'top',
          caption = 'Ratio_No_Diff',
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            scrollY = TRUE,
            autoWidth = TRUE,
            columnDefs = list(
              list(className = 'dt-center ', targets = "_all"),
              list(
                targets = 0,
                render = JS(
                  "function(data, type, row, meta) {",
                  "return type === 'display' && data.length > 44 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                  "}"
                )
              )
            )
          )
        )
      )
  } else {
    output$ratio3table <-
      DT::renderDataTable(
        datatable(LIST_DATA$gene_file[[1]]$empty, 
                  rownames = FALSE,
                  colnames = "Ratio_No_Diff n = 0",
                  options = list(searching = FALSE)))
  }
    enable('ratiotooltab')
  })
  
  # cluster tool picker control ----
  observeEvent(input$selectclusterfile, ignoreInit = TRUE, {
    print("cluster picker update")
    updatePickerInput(
      session,
      "pickerclusterfile",
      choices = names(LIST_DATA$table_file),
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
      ), sep = ":"))
    )
  })
  
  # cluster tool picker enable/disable ----
  observeEvent(input$pickerclusterfile, ignoreInit = TRUE,
               ignoreNULL = FALSE,
               {
                 if (!is.null(input$pickerclusterfile)) {
                   enable("enablemaincluster")
                 } else {
                   disable("enablemaincluster")
                   disable('hideclustertable')
                 }
               })
  
  # cluster tool gene lists $use ----
  observeEvent(input$cluster1table_rows_all, ignoreInit = TRUE, {
    newname <-  paste("Cluster_1\nn =", length(input$cluster1table_rows_all))
    oldname <- grep("Cluster_1\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("cluster1 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster1table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  observeEvent(input$cluster2table_rows_all, ignoreInit = TRUE, {
    newname <-  paste("Cluster_2\nn =", length(input$cluster2table_rows_all))
    oldname <- grep("Cluster_2\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("cluster2 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster2table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  observeEvent(input$cluster3table_rows_all, ignoreInit = TRUE, {
    newname <-  paste("Cluster_3\nn =", length(input$cluster3table_rows_all))
    oldname <- grep("Cluster_3\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("cluster3 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster3table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  observeEvent(input$cluster4table_rows_all, ignoreInit = TRUE, {
    newname <-  paste("Cluster_4\nn =", length(input$cluster4table_rows_all))
    oldname <- grep("Cluster_4\nn =", names(LIST_DATA$gene_info))
    if(newname != oldname){
    print("cluster4 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster4table_rows_all])
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    reactive_values$Picker_controler <- names(LIST_DATA$gene_info)
    }
  })
  
  # Cluster tool action ----
  observeEvent(input$actionclustertool, ignoreInit = TRUE, {
    print("cluster tool action")
    # withProgress(message = 'Calculation in progress',
    #              detail = 'This may take a while...',
    #              value = 0,
    #              {
    #                LD <-
    #                  FindClusters(
    #                    LIST_DATA,
    #                    input$selectclusterfile,
    #                    input$pickerclusterfile,
    #                    input$sliderbincluster[1],
    #                    input$sliderbincluster[2],
    #                    input$selectclusternumber
    #                  )
    #              })
    # if (!is_empty(LD$table_file)) {
    #   LIST_DATA <<- LD
    # } else {
    #   return()
    # }
    newnames1 <- gsub("(.{20})", "\\1... ", input$pickercluster1file)
    if(any(grep("cluster_Up_file1\nn =",names(LIST_DATA$gene_info))>0)){
    output$cluster1table <-
      DT::renderDataTable(
        datatable(
          LIST_DATA$gene_file[[grep("cluster_Up_file1\nn =",
                                    names(LIST_DATA$gene_info))]]$full,
          rownames = FALSE,
          colnames = strtrim(newnames1, 24),
          class = 'cell-border stripe compact',
          filter = 'top',
          caption = 'cluster_Up_file1',
          options = list(
            pageLength = 15,
            scrollX = TRUE,
            scrollY = TRUE,
            autoWidth = TRUE,
            columnDefs = list(
              list(className = 'dt-center ', targets = "_all"),
              list(
                targets = 0,
                render = JS(
                  "function(data, type, row, meta) {",
                  "return type === 'display' && data.length > 44 ?",
                  "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                  "}"
                )
              )
            )
          )
        )
      )
    } else {
      output$cluster1table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames1, 24),
                    options = list(searching = FALSE)))
    }
    
    newnames2 <- gsub("(.{20})", "\\1... ", input$pickercluster2file)
    if(any(grep("cluster_Up_file2\nn =",names(LIST_DATA$gene_info))>0)){
      output$cluster2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep("cluster_Up_file2\nn =",
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames2, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = 'cluster_Up_file2',
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              scrollY = TRUE,
              autoWidth = TRUE,
              columnDefs = list(
                list(className = 'dt-center ', targets = "_all"),
                list(
                  targets = 0,
                  render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 44 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                    "}"
                  )
                )
              )
            )
          )
        )
    } else {
      output$cluster2table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames2, 24),
                    options = list(searching = FALSE)))
    }
    
    newnames3 <- gsub("(.{20})", "\\1... ", input$pickercluster3file)
    if(any(grep("cluster_Up_file3\nn =",names(LIST_DATA$gene_info))>0)){
      output$cluster3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep("cluster_Up_file3\nn =",
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames3, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = 'cluster_Up_file3',
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              scrollY = TRUE,
              autoWidth = TRUE,
              columnDefs = list(
                list(className = 'dt-center ', targets = "_all"),
                list(
                  targets = 0,
                  render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 44 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                    "}"
                  )
                )
              )
            )
          )
        )
    } else {
      output$cluster3table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames3, 24),
                    options = list(searching = FALSE)))
    }
    
    newnames4 <- gsub("(.{20})", "\\1... ", input$pickercluster4file)
    if(any(grep("cluster_Up_file4\nn =",names(LIST_DATA$gene_info))>0)){
      output$cluster4table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep("cluster_Up_file4\nn =",
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames4, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = 'cluster_Up_file4',
            options = list(
              pageLength = 15,
              scrollX = TRUE,
              scrollY = TRUE,
              autoWidth = TRUE,
              columnDefs = list(
                list(className = 'dt-center ', targets = "_all"),
                list(
                  targets = 0,
                  render = JS(
                    "function(data, type, row, meta) {",
                    "return type === 'display' && data.length > 44 ?",
                    "'<span title=\"' + data + '\">' + data.substr(0, 39) + '...</span>' : data;",
                    "}"
                  )
                )
              )
            )
          )
        )
    } else {
      output$cluster4table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames4, 24),
                    options = list(searching = FALSE)))
    }
    enable('hideclustertable')
  })
  
  # hides sidebar on start up
  shinyjs::addClass(selector = "body", class = "sidebar-collapse")
}

# UI -----
ui <- dashboardPage(
  dashboardHeader(title = "Ben Tools"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Load Data", tabName = "loaddata", icon = icon("file")),
      menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
      hidden(div(
        style = "padding-left: 15%;",
        id = "showpicker",
        uiOutput("DynamicGenePicker")
      )),
      menuItem("Sort Tool", tabName = "sorttool", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showsorttoolpicker",
          selectInput(
            inputId = "selectsortfile",
            label = "Select gene list to sort on",
            choices = "Load data file",
            width = "99%"
          ),
          pickerInput(
            inputId = "pickersortfile",
            width = "99%",
            label = "Select files to sort on",
            choices = "Load data file",
            multiple = T,
            options = list(`actions-box` = TRUE, `selected-text-format` = "count > 1")
          )
        )
      ),
      menuItem("Ratio Tool", tabName = "ratiotool", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showratiotoolpicker",
          selectInput(
            inputId = "selectratiofile",
            label = "Select gene list to sort on",
            choices = "Load data file",
            width = "99%"
          ),
          pickerInput(
            inputId = "pickerratio1file",
            width = "99%",
            label = "Select first file",
            choices = "Load data file",
            multiple = F,
            options = list(title = "Select first file")
          ),
          pickerInput(
            inputId = "pickerratio2file",
            width = "99%",
            label = "Select second file",
            choices = "Load data file",
            multiple = F,
            options = list(title = "Select second file")
          )
        )
      ),
      
      menuItem("Cluster Tools", tabName = "clustertool", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showclustertoolpicker",
          selectInput(
            inputId = "selectclusterfile",
            label = "Select gene list to sort on",
            choices = "Load data file",
            width = "99%"
          ),
          pickerInput(
            inputId = "pickerclusterfile",
            width = "99%",
            label = "Select first file",
            choices = "Load data file",
            multiple = F,
            options = list(title = "Select first file")
          )
        )
      ),
      
      hr(style = "color: #fff; background-color: #337ab7; border-color: #2e6da4;"),
      hidden(
        div(
          id = "selectlineslablesshow",
          style = "padding-left: 10%;",
          selectInput(
            "selectlineslables",
            width = "85%",
            label = "quick set lines and lables",
            choices = c(
              "543 bins 20,20,40",
              "543 bins 10,10,10",
              "5' 1k 1k 80bins",
              "5' .25k 10k 205bins",
              "3'",
              "4"
            )
          )
        )
      )
    )
  ),
  dashboardBody(useShinyjs(),
                tabItems(
                  # load data tab
                  tabItem(tabName = "loaddata",
                          fluidRow(
                            tags$style(
                              ".nav-tabs-custom .nav-tabs li.active a { background-color: transparent; border-color: #2e6da4; } "
                            ),
                            tabBox(
                              title = "Load files",
                              id = "loadfiles",
                              width = 4,
                              tabPanel("Table",
                                       div(
                                         style = "height: 200px;",
                                         fluidRow(div(
                                           style = "padding: 5px 15px;",
                                           fileInput(
                                             "filetable",
                                             label = "Load table file",
                                             accept = c('.table'),
                                             multiple = TRUE
                                           )
                                         )),
                                         helpText("load windowed bedGraph file(s)")
                                       )),
                              tabPanel(title = "Gene",
                                       div(style = "height: 200px;",
                                           fluidRow(
                                             div(style = "padding: 5px 15px;",
                                                 hidden(
                                                   fileInput("filegene1",
                                                             label = "Load 1st gene list",
                                                             accept = c('.txt')),
                                                   checkboxInput("checkboxconvert",
                                                                 "gene list partial matching", value = FALSE)
                                                   
                                                 ))
                                           ))),
                              
                              tabPanel(title = "Color",
                                       div(style = "height: 200px;",
                                           fluidRow(
                                             div(
                                               style = "padding: 5px 15px;",
                                               hidden(fileInput(
                                                 "filecolor",
                                                 label = "Load color list",
                                                 accept = c('.color.txt')
                                               )),
                                               helpText("Applies to current gene list")
                                             )
                                           )))
                            ),
                            
                            hidden(div(
                              id = "startoff",
                              box(
                                title =  "Select Gene list",
                                width = 8,
                                status = "primary",
                                solidHeader = T,
                                selectInput("selectgenelistoptions", "", choices = "common"),
                                div(style = "height: 170px; overflow: auto;", fluidRow(
                                  div(
                                    style = "padding: 0px 30px;",
                                    awesomeRadio("radiodataoption", "select data", choices = "Load data file")
                                  )
                                ))
                              ),
                              box(
                                title =  "Set Plot Color Options",
                                width = 4,
                                status = "primary",
                                solidHeader = T,
                                selectInput("selectdot", "Select dot type", choices = kDotOptions),
                                selectInput("selectline", "Select line type", choices = kLineOptions),
                                
                                fluidRow(
                                  # column(1,tags$br()),
                                  box(
                                    width = 12,
                                    status = "info",
                                    background = "light-blue",
                                    colourInput("colourhex", "Select color HEX"),
                                    tags$hr(),
                                    textInput("textrgbtohex", "RGB"),
                                    actionButton("actionmyrgb", "Update HEX color")
                                  )
                                )
                              ),
                              box(
                                title =  "Set File Plot Options",
                                width = 4,
                                status = "primary",
                                solidHeader = T,
                                textInput("textnickname", "Update Nickname"),
                                numericInput("normfactor", "Set norm factor, score/rpm", value = 1),
                                actionButton("actionoptions", "Set Nickname"),
                                helpText("Need to update Nickname and/or nrom factor"),
                                actionButton("actionremovefile", "Remove File")
                              ),
                              box(
                                title = "Save/Remove gene list",
                                width = 4,
                                status = "primary",
                                solidHeader = T,
                                downloadButton("downloadGeneList", "Save Gene/Color List"),
                                checkboxInput("checkboxsavesplit", "split location and name"),
                                helpText("Switch to Color tab to save Color list"),
                                actionButton("actionremovegene", "Remove Gene list")
                              )
                            ))
                          )),
                  
                  # main plot tab
                  tabItem(tabName = "mainplot",
                          fluidRow(box(
                            width = 12, withSpinner(plotOutput("plot"), type = 4),
                            hidden(
                              div(
                                id = "actionmyplotshow",
                                style = "position: absolute; z-index: 1; left: 45%; top: 50%;",
                                actionButton(
                                  "actionmyplot",
                                  "Update Plot",
                                  icon = icon("area-chart"),
                                  style = "color: #fff; background-color: #337ab7; border-color: #2e6da4;"
                                )
                              )
                            )
                          )),
                          hidden(div(
                            id = "showmainplot",  fluidRow(
                              box(
                                title = "Sliders",
                                status = "primary",
                                solidHeader = T,
                                width = 7,
                                collapsible = TRUE,
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
                                style = 'padding:2px;',
                                title = "math",
                                collapsed = F,
                                collapsible = T,
                                width = 2,
                                status = "primary",
                                solidHeader = T,
                                awesomeRadio(
                                  "myMath",
                                  label =
                                    " ",
                                  choices = c("mean", "sum", "median", "var"),
                                  selected = "mean"
                                )
                              ),
                              box(
                                title = "Normalization",
                                status = "primary",
                                solidHeader = T,
                                width = 3,
                                collapsible = TRUE,
                                checkboxInput("checkboxrf", label = "relative frequency"),
                                checkboxInput("checkboxrgf", label = "relative gene frequency"),
                                checkboxInput("checkboxsmooth", label = "smooth"),
                                numericInput("numericnormbin", "Norm to bin", value = 0)
                              ),
                              box(
                                title = "Lines and Labels",
                                width = 9,
                                status = "primary",
                                solidHeader = T,
                                collapsible = TRUE,
                                collapsed = TRUE,
                                div(
                                  style = "padding:2px; display:inline-block;",
                                  numericInput(
                                    "numerictss",
                                    "TSS bin",
                                    value = 15,
                                    min = 0,
                                    max = 100
                                  )
                                ),
                                div(
                                  style = "padding:2px; display:inline-block;",
                                  numericInput(
                                    "numerictes",
                                    "TES bin",
                                    value = 45,
                                    min = 0,
                                    max = 100
                                  )
                                ),
                                div(
                                  style = "padding:2px; display:inline-block;",
                                  numericInput(
                                    "numericbinsize",
                                    "bp/bin",
                                    value = 100,
                                    min = 20,
                                    max = 1000,
                                    step = 5
                                  )
                                ),
                                div(
                                  style = "padding:2px; display:inline-block;",
                                  numericInput(
                                    "numericbody1",
                                    "5|4 bin",
                                    value = 20,
                                    min = 0,
                                    max = 100
                                  )
                                ),
                                div(
                                  style = "padding:2px; display:inline-block;",
                                  numericInput(
                                    "numericbody2",
                                    "4|3 bin",
                                    value = 40,
                                    min = 0,
                                    max = 100
                                  )
                                ),
                                div(
                                  style = "padding:2px; display:inline-block;",
                                  numericInput(
                                    "numericlabelspaceing",
                                    "every bin",
                                    value = 5,
                                    min = 0,
                                    max = 100
                                  )
                                ),
                                div(
                                  style = "padding-left:33%;",
                                  actionButton("actionlineslabels", "Update Lines and Lables")
                                )
                              ),
                              box(
                                title = "Quick Color Change",
                                width = 3,
                                status = "primary",
                                solidHeader = T,
                                selectInput(
                                  "kbrewer",
                                  "color brewer",
                                  choices = kBrewerList,
                                  selected = kBrewerList[6]
                                )
                              )
                            )
                          ))),
                  # main sort tab
                  tabItem(tabName = "sorttool",
                          div(
                            id = "enablemainsort",
                            box(
                              title = "Sort tool",
                              status = "primary",
                              solidHeader = T,
                              width = 12,
                              fluidRow(
                                column(
                                  2,
                                  selectInput(
                                    "selectsorttop",
                                    "Sort Options",
                                    choices = c("Top%", "Bottom%"),
                                    selected = "Top%"
                                  )
                                ),
                                column(
                                  5,
                                  sliderInput(
                                    "slidersortpercent",
                                    label = "% select:",
                                    post = "%",
                                    min = 1,
                                    max = 100,
                                    value = 80
                                  )
                                ),
                                column(
                                  5,
                                  sliderInput(
                                    "slidersortbinrange",
                                    label = "Select Bin Range:",
                                    min = 0,
                                    max = 80,
                                    value = c(0, 80)
                                  )
                                )
                              ),
                              actionButton("actionsorttool", "Sort"),
                              actionButton("actionsortquick", "Quick Sort")
                            ),
                            div(
                              id = "hidesorttable",
                              box(
                                title = "Sort Table",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                DT::dataTableOutput('sorttable')
                              )
                            )
                          )
                          ),
                  # main ratio tab
                  tabItem(tabName = "ratiotool",
                          div(
                            id = "enablemainratio",
                            box(
                              title = "Ratio tool",
                              status = "primary",
                              solidHeader = T,
                              width = 12,
                              fluidRow(
                                column(
                                  2,
                                  numericInput(
                                    "numericratio",
                                    "Fold Change",
                                    value = 2,
                                    min = 0,
                                    max = 10,
                                    step = 0.1
                                  )
                                ),
                                column(
                                  5,
                                  sliderInput(
                                    "sliderbinratio1",
                                    label = "Select 5' Bin Range:",
                                    min = 0,
                                    max = 80,
                                    value = c(0, 80)
                                  )
                                ),
                                column(
                                  5,
                                  sliderInput(
                                    "sliderbinratio2",
                                    label = "Optinal 3' Bin Range:",
                                    min = 0,
                                    max = 80,
                                    value = c(0, 0)
                                  )
                                )
                              ),
                              actionButton("actionratiotool", "Get fold changes")
                              
                            ),
                            div(
                              id = "hideratiotable",
                              box(
                                title = "Ratio Tables",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                tabBox(
                                  id = "ratiotooltab",
                                  width = 12,
                                  tabPanel("Up Fold Change file 1",
                                           DT::dataTableOutput('ratio1table')),
                                  tabPanel("Up Fold Change file 2",
                                           DT::dataTableOutput('ratio2table')),
                                  tabPanel("No Fold Change",
                                           DT::dataTableOutput('ratio3table'))
                                )
                              )
                            )
                          )),
                  # main cluster tool tab
                  tabItem(tabName = "clustertool",
                          div(
                            id = "enablemaincluster",
                            box(
                              title = "Cluster tools",
                              status = "primary",
                              solidHeader = T,
                              width = 12,
                              fluidRow(
                                column(
                                  2,
                                  selectInput(
                                    inputId = "selectclusternumber",
                                    label = "Select number of clusters",
                                    choices = c(1:4),
                                    selected = 4,
                                    width = "99%"
                                  )
                                ),
                                column(
                                  5,
                                  sliderInput(
                                    "sliderbincluster",
                                    label = "Select Bin Range:",
                                    min = 0,
                                    max = 80,
                                    value = c(0, 80)
                                  )
                                )
                              ),
                              actionButton("actionclustertool", "Get clusters"),
                              actionButton("actiongroupstool", "Get groups")
                              
                            ),
                            div(
                              id = "hideclustertable",
                              box(
                                title = "Cluster Tables",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                tabBox(
                                  id = "clustertooltab",
                                  width = 12,
                                  tabPanel("Cluster 1",
                                           DT::dataTableOutput('cluster1table')),
                                  tabPanel("Cluster 2",
                                           DT::dataTableOutput('cluster2table')),
                                  tabPanel("Cluster 3",
                                           DT::dataTableOutput('cluster3table')),
                                  tabPanel("Cluster 4",
                                           DT::dataTableOutput('cluster4table'))
                                )
                              )
                            )
                          ))
                ))
)

# exicute ----
shinyApp(ui = ui, server = server)