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
    pickerfile_controler = "",
    Y_Axis_Lable = NULL,
    Y_Axis_numbers = NULL,
    Lines_Lables_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL,
    Picker_controler = NULL,
    Y_Axis_plot = NULL
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
    toggle(
      "showcdftoolpicker",
      condition = (input$tabs == "cdftool" &
                     LIST_DATA$STATE[1] != 0)
    )
    toggle(
      "showgenelistspicker",
      condition = (input$tabs == "genelists" &
                     length(LIST_DATA$gene_file) > 1)
    )
    if(input$tabs == "filenorm" & LIST_DATA$STATE[1] != 0){
      updatePickerInput(
        session,
        "pickernumerator", choices = names(LIST_DATA$table_file),
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[1]], "[[", 4)
        ), sep = ":")))
      updatePickerInput(
        session,
        "pickerdenominator", choices = names(LIST_DATA$table_file),
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[1]], "[[", 4)
        ), sep = ":")))
    }
    if(input$tabs == "genelists" & length(LIST_DATA$gene_file) != 0){
      hide('actiongenelistsdatatable')
      if(length(LIST_DATA$gene_file) > 1){
      og <- input$pickergenelists
      if(!all(NULL %in% names(LIST_DATA$gene_file))){
        og <- NULL
      }
      updatePickerInput(
        session,
        "pickergenelists", choices = names(LIST_DATA$gene_file),
        selected = og,
        choicesOpt = list(style = rep("color:black", length(names(LIST_DATA$gene_file))))
      )
      }
    }
    if (input$tabs == "sorttool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectsortfile
      og <- input$pickersortfile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- names(LIST_DATA$gene_file)[1]
      } else if(!all(og %in% names(LIST_DATA$table_file)) | 
                LIST_DATA$STATE[3] == "Sort\nn"){
        og <- NULL
      }
      updateSelectInput(
        session,
        "selectsortfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      updatePickerInput(
        session,
        "pickersortfile",
        choices = names(LIST_DATA$table_file),
        selected = og,
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[input$selectsortfile]], "[[", 4)
        ), sep = ":"))
      )
      if (sum(grepl("Sort\nn =", names(LIST_DATA$gene_file))) == 0) {
        updateSliderInput(
          session,
          "slidersortbinrange",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
        hide('actionsortdatatable')
      }
    }
    
    if (input$tabs == "ratiotool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectratiofile
      og <- c(input$pickerratio1file, input$pickerratio2file)
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- names(LIST_DATA$gene_file)[1]
      } else if(!all(og %in% names(LIST_DATA$table_file)) | 
                LIST_DATA$STATE[3] == "Ratio_"){
        og <- c("","")
      }
      updateSelectInput(session,
                        "selectratiofile",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol)
      updatePickerInput(
        session,
        "pickerratio1file",
        choices = names(LIST_DATA$table_file),
        selected = og[1],
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
        ), sep = ":"))
      )
      updatePickerInput(
        session,
        "pickerratio2file",
        choices = names(LIST_DATA$table_file),
        selected = og[2],
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
        ), sep = ":"))
      )
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
        hide('actionratiodatatable')
      }
    }
    
    if (input$tabs == "clustertool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectclusterfile
      og <- input$pickerclusterfile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- names(LIST_DATA$gene_file)[1]
      }else if(!all(og %in% names(LIST_DATA$table_file)) |
               LIST_DATA$STATE[3] == "Cluste"){
        og <- NULL
      }
      updateSelectInput(session,
                        "selectclusterfile",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol)
      updatePickerInput(
        session,
        "pickerclusterfile",
        choices = names(LIST_DATA$table_file),
        selected = og,
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[input$selectclusterfile]], "[[", 4)
        ), sep = ":"))
      )
      if (sum(grepl("Cluster_1\nn =", names(LIST_DATA$gene_file))) == 0) {
        output$plotcluster <- renderPlot({
          NULL
        })
        updateSliderInput(
          session,
          "sliderbincluster",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = LIST_DATA$x_plot_range
        )
        hide('actionclusterplot')
        hide('actionclusterdatatable')
      }
    }
    
    if (input$tabs == "cdftool" & LIST_DATA$STATE[1] != 0) {
      ol1 <- input$selectcdffile1
      og1 <- input$pickercdffile1
      if(!ol1 %in% names(LIST_DATA$gene_file)){
        ol1 <- names(LIST_DATA$gene_file)[1]
      } else if(!all(og1 %in% names(LIST_DATA$table_file)) | 
                LIST_DATA$STATE[3] == "CDF\nn"){
        og1 <- ""
      }
      updateSelectInput(session,
                        "selectcdffile1",
                        choices = names(LIST_DATA$gene_file),
                        selected = ol1)
      updatePickerInput(
        session,
        "pickercdffile1",
        choices = names(LIST_DATA$table_file),
        selected = og1,
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[input$selectcdffile1]], "[[", 4)
        ), sep = ":"))
      )
      if (sum(grepl("CDF\nn", names(LIST_DATA$gene_file))) == 0) {
        output$plotcdf <- renderPlot({
          NULL
        })
        hide('plotcdf')
        updateSliderInput(
          session,
          "sliderbincdf1",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = c(LIST_DATA$x_plot_range[1],floor(LIST_DATA$x_plot_range[2]/3))
        )
        updateSliderInput(
          session,
          "sliderbincdf2",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = c(ceiling(LIST_DATA$x_plot_range[2]/3),LIST_DATA$x_plot_range[2])
        )
        
        hide('actioncdfdatatable')
      }
    }
    
    toggle(
      "selectlineslablesshow",
      condition = (input$tabs == "mainplot" &
                     LIST_DATA$STATE[1] != 0)
    )
    # first time switch tab auto plot
    if (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0){
      reactive_values$Picker_controler <- sapply(LIST_DATA, function(i) (names(i)))
      if(LIST_DATA$STATE[4] == 0) {
      reactive_values$Apply_Math <-
        ApplyMath(
          LIST_DATA,
          input$myMath,
          input$radioplotnrom,
          as.numeric(input$sliderplotBinNorm)
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
      updateSliderInput(
        session,
        "sliderplotBinNorm",
        min = 0,
        max = LIST_DATA$x_plot_range[2],
        value = 0
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
                 if(LIST_DATA$STATE[1] == 0){
                   return()
                 }
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
            input$radioplotnrom,
            as.numeric(input$sliderplotBinNorm)
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
  observeEvent(reactiveValuesToList(input)[gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", names(LIST_DATA$gene_info)))],
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 reactive_values$picker <-
                   reactiveValuesToList(input)[gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", names(LIST_DATA$gene_info)))]
                 
                 
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
                       selectgenelistonoff <- gsub("-bensspace2-", " ",gsub("-bensspace1-", "\n", i))
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
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
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
  
  # updates norm applymath ----
  observeEvent(
    c(input$myMath,
      input$sliderplotBinNorm,
      input$radioplotnrom
    ),
    ignoreInit = TRUE,
    {
      
      reactive_values$Y_Axis_Lable <-
        YAxisLable(
          input$myMath,
          input$radioplotnrom,
          as.numeric(input$sliderplotBinNorm),
          input$checkboxsmooth
        )
      if (LIST_DATA$STATE[1] == 1) {
        print("apply math")
        reactive_values$Apply_Math <-
          ApplyMath(
            LIST_DATA,
            input$myMath,
            input$radioplotnrom,
            as.numeric(input$sliderplotBinNorm)
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
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm),
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
  
  # update plot picker ----
  observeEvent(reactive_values$Picker_controler, ignoreNULL = FALSE, ignoreInit = TRUE,{
    print("plot pickers update")
    output$DynamicGenePicker <- renderUI({
      pickerlist <- list()
      for (i in names(LIST_DATA$gene_info)) {
        pickerlist[[i]] <-
          list(
            pickerInput(
              inputId = gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", i)),
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
    if (LIST_DATA$STATE[1] > 0) {
    ff <- names(LIST_DATA$table_file)
    updateAwesomeRadio(session,
                       "radiodataoption",
                       choices = ff,
                       selected = last(ff))
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    } else{
      updateAwesomeRadio(session,
                         "radiodataoption",
                         choices = '')
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = 'Load Data File'
      )
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
    
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
     
  })
  
  # create norm file ----
  observeEvent(input$actionnorm, ignoreInit = TRUE, {
    LIST_DATA <<- MakeNormFile(LIST_DATA, input$pickernumerator, 
                               input$pickerdenominator, input$checkboxnormmean,
                               input$checkboxnormzero)
    updatePickerInput(
      session,
      "pickernumerator", selected = "")
    updatePickerInput(
      session,
      "pickerdenominator", selected = "")
    ff <- names(LIST_DATA$table_file)
    updateAwesomeRadio(session,
                       "radiodataoption",
                       choices = ff,
                       selected = last(ff))
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
  })
  
  # Gene lists picker enable/disable ----
  observeEvent(input$actiongenelists, ignoreNULL = FALSE, ignoreInit = TRUE,{
    if(is.null(input$actiongenelists)){
      hide('genelists1table')
      hide('genelists2table')
      hide('genelists3table')
    } 
  })
  
  
  # Gene action ----
  observeEvent(input$actiongenelists,{
    print("gene lists action")
    hide('actiongenelistsdatatable')
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- IntersectGeneLists(
                     LIST_DATA,
                     input$pickergenelists
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      glo <- input$selectgenelistoptions
      if(!glo %in% names(LIST_DATA$gene_file)){
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo)
      ol <- input$pickergenelists
      if(!any(ol %in% names(LIST_DATA$gene_file))){
        ol <- grep("Gene_List_", names(LIST_DATA$gene_file), value = TRUE)
      }else {
      }
      updateSelectInput(
        session,
        "selectsortfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      show('actiongenelistsdatatable')
    } else {
      return()
    }
  })
  
  # Gene lists show gene list ----
  observeEvent(input$actiongenelistsdatatable, ignoreInit = TRUE,{
    print("generiate gene lists table")
    hide('actiongenelistsdatatable')
    if(any(grep("Gene_List_intersect\nn =",names(LIST_DATA$gene_info))>0)){
      newnames1 <- gsub("\n", " ", grep("Gene_List_intersect\nn =",names(LIST_DATA$gene_info), value = TRUE))
      mytab <- "Intersected Gene Lists"
      output$genelists1table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep("Gene_List_intersect\nn =",
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = newnames1,
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = 'Gene_List_intersect',
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
      show('genelists1table')
    } else {
      output$genelists1table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames1, 24),
                    options = list(searching = FALSE)))
      mytab <- "Inclusive Gene Lists"
    }
    if(any(grep("Gene_List_inclusive\nn =",names(LIST_DATA$gene_info))>0)){
      newnames2 <- gsub("\n", " ", grep("Gene_List_inclusive\nn =",names(LIST_DATA$gene_info), value = TRUE))
      output$genelists2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep("Gene_List_inclusive\nn =",
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = newnames2,
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = 'Gene_List_inclusive',
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
      show('genelists2table')
    } else {
      output$genelists2table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames2, 24),
                    options = list(searching = FALSE)))
      if(mytab == "Inclusive Gene Lists"){
        mytab <- "Exclusive Gene Lists"
      }
    }
    if(any(grep("Gene_List_exclusive\nn =",names(LIST_DATA$gene_info))>0)){
      newnames3 <- gsub("\n", " ", grep("Gene_List_exclusive\nn =", names(LIST_DATA$gene_info), value = T))
      output$genelists3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep("Gene_List_exclusive\nn =",
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = newnames3,
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = 'Gene_List_exclusive',
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
      show('genelists3table')
    } else {
      output$genelists3table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = "Gene_List_exclusive n = 0",
                    options = list(searching = FALSE)))
      if(mytab == "Exclusive Gene Lists"){
        mytab <- "Inclusive Gene Lists"
      }
    }
    updateTabItems(session, "geneliststooltab", mytab)
  })
  
  
  # sort tool picker control ----
  observeEvent(input$selectsortfile, ignoreInit = TRUE, {
    print("sort picker update")
    updatePickerInput(
      session,
      "pickersortfile",
      choices = names(LIST_DATA$table_file),
      selected = reactive_values$pickerfile_controler,
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectsortfile]], "[[", 4)
      ), sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  # sort tool picker enable/disable ----
  observeEvent(input$pickersortfile, ignoreNULL = FALSE, ignoreInit = TRUE,{
    if(is.null(input$pickersortfile)){
      hide('sorttable')
    } else if (input$pickersortfile[1] == "") {
      hide('sorttable')
    }
  })
  
  # sort tool action ----
  observeEvent(input$actionsorttool, {
    print("sort tool")
    hide('actionsortdatatable')
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
      show('actionsortdatatable')
      glo <- input$selectgenelistoptions
      if(!glo %in% names(LIST_DATA$gene_file)){
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo)
      ol <- input$selectsortfile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- grep("Sort\nn", names(LIST_DATA$gene_file), value = TRUE)
        reactive_values$pickerfile_controler <- input$pickersortfile
        }else {
        reactive_values$pickerfile_controler <- ""
      }
      updateSelectInput(
        session,
        "selectsortfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
    } else {
      return()
    }
  })
  
  # sort quick tool action ----
  observeEvent(input$actionsortquick, ignoreInit = TRUE, {
    print("quick sort")
    hide('actionsortdatatable')
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
      show('actionsortdatatable')
    } else {
      return()
    }
  })
  
  # Sort generate gene list ----
  observeEvent(input$actionsortdatatable, ignoreInit = TRUE,{
    if(any(grep("Sort\nn", names(LIST_DATA$gene_info)) > 0)){
      newnames <-
        gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$full))
      dt <- datatable(
        LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$full,
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
      ) %>% formatPercentage(names(LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$full)[-1])
    } else {
      dt <- datatable(LIST_DATA$gene_file[[1]]$empty, 
                      rownames = FALSE,
                      colnames = strtrim(newnames, 24),
                      options = list(searching = FALSE))
    }
    output$sorttable <- DT::renderDataTable(dt)
    hide('actionsortdatatable')
    show('sorttable')
  })
  
  # sort tool gene list $use ----
  observeEvent(input$sorttable_rows_all, ignoreInit = TRUE, {
    newname <- strtrim(gsub("(.{30})", "\\1... ", paste0(sub("([0-9]+)", length(input$sorttable_rows_all), LIST_DATA$STATE[2]))), 33)
    oldname <- grep("Sort\nn =", names(LIST_DATA$gene_file))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("sort filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$sorttable_rows_all])
    print("sort $use picker")
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_file),
      selected = glo
    )
    ol <- input$selectsortfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- input$pickersortfile
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectsortfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  # ratio tool picker control ----
  observeEvent(input$selectratiofile, ignoreInit = TRUE, {
    print("ratio picker update")
    if(reactive_values$pickerfile_controler[1] == ""){
      reactive_values$pickerfile_controler <- c("", "")
    }
    updatePickerInput(
      session,
      "pickerratio1file",
      choices = names(LIST_DATA$table_file),
      selected = reactive_values$pickerfile_controler[1],
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
      ), sep = ":"))
    )
    updatePickerInput(
      session,
      "pickerratio2file",
      choices = names(LIST_DATA$table_file),
      selected = reactive_values$pickerfile_controler[2],
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectratiofile]], "[[", 4)
      ), sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  # ratio tool picker enable/disable ----
  observeEvent(c(input$pickerratio1file, input$pickerratio2file), ignoreInit = TRUE,
               ignoreNULL = FALSE,
               {
                 if (input$pickerratio1file != "" | input$pickerratio2file != "") {
                 } else {
                   hide('ratio1table')
                   hide('ratio2table')
                   hide('ratio3table')
                 }
               })
  
  # ratio tool gene lists $use ----
  observeEvent(input$ratio1table_rows_all, ignoreInit = TRUE, {
    newname <- paste("Ratio_Up_file1\nn =", length(input$ratio1table_rows_all))
    oldname <- grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("ratio1 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio1table_rows_all])
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo
    )
    ol <- input$selectratiofile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- c(input$pickerratio1file, input$pickerratio2file)
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectratiofile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  observeEvent(input$ratio2table_rows_all, ignoreInit = TRUE, {
    newname <- paste("Ratio_Up_file2\nn =", length(input$ratio2table_rows_all))
    oldname <- grep("Ratio_Up_file2\nn =", names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("ratio2 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio2table_rows_all])
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo)
    ol <- input$selectratiofile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- c(input$pickerratio1file, input$pickerratio2file)
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectratiofile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  observeEvent(input$ratio3table_rows_all, ignoreInit = TRUE, {
    newname <-  paste("Ratio_No_Diff\nn =", length(input$ratio3table_rows_all))
    oldname <- grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("no ratio filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio3table_rows_all])
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo)
    ol <- input$selectratiofile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- c(input$pickerratio1file, input$pickerratio2file)
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectratiofile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  # Ratio tool action ----
  observeEvent(input$actionratiotool, ignoreInit = TRUE, {
    print("ratio tool action")
    hide('ratio1table')
    hide('ratio2table')
    hide('ratio3table')
    if(is.numeric(input$numericratio)){
      if(input$numericratio < 0){
        updateNumericInput(session, "numericratio", value = 2)
      }
    } else {
      updateNumericInput(session, "numericratio", value = 2)
    }
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
                       input$numericratio,
                       input$checkboxnodivzero
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      show('actionratiodatatable')
      glo <- input$selectgenelistoptions
      if(!glo %in% names(LIST_DATA$gene_file)){
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo)
      ol <- input$selectratiofile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- grep(strsplit(ol, "\nn = ")[[1]][1], names(LIST_DATA$gene_file), value = TRUE)
        reactive_values$pickerfile_controler <- c(input$pickerratio1file, input$pickerratio2file)
      }else {
        reactive_values$pickerfile_controler <- ""
      }
      updateSelectInput(
        session,
        "selectratiofile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
    } else {
      return()
    }
  })
  
  # Ratio show gene list ----
  observeEvent(input$actionratiodatatable, ignoreInit = TRUE,{
    print("generiate ratio table")
    hide('actionratiodatatable')
    if(any(grep("Ratio_Up_file1\nn =",names(LIST_DATA$gene_info))>0)){
      newnames1 <- gsub("\n", " ", grep("Ratio_Up_file1\nn =",names(LIST_DATA$gene_info), value = TRUE))
      mytab <- "Up Fold Change file 1"
    output$ratio1table <-
      DT::renderDataTable(
        datatable(
          LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn =",
                                    names(LIST_DATA$gene_info))]]$full,
          rownames = FALSE,
          colnames = newnames1,
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
    show('ratio1table')
    } else {
      output$ratio1table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames1, 24),
                    options = list(searching = FALSE)))
      mytab <- "Up Fold Change file 2"
    }
    if(any(grep("Ratio_Up_file2\nn =",names(LIST_DATA$gene_info))>0)){
      newnames2 <- gsub("\n", " ", grep("Ratio_Up_file2\nn =",names(LIST_DATA$gene_info), value = TRUE))
    output$ratio2table <-
      DT::renderDataTable(
        datatable(
          LIST_DATA$gene_file[[grep("Ratio_Up_file2\nn =",
                                    names(LIST_DATA$gene_info))]]$full,
          rownames = FALSE,
          colnames = newnames2,
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
    show('ratio2table')
    } else {
      output$ratio2table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames2, 24),
                    options = list(searching = FALSE)))
      if(mytab == "Up Fold Change file 2"){
        mytab <- "No Fold Change"
      }
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
    show('ratio3table')
  } else {
    output$ratio3table <-
      DT::renderDataTable(
        datatable(LIST_DATA$gene_file[[1]]$empty, 
                  rownames = FALSE,
                  colnames = "Ratio_No_Diff n = 0",
                  options = list(searching = FALSE)))
    if(mytab == "No Fold Change"){
      mytab <- "Up Fold Change file 1"
    }
  }
    updateTabItems(session, "ratiotooltab", mytab)
  })
  
  # cluster tool picker control ----
  observeEvent(input$selectclusterfile, ignoreInit = TRUE, {
    print("cluster picker update")
    updatePickerInput(
      session,
      "pickerclusterfile",
      choices = names(LIST_DATA$table_file),
      selected = reactive_values$pickerfile_controler,
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectclusterfile]], "[[", 4)
      ), sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  # cluster tool picker enable/disable ----
  observeEvent(input$pickerclusterfile, ignoreInit = TRUE,
               ignoreNULL = FALSE,
               {
                 if (input$pickerclusterfile != "") {
                 } else {
                   hide("cluster1table")
                   hide("cluster2table")
                   hide("cluster3table")
                   hide("cluster4table")
                   reactive_values$clustergroups <- NULL
                 }
               })
  
  # cluster tool gene lists $use ----
  observeEvent(input$cluster1table_rows_all, ignoreInit = TRUE, {
    newname <-  paste0(reactive_values$clustergroups, "1\nn = ", length(input$cluster1table_rows_all))
    oldname <- grep(paste0(reactive_values$clustergroups, "1\nn ="), names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("cluster1 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster1table_rows_all])
    LD <- LIST_DATA$gene_info
    
    sapply(names(LD), function(i)
      sapply(names(LD[[i]]), function(j)
        if(i %in% grep(reactive_values$clustergroups, names(LD),value = T) & j == input$pickerclusterfile){
          LD[[i]][[j]][5] <<- input$pickerclusterfile
        } else{
          LD[[i]][[j]][5] <<- 0
        })
    )
    LIST_DATA$gene_info <- LD
    reactive_values$Apply_Cluster_Math <- ApplyMath(
      LIST_DATA,
      input$myMath,
      input$radioplotnrom,
      as.numeric(input$sliderplotBinNorm)
    )
    if (!is.null(reactive_values$Apply_Cluster_Math)) {
      reactive_values$Plot_Cluster_Options <- MakePlotOptionFrame(LIST_DATA)
      Y_Axis_Cluster_numbers <-
        MyXSetValues(reactive_values$Apply_Cluster_Math,
                     input$sliderplotBinRange)
      output$plotcluster <- renderPlot({
        GGplotLineDot(
          reactive_values$Apply_Cluster_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Cluster_Options,
          Y_Axis_Cluster_numbers,
          reactive_values$Lines_Lables_List,
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
      })
    }
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- input$pickerclusterfile
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  observeEvent(input$cluster2table_rows_all, ignoreInit = TRUE, {
    newname <-  paste0(reactive_values$clustergroups, "2\nn = ", length(input$cluster2table_rows_all))
    oldname <- grep(paste0(reactive_values$clustergroups, "2\nn ="), names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("cluster2 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster2table_rows_all])
    LD <- LIST_DATA$gene_info
    sapply(names(LD), function(i)
      sapply(names(LD[[i]]), function(j)
        if(i %in% grep(reactive_values$clustergroups, names(LD),value = T) & j == input$pickerclusterfile){
          LD[[i]][[j]][5] <<- input$pickerclusterfile
        } else{
          LD[[i]][[j]][5] <<- 0
        })
    )
    LIST_DATA$gene_info <- LD
    reactive_values$Apply_Cluster_Math <- ApplyMath(
      LIST_DATA,
      input$myMath,
      input$radioplotnrom,
      as.numeric(input$sliderplotBinNorm)
    )
    if (!is.null(reactive_values$Apply_Cluster_Math)) {
      reactive_values$Plot_Cluster_Options <- MakePlotOptionFrame(LIST_DATA)
      Y_Axis_Cluster_numbers <-
        MyXSetValues(reactive_values$Apply_Cluster_Math,
                     input$sliderplotBinRange)
      output$plotcluster <- renderPlot({
        GGplotLineDot(
          reactive_values$Apply_Cluster_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Cluster_Options,
          Y_Axis_Cluster_numbers,
          reactive_values$Lines_Lables_List,
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
      })
    }
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- input$pickerclusterfile
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  observeEvent(input$cluster3table_rows_all, ignoreInit = TRUE, {
    newname <-  paste0(reactive_values$clustergroups, "3\nn = ", length(input$cluster3table_rows_all))
    oldname <- grep(paste0(reactive_values$clustergroups, "3\nn ="), names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("cluster3 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster3table_rows_all])
    LD <- LIST_DATA$gene_info
    sapply(names(LD), function(i)
      sapply(names(LD[[i]]), function(j)
        if(i %in% grep(reactive_values$clustergroups, names(LD),value = T) & j == input$pickerclusterfile){
          LD[[i]][[j]][5] <<- input$pickerclusterfile
        } else{
          LD[[i]][[j]][5] <<- 0
        })
    )
    LIST_DATA$gene_info <- LD
    reactive_values$Apply_Cluster_Math <- ApplyMath(
      LIST_DATA,
      input$myMath,
      input$radioplotnrom,
      as.numeric(input$sliderplotBinNorm)
    )
    if (!is.null(reactive_values$Apply_Cluster_Math)) {
      reactive_values$Plot_Cluster_Options <- MakePlotOptionFrame(LIST_DATA)
      Y_Axis_Cluster_numbers <-
        MyXSetValues(reactive_values$Apply_Cluster_Math,
                     input$sliderplotBinRange)
      output$plotcluster <- renderPlot({
        GGplotLineDot(
          reactive_values$Apply_Cluster_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Cluster_Options,
          Y_Axis_Cluster_numbers,
          reactive_values$Lines_Lables_List,
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
      })
    }
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- input$pickerclusterfile
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  observeEvent(input$cluster4table_rows_all, ignoreInit = TRUE, {
    newname <-  paste0(reactive_values$clustergroups, "4\nn = ", length(input$cluster4table_rows_all))
    oldname <- grep(paste0(reactive_values$clustergroups, "4\nn ="), names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_file)[oldname]){
    print("cluster4 filter $use")
    LIST_DATA$STATE[2] <<- newname
    names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
    names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
    LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
      tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster4table_rows_all])
    LD <- LIST_DATA$gene_info
    
    sapply(names(LD), function(i)
      sapply(names(LD[[i]]), function(j)
        if(i %in% grep(reactive_values$clustergroups, names(LD),value = T) & j == input$pickerclusterfile){
          LD[[i]][[j]][5] <<- input$pickerclusterfile
        } else{
          LD[[i]][[j]][5] <<- 0
        })
    )
    LIST_DATA$gene_info <- LD
    reactive_values$Apply_Cluster_Math <- ApplyMath(
      LIST_DATA,
      input$myMath,
      input$radioplotnrom,
      as.numeric(input$sliderplotBinNorm)
    )
    if (!is.null(reactive_values$Apply_Cluster_Math)) {
      reactive_values$Plot_Cluster_Options <- MakePlotOptionFrame(LIST_DATA)
      Y_Axis_Cluster_numbers <-
        MyXSetValues(reactive_values$Apply_Cluster_Math,
                     input$sliderplotBinRange)
      output$plotcluster <- renderPlot({
        GGplotLineDot(
          reactive_values$Apply_Cluster_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Cluster_Options,
          Y_Axis_Cluster_numbers,
          reactive_values$Lines_Lables_List,
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
      })
    }
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo
    )
    ol <- input$selectclusterfile
    if(!ol %in% names(LIST_DATA$gene_file)){
      ol <- newname
      reactive_values$pickerfile_controler <- input$pickerclusterfile
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectclusterfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
    }
  })
  
  # Cluster tool action ----
  observeEvent(input$actionclustertool, ignoreInit = TRUE, {
    print("cluster tool action")
    hide('plotcluster')
    hide("cluster1table")
    hide("cluster2table")
    hide("cluster3table")
    hide("cluster4table")
    reactive_values$clustergroups <- NULL
    if(n_distinct(LIST_DATA$gene_file[[input$selectclusterfile]]$use) < 4){
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   LD <-
                     FindClusters(
                       LIST_DATA,
                       input$selectclusterfile,
                       input$pickerclusterfile,
                       input$sliderbincluster[1],
                       input$sliderbincluster[2]
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      reactive_values$clustergroups <- "Cluster_"
    } else {
      return()
    }
  })
  
  # Group tool action ----
  observeEvent(input$actiongroupstool, ignoreInit = TRUE, {
    print("group tool action")
    hide('plotcluster')
    hide("cluster1table")
    hide("cluster2table")
    hide("cluster3table")
    hide("cluster4table")
    reactive_values$clustergroups <- NULL
    if(n_distinct(LIST_DATA$gene_file[[input$selectclusterfile]]$use) < 4){
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   LD <-
                     FindGroups(
                       LIST_DATA,
                       input$selectclusterfile,
                       input$pickerclusterfile,
                       input$sliderbincluster[1],
                       input$sliderbincluster[2],
                       input$selectclusternumber
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      reactive_values$clustergroups <- "Group_"
    } else {
      return()
    }
  })
  
  # Cluster tool numbers ----
  observeEvent(c(input$selectclusternumber, reactive_values$clustergroups), ignoreInit = TRUE, {
    print("cluster tool number")
    if(is.null(reactive_values$clustergroups)){
      return()
    }
    hide('actionclusterdatatable')
    hide('actionclusterplot')
    hide('plotcluster')
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0, {
                   LD <-
                     ClusterNumList(
                       LIST_DATA,
                       input$selectclusterfile,
                       input$pickerclusterfile,
                       input$sliderbincluster[1],
                       input$sliderbincluster[2],
                       input$selectclusternumber,
                       reactive_values$clustergroups
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      show('actionclusterdatatable')
      show('actionclusterplot')
      glo <- input$selectgenelistoptions
      if(!glo %in% names(LIST_DATA$gene_file)){
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo)
      ol <- input$selectclusterfile
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- grep(strsplit(ol, "\nn")[[1]][1], names(LIST_DATA$gene_file), value = TRUE)
        reactive_values$pickerfile_controler <- input$pickerclusterfile
      }else {
        reactive_values$pickerfile_controler <- ""
      }
      updateSelectInput(
        session,
        "selectclusterfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
    } else {
      return()
    }
  })
  
  # Create and Show cluster data table ----
  observeEvent(input$actionclusterdatatable, ignoreInit = TRUE, {
    updateTabItems(session, "clustertooltab", "Cluster 1")
    newnames <- gsub("(.{20})", "\\1... ", input$pickerclusterfile)
    if(any(grep(paste0(reactive_values$clustergroups, "1\nn ="),names(LIST_DATA$gene_info))>0)){
      output$cluster1table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep(paste0(reactive_values$clustergroups, "1\nn ="),
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = paste0(reactive_values$clustergroups, "1"),
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
      show("cluster1table")
    } else {
      output$cluster1table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames, 24),
                    options = list(searching = FALSE)))
    }
    
    if(any(grep(paste0(reactive_values$clustergroups, "2\nn ="),names(LIST_DATA$gene_info))>0)){
      output$cluster2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep(paste0(reactive_values$clustergroups, "2\nn ="),
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = paste0(reactive_values$clustergroups, "2"),
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
      show("cluster2table")
    } else {
      output$cluster2table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames, 24),
                    options = list(searching = FALSE)))
    }
    
    if(any(grep(paste0(reactive_values$clustergroups, "3\nn ="),names(LIST_DATA$gene_info))>0)){
      output$cluster3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep(paste0(reactive_values$clustergroups, "3\nn ="),
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = paste0(reactive_values$clustergroups, "3"),
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
      show("cluster3table")
    } else {
      output$cluster3table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames, 24),
                    options = list(searching = FALSE)))
    }
    
    if(any(grep(paste0(reactive_values$clustergroups, "4\nn ="),names(LIST_DATA$gene_info))>0)){
      output$cluster4table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[grep(paste0(reactive_values$clustergroups, "4\nn ="),
                                      names(LIST_DATA$gene_info))]]$full,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            class = 'cell-border stripe compact',
            filter = 'top',
            caption = paste0(reactive_values$clustergroups, "4"),
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
      show("cluster4table")
    } else {
      output$cluster4table <-
        DT::renderDataTable(
          datatable(LIST_DATA$gene_file[[1]]$empty, 
                    rownames = FALSE,
                    colnames = strtrim(newnames, 24),
                    options = list(searching = FALSE)))
    }
    hide('actionclusterdatatable')
  })
  
  # creat and show cluster plot ----
  observeEvent(input$actionclusterplot,{
    print("cluster plot button")
    show('plotcluster')
    LD <- LIST_DATA$gene_info
    sapply(names(LIST_DATA$gene_info), function(i)
      sapply(names(LIST_DATA$gene_info[[i]]), function(j)
        if(i %in% grep(reactive_values$clustergroups, names(LIST_DATA$gene_info),value = T) & j == input$pickerclusterfile){
          LIST_DATA$gene_info[[i]][[j]][5] <<- input$pickerclusterfile
        } else{
          LIST_DATA$gene_info[[i]][[j]][5] <<- 0
        })
    )
    
    reactive_values$Apply_Cluster_Math <- ApplyMath(
        LIST_DATA,
        input$myMath,
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
      )
    if (!is.null(reactive_values$Apply_Cluster_Math)) {
      reactive_values$Plot_Cluster_Options <- MakePlotOptionFrame(LIST_DATA)
      Y_Axis_Cluster_numbers <-
        MyXSetValues(reactive_values$Apply_Cluster_Math,
                     input$sliderplotBinRange)
      output$plotcluster <- renderPlot({
        GGplotLineDot(
          reactive_values$Apply_Cluster_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Cluster_Options,
          Y_Axis_Cluster_numbers,
          reactive_values$Lines_Lables_List,
          input$checkboxsmooth,
          reactive_values$Y_Axis_Lable
        )
      })
    }
    LIST_DATA$gene_info <<- LD
    hide('actionclusterplot')
  })
  
  # CDF tool picker control ----
  observeEvent(c(input$selectcdffile1), ignoreInit = TRUE, {
    print("cdf picker update")
    updatePickerInput(
      session,
      "pickercdffile1",
      choices = names(LIST_DATA$table_file),
      selected = reactive_values$pickerfile_controler,
      choicesOpt = list(style = paste("color", c(
        sapply(LIST_DATA$gene_info[[input$selectcdffile1]], "[[", 4)
      ), sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  
  # CDF tool picker enable/disable ----
  observeEvent(c(input$pickercdffile1), ignoreInit = TRUE,
               ignoreNULL = FALSE, {
                 if(is.null(input$pickercdffile1)){
                   hide('cdftable')
                   hide('plotcdf')
                 } else if (input$pickercdffile1[1] == "") {
                   hide('cdftable')
                   hide('plotcdf')
                 }
               })
  
  # CDF percent reactive
  observeEvent(input$slidercdfper, ignoreInit = TRUE, {
    oldname <- grep("CDF\nn =", names(LIST_DATA$gene_info))
    if(is_empty(oldname)){
      return()
    }
    gene_count <- n_distinct(LIST_DATA$gene_file[[oldname]]$full$gene)
    num <- c(ceiling(gene_count * input$slidercdfper[1]/100), 
             ceiling(gene_count * input$slidercdfper[2]/100))
    gene_list <- group_by(LIST_DATA$gene_file[[oldname]]$full, gene) %>% 
      filter(all(between(bin, num[1], num[2]))) %>% 
      distinct(gene) %>%
      ungroup()
    newname <- paste("CDF\nn =", n_distinct(gene_list$gene))
    if(newname != names(LIST_DATA$gene_info)[oldname]){
      hide('cdftable')
      show('actioncdfdatatable')
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<- gene_list
    df_options <- semi_join(bind_rows(LIST_DATA$gene_info[[newname]]), 
                            distinct(LIST_DATA$gene_file[[newname]]$full, set), 
                            by = "set") %>% 
      mutate(set = paste(sub("\n", " ", newname), 
                         gsub("(.{17})", "\\1\n", set), sep = '\n'))
    df <- inner_join(LIST_DATA$gene_file[[newname]]$full, 
                     LIST_DATA$gene_file[[newname]]$use, by = "gene") %>% 
      mutate(set = paste(sub("\n", " ", newname), 
                         gsub("(.{17})", "\\1\n", set), sep = '\n')) 
    
    use_header <- pull(distinct(df_options, myheader))
    if(n_groups(group_by(df_options, set)) == 2 & n_distinct(df$gene) > 1){
      tt_name <- pull(distinct(df_options, set))
      tt <- suppressWarnings(ks.test(pull(filter(df, set == tt_name[1]),value),
                                     pull(filter(df, set == tt_name[2]),value)))
      if (tt[[2]] == 0) {
        ttt <- "< 2.2e-16"
      } else {
        ttt <- tt[[2]]
        
        use_header <-
          paste(use_header, paste("  p-value = ", format(ttt, scientific = TRUE)))
      }
    }
    output$plotcdf <- renderPlot({
      GGplotC(df, df_options, use_header)
    })
    show('plotcdf')
    }
  })
  
  # CDF generate gene list ----
  observeEvent(input$actioncdfdatatable, ignoreInit = TRUE,{
    if(any(grep("CDF\nn", names(LIST_DATA$gene_info)) > 0)){
      dt <- datatable(
        LIST_DATA$gene_file[[grep("CDF\nn", names(LIST_DATA$gene_info))]]$use,
        rownames = FALSE,
        class = 'cell-border stripe compact',
        filter = 'top',
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          scrollY = TRUE
        )
      )
    } else {
      dt <- datatable(LIST_DATA$gene_file[[1]]$empty, 
                      rownames = FALSE,
                      options = list(searching = FALSE))
    }
    output$cdftable <- DT::renderDataTable(dt)
    hide('actioncdfdatatable')
    show('cdftable')
  })
  

  # cdf tool gene lists $use ----
  observeEvent(input$cdftable_rows_all, ignoreInit = TRUE, {
    newname <- paste("CDF\nn =", length(input$cdftable_rows_all))
    oldname <- grep("CDF\nn =", names(LIST_DATA$gene_info))
    if(newname != names(LIST_DATA$gene_info)[oldname]){
      print("cdf filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cdftable_rows_all])
      df_options <- semi_join(bind_rows(LIST_DATA$gene_info[[newname]]), 
                              distinct(LIST_DATA$gene_file[[newname]]$full, set), 
                              by = "set") %>% 
        mutate(set = paste(sub("\n", " ", newname), 
                           gsub("(.{17})", "\\1\n", set), sep = '\n'))
      df <- inner_join(LIST_DATA$gene_file[[newname]]$full, 
                       LIST_DATA$gene_file[[newname]]$use, by = "gene") %>% 
        mutate(set = paste(sub("\n", " ", newname), 
                           gsub("(.{17})", "\\1\n", set), sep = '\n')) 
      
      use_header <- pull(distinct(df_options, myheader))
      if(n_groups(group_by(df_options, set)) == 2 & n_distinct(df$gene) > 1){
        tt_name <- pull(distinct(df_options, set))
        tt <- suppressWarnings(ks.test(pull(filter(df, set == tt_name[1]),value),
                                       pull(filter(df, set == tt_name[2]),value)))
        if (tt[[2]] == 0) {
          ttt <- "< 2.2e-16"
        } else {
          ttt <- tt[[2]]
          
          use_header <-
            paste(use_header, paste("  p-value = ", format(ttt, scientific = TRUE)))
        }
      }
      output$plotcdf <- renderPlot({
        GGplotC(df, df_options, use_header)
      })
    glo <- input$selectgenelistoptions
    if(!glo %in% names(LIST_DATA$gene_file)){
      glo <- names(LIST_DATA$gene_file)[1]
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = glo
    )
    ol1 <- input$selectcdffile1
    if(!ol1 %in% names(LIST_DATA$gene_file)){
      ol1 <- newname
      reactive_values$pickerfile_controler <- input$pickercdffile1
    }else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectcdffile1",
      choices = names(LIST_DATA$gene_file),
      selected = ol1
    )
    }
  })

  # CDF tool action ----
  observeEvent(input$actioncdftool, ignoreInit = TRUE, {
    print("CDF tool action")
    hide('cdftable')
    hide('plotcdf')
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <<-
                     CumulativeDistribution(
                       LIST_DATA,
                       input$selectcdffile1,
                       input$pickercdffile1,
                       input$sliderbincdf1[1],
                       input$sliderbincdf1[2],
                       input$sliderbincdf2[1],
                       input$sliderbincdf2[2],
                       input$slidercdfper[1],
                       input$slidercdfper[2],
                       input$checkboxnodivzero
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      show('actioncdfdatatable')
      show('plotcdf')
      newname <- grep("CDF\nn =", names(LIST_DATA$gene_info), value = TRUE)
      df_options <- semi_join(bind_rows(LIST_DATA$gene_info[[newname]]), 
                              distinct(LIST_DATA$gene_file[[newname]]$full, set), 
                              by = "set") %>% 
        mutate(set = paste(sub("\n", " ", newname), 
                           gsub("(.{17})", "\\1\n", set), sep = '\n'))
      df <- inner_join(LIST_DATA$gene_file[[newname]]$full, 
                       LIST_DATA$gene_file[[newname]]$use, by = "gene") %>% 
        mutate(set = paste(sub("\n", " ", newname), 
                           gsub("(.{17})", "\\1\n", set), sep = '\n')) 
      
      use_header <- pull(distinct(df_options, myheader))
      if(n_groups(group_by(df_options, set)) == 2 & n_distinct(df$gene) > 1){
        tt_name <- pull(distinct(df_options, set))
        tt <- suppressWarnings(ks.test(pull(filter(df, set == tt_name[1]),value),
                                       pull(filter(df, set == tt_name[2]),value)))
        if (tt[[2]] == 0) {
          ttt <- "< 2.2e-16"
        } else {
          ttt <- tt[[2]]

          use_header <-
            paste(use_header, paste("  p-value = ", format(ttt, scientific = TRUE)))
        }
      }
      output$plotcdf <- renderPlot({
        GGplotC(df, df_options, use_header)
      })
      glo <- input$selectgenelistoptions
      if(!glo %in% names(LIST_DATA$gene_file)){
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo)
      ol <- input$selectcdffile1
      if(!ol %in% names(LIST_DATA$gene_file)){
        ol <- grep("CDF\nn", names(LIST_DATA$gene_file), value = TRUE)
        reactive_values$pickerfile_controler <- input$pickercdffile1
      }else {
        reactive_values$pickerfile_controler <- ""
      }
      updateSelectInput(
        session,
        "selectcdffile1",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
    } else {
      return()
    }
  })


  # hides sidebar on start up ----
  shinyjs::addClass(selector = "body", class = "sidebar-collapse")
}

# UI -----
ui <- dashboardPage(
  dashboardHeader(title = "Ben Tools"),
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Load Data", tabName = "loaddata", icon = icon("file")),
      menuItem("Norm data", tabName = "filenorm", icon = icon("files-o")),
      menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
      hidden(div(
        style = "padding-left: 15%;",
        id = "showpicker",
        uiOutput("DynamicGenePicker")
      )),
      menuItem("Gene Lists", tabName = "genelists", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showgenelistspicker",
          pickerInput(
            inputId = "pickergenelists",
            width = "99%",
            label = "Select Gene lists",
            choices = "Load data file",
            multiple = T,
            options = list(`actions-box` = TRUE, `selected-text-format` = "count > 1")
          )
        )
      ),
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
      menuItem("CDF Tools", tabName = "cdftool", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showcdftoolpicker",
          selectInput(
            inputId = "selectcdffile1",
            label = "Select gene list",
            choices = "Load data file",
            width = "99%"
          ),
          pickerInput(
            inputId = "pickercdffile1",
            width = "99%",
            label = "Select file(s)",
            choices = "Load data file",
            multiple = T,
            options = list(`actions-box` = TRUE, `selected-text-format` = "count > 1")
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
                  # load data tab ----
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
                              ),
                              box(
                                title = "Quick Color Change",
                                width = 4,
                                status = "primary",
                                solidHeader = T,
                                selectInput(
                                  "kbrewer",
                                  "color brewer",
                                  choices = kBrewerList,
                                  selected = kBrewerList[6]
                                )
                              )
                            ))
                          )),
                  # Norm file tab ----
                  tabItem(tabName = "filenorm",
                          
                          box(title = "Select files for normalization",
                            width = 12, status = "primary",
                               solidHeader = T,
                            div(
                              style = "padding-left: 15%;",
                              fluidRow(
                            pickerInput("pickernumerator", label = "numerator",
                                        width = "90%", choices = "Load data file",
                                        multiple = F,
                                        options = list(title = "Select first file")))),
                            div(
                              style = "padding-left: 15%;",
                              fluidRow(
                              pickerInput("pickerdenominator", label = "denominator",
                                          width = "90%", choices = "Load data file",
                                          multiple = F,
                                          options = list(title = "Select second file")))),
                            actionButton("actionnorm",label = "create norm file"),
                            checkboxInput("checkboxnormmean", label = "gene by gene", value = TRUE),
                            checkboxInput("checkboxnormzero", label = "denom 0 -> min/2", value = TRUE),
                            helpText("if 0's are not converted genes containing will be removed from all gene lists")
                            
                          )),
                  
                  # main plot tab ----
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
                                width = 6,
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
                                title = "Sliders",
                                status = "primary",
                                solidHeader = T,
                                width = 6,
                                collapsible = TRUE,
                                sliderInput(
                                  "sliderplotBinNorm",
                                  label = "Bin Norm:",
                                  min = 0,
                                  max = 80,
                                  value = 0
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
                                style = 'padding:2px;',
                                title = "Normalization",
                                status = "primary",
                                solidHeader = T,
                                width = 3,
                                collapsible = TRUE,
                                awesomeRadio("radioplotnrom", label = "Set Y Normalization", 
                                             choices = c("none", "relative frequency","rel gene frequency"),
                                             selected = "none"),
                                checkboxInput("checkboxsmooth", label = "smooth")
                              ),
                              box(
                                title = "Lines and Labels",
                                width = 12,
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
                              )
                            )
                          ))),
                  # main gene lists tab ----
                  tabItem(tabName = "genelists",
                          div(
                            id = "enablemaingenelists",
                            box(
                              title = "Gene Lists",
                              status = "primary",
                              solidHeader = T,
                              width = 12,
                              actionButton("actiongenelists", "Get fold changes"),
                              helpText("Shows intersected, exlusive, and inclusive gene lists")
                            ),
                              box(
                                title = "Gene List Tables",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                actionButton("actiongenelistsdatatable", "Show gene list"),
                                tabBox(
                                  id = "geneliststooltab",
                                  width = 12,
                                  tabPanel("Intersected Gene Lists",
                                           DT::dataTableOutput('genelists1table')),
                                  tabPanel("Inclusive Gene Lists",
                                           DT::dataTableOutput('genelists2table')),
                                  tabPanel("Exclusive Gene Lists",
                                           DT::dataTableOutput('genelists3table'))
                                )
                              )
                          )),
                  # main sort tab ----
                  tabItem(tabName = "sorttool",
                          div(id = "enablemainsort",
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
                                actionButton("actionsortdatatable", "Show gene list"),
                                DT::dataTableOutput('sorttable')
                              )
                            )
                          )
                          ),
                  # main ratio tab ----
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
                              actionButton("actionratiotool", "Get fold changes"),
                              checkboxInput("checkboxnodivzero", label = "0 to min/2", value = TRUE),
                              helpText("if 0's are not converted gene's containing region with sum(0) will be removed from results")
                            ),
                            div(
                              id = "hideratiotable",
                              box(
                                title = "Ratio Tables",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                actionButton("actionratiodatatable", "Show gene list(s)"),
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
                  # main cluster tool tab ----
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
                                    choices = c(4:2),
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
                            box(title = "Cluster Plot", status = "primary", solidHeader = TRUE,
                                width = 12, collapsible = TRUE, collapsed = TRUE,
                                actionButton("actionclusterplot", "plot"),
                                withSpinner(plotOutput("plotcluster"), type = 4)),
                            div(
                              id = "hideclustertable",
                              box(
                                title = "Cluster Tables",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                actionButton("actionclusterdatatable", "Show gene list(s)"),
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
                          )),
                  # main CDF tool tab ----
                  tabItem(tabName = "cdftool",
                          div(id = "enablemaincdf",
                            box(
                              title = "CDF tool",
                              status = "primary",
                              solidHeader = T,
                              width = 12,
                              fluidRow(
                                column(
                                  5,
                                  sliderInput(
                                    "sliderbincdf1",
                                    label = "Select 5' Bin Range:",
                                    min = 0,
                                    max = 80,
                                    value = c(0, 80)
                                  )
                                ),
                                column(
                                  5,
                                  sliderInput(
                                    "sliderbincdf2",
                                    label = "Select 3' Bin Range:",
                                    min = 0,
                                    max = 80,
                                    value = c(0, 0)
                                  )
                                )
                              ),
                              column(
                                5,
                                sliderInput(
                                  "slidercdfper",
                                  label = "Select upper and lower %s",
                                  post = "%",
                                  min = 0,
                                  max = 100,
                                  value = c(0, 100)
                                )
                              ),
                              actionButton("actioncdftool", "Plot CDF"),
                              checkboxInput("checkboxnodivzerocdf", label = "0 to min/2", value = TRUE),
                              helpText("if 0's are not converted gene's containing region with sum(0) will be removed from results")
                            ),
                            box(title = "CDF Plot", status = "primary", solidHeader = TRUE,
                                width = 12, collapsible = TRUE,
                                withSpinner(plotOutput("plotcdf"), type = 4)),
                            div(
                              id = "hidecdftable",
                              box(
                                title = "CDF Tables",
                                status = "primary",
                                solidHeader = T,
                                width = 12,
                                actionButton("actioncdfdatatable", "Show gene list(s)"),
                                DT::dataTableOutput('cdftable')
                              )
                            )
                          ))
                ))
)

# exicute ----
shinyApp(ui = ui, server = server)