# Created by Benjamin Erickson BBErickson@gmail.com

# load/save data tab
#   load other gene list
#   remove gene list (exclued Main list)
#   reset app
#   save gene lists, options for list, bed, tool
# data options and QC tab
# plot tab
# set observeEvents for all items, icon size/color?
# work on adding other boxes and sizes and functions

# program for loading packages ----
my_packages <- function(x) {
  for (i in x) {
    #  require returns TRUE invisibly if it was able to load package
    if (!require(i , character.only = TRUE)) {
      #  If package was not able to be loaded then re-install
      install.packages(i , dependencies = TRUE,)
      print(paste("installing ", i, " : please wait"))
    }
    #  Load package after installing
    require(i , character.only = TRUE)
  }
}

# run load needed packages using my_packages(x) ----
suppressPackageStartupMessages(my_packages(
  c("tidyverse",
    "shiny",
    "shinydashboard",
    "shinydashboardPlus",
    "shinyWidgets",
    "shinyjs",
    "RColorBrewer")
))

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 500MB. ----
options(shiny.maxRequestSize = 500 * 1024 ^ 2)

LIST_DATA <<- list(
  table_file = NULL,
  # gene bin score set
  gene_file = NULL,
  # holds $Compleat genes from files and $gene file(s)
  gene_info = NULL,
  # for holding meta data gene file(s) [c("gene_list", "set", "color", plot?, legend)]
  ttest = NULL,
  # t.test results $use is for numbers $gene_info for holding plotting options
  clust = NULL,
  # Cluster holder
  x_plot_range = c(0, 0),
  STATE = c(0, 0, 5) # flow control,
  # [1] 1 = at least one file has been loaded and lets reactive fill in info
  #
  # [2] 0 = first time switching tab auto plotting
  #     1 = hidden plot button, reactive for plot enabled
  #     2 = on/off reactive picker changed, shows plot button, reactive for plot disabled
  # [3] line and label type i.e. 5, 4, 3, 543
)

# Brewer color sets to be available ----
kBrewerList <-
  c("Set1","Paired","Dark2","Spectral")

# lines and labels types ----
kLinesandlabels <- c(
  "543 bins 20,20,40",
  "5' 1k 1k 80bins",
  "3' 1k 9k 100bins",
  "5' .25k 10k 205bins",
  "5'",
  "4",
  "3",
  "generic 543"
)

# color Brewer set that is active to use in plot remove all yellows ----
kListColorSet <- brewer.pal(11, kBrewerList[2]) %>% grep("#FF",.,value = T,invert = T)

# math options available ----
kMathOptions <- c("mean", "sum", "median")

# functions ----

# basic test for valid rgb color and type and
# switches hex <--> rgb, able to apply tint colors
RgbToHex <- function(x, 
                     convert = "hex", 
                     tint = FALSE) {
  if(str_count(x, ",") == 2){
    # hex <- rgb
    myhex <- try(rgb(str_split_fixed(x,",",n=3),
                     maxColorValue = 255), silent = TRUE)
    if("try-error" %in% class(myhex)){
      myhex <- "#000000"
      mygrb <- "0,0,0"
    } else {
      myrgb <- x
      if (is.numeric(tint) & between(tint,0,1)) {
        myrgb <- as.numeric(str_split_fixed(myrgb,",",n=3))
        myrgb <-
          paste(round(myrgb + (255 - myrgb) * tint), collapse = ",")
        myhex <- rgb(str_split_fixed(myrgb,",",n=3),
                     maxColorValue = 255)
      } 
    }
  } else {
    # rgb <- hex (tint)
    myrgb <- try(col2rgb(x), silent = TRUE)
    if("try-error" %in% class(myrgb)){
      myhex <- "#000000"
      myrgb <- "0,0,0"
    } else {
      if (is.numeric(tint) & between(tint,0,1)) {
        myrgb <-
          paste(round(as.numeric(myrgb) + (255 - as.numeric(myrgb)) * tint), collapse = ",")
        myhex <- rgb(str_split_fixed(myrgb,",",n=3),
                     maxColorValue = 255)
        
      } else {
        myhex <- x
        myrgb <- paste(myrgb, collapse = ",")
      }
    }
  }
  if(convert == "hex"){
    return(myhex)
  } else {
    return(myrgb) 
  }
}

# reads in table file(s), tests, fills out info and returns list_data
LoadTableFile <-
  function(file_path,
           file_name,
           list_data) {
    my_color <- NULL
    my_landl <- NA
    legend_nickname <- NULL
    my_remote_file <- NULL
    file_count <- length(list_data$table_file)
    # shiny progress bar
    setProgress(1, detail = "start gathering information on file(s)")
    # tests if loading a file with a list of address to remote files, requires .url.txt in file name
    for (i in seq_along(file_name)) {
      if (str_detect(file_name[i], ".url.txt")) {
        num_col <-
          try(count_fields(file_path[i],
                           n_max = 1,
                           skip = 1,
                           tokenizer = tokenizer_delim(" ")),silent = T)
        if ("try-error" %in% class(num_col)) {
          showModal(modalDialog(
            title = "Information message",
            paste(file_name[i], "cant find file to load"),
            size = "s",
            easyClose = TRUE
          ))
          next()
        }
        if(num_col > 1){
          col_names <- c("file","type","nick","color")
          col_types <- cols(file=col_character(),
                           type=col_character(),
                           nick=col_character(),
                           color=col_character())
        } else {
          col_names <- c("file")
          col_types <- cols(file=col_character())
        }
        meta_data <- suppressMessages(read_delim(file_path[i],delim = " ", 
                                                 col_names = col_names,
                                                 col_types = col_types))
        if(num_col == 1){
          meta_data <-  meta_data %>% mutate(type = NA,
                                             nick = NA,
                                             color = NA)
        }
        my_remote_file <- c(my_remote_file, meta_data$file)
        my_color <- c(my_color, meta_data$color) 
        legend_nickname <- c(legend_nickname, str_replace(meta_data$nick,"\\.","_")) 
        my_landl <- last(meta_data$type)
        if(i > 1){
          file_path[i] <- NULL
        } else {
          file_path <- NULL
        }
      } else {
        my_color <- c(my_color, NA)
        legend_nickname <- c(legend_nickname, last(str_split(file_name[i],"/",simplify = T)) %>% 
                               str_remove(., ".table")) %>% str_replace("\\.","_")
      }
    }
    if (!is.null(my_remote_file)) {
      file_path <- c(file_path, my_remote_file)
    }
    # shiny progress bar
    setProgress(2, detail = "getting meta data")
    # loop each item in file_path
    for (x in seq_along(file_path)) {
      # gets number of columns in file, used to guess how to deal with file
      #  and checks if file exits
      num_col <-
        try(count_fields(file_path[x],
                         n_max = 1,
                         skip = 1,
                         tokenizer = tokenizer_tsv()),silent = T)
      if ("try-error" %in% class(num_col)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        next()
      }
      # checking/creating meta data
      if(is.na(legend_nickname)[x]){
        legend_nickname[x] <- last(str_split(file_path[x],"/",simplify = T)) %>% 
          str_remove(., ".table") %>% str_replace("\\.","_")
      }
      if(file_count > 0 ){
        # checks if file with same name is in master list of lists
        if (legend_nickname[x] %in% list_data$gene_info$set) {
          showModal(modalDialog(
            title = "Information message",
            paste(legend_nickname[x], "has already been loaded"),
            size = "s",
            easyClose = TRUE
          ))
          next()
        }
      }
      if(is.na(my_landl)){
        my_landl <- str_extract(legend_nickname[x], "^(\\d)+") %>% 
          replace_na("none") 
      }
      
      if (is.na(my_color[x])) {
        my_color[x] <- sample(suppressWarnings(brewer.pal(11, sample(kBrewerList, size=1))) %>% 
                                grep("#FF",.,value = T,invert = T),size = 1)
      } else {
        my_color[x] <- RgbToHex(my_color[x], convert = "hex")
      }
      
      if (num_col > 6) {
        # guesses is in wide format
        # shiny progress bar
        setProgress(3, detail = "loading wide file and converting")
        tablefile <- suppressMessages(
          read_tsv(
            file_path[x],
            comment = "#",
            col_names = c("gene", 1:(num_col - 1)),
            skip = 1
          ) %>%
            gather(., bin, score, 2:(num_col))
        ) %>%
          dplyr::select(gene, bin, score) %>%
          dplyr::mutate(set = legend_nickname[x],
                        bin = as.numeric(bin),
                        score = as.numeric(score)) %>%
          na_if(Inf) %>%
          replace_na(list(score = 0)) %>% 
          distinct(gene,bin,.keep_all = T)
        # guesses is in long bedtools from (bed or bedGraph)
      } else {
        if (num_col == 4) {
          # settings for new style with meta data info
          col_names <- c("gene", "bin", "score", "set")
          # settings for reading in bedGraph file style
        } else if (num_col == 3) {
          col_names <- c("gene", "bin", "score")
        } else {
          showModal(
            modalDialog(
              title = "I dont know how to load this file",
              "I use binned coverage files: gene bin score optinal(set) ",
              size = "s",
              easyClose = TRUE
            )
          )
          next()
        }
        # shiny progress bar
        setProgress(3, detail = "downloading/reading in file")
        # reads in file
        tablefile <-
          suppressMessages(read_tsv(file_path[x],
                                    comment = "#",
                                    col_names = col_names)) %>%
          dplyr::mutate(set = legend_nickname[x]) %>% na_if(Inf) %>%
          replace_na(list(score = 0)) %>% distinct(gene,bin,.keep_all = T)
      }
      num_bins <- n_distinct(tablefile$bin)
      # shiny progress bar
      setProgress(4, detail = "Checking form problems")
      # checks the number of bins and gene naming skeam is consistent
      if (file_count > 0) {
        if (num_bins != list_data$x_plot_range[2]) {
          showModal(
            modalDialog(
              title = "Information message",
              "Can't load file, different number of bins",
              size = "s",
              easyClose = TRUE
            )
          )
          next()
        }
        # test data is compatible with already loaded data
        gene_names <-
          semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
        if (n_distinct(gene_names$gene) == 0) {
          showModal(
            modalDialog(
              title = "Information message",
              " No genes in common ",
              size = "s",
              easyClose = TRUE
            )
          )
          break()
        } else {
          # make complete gene list
          gene_names <-
            full_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% 
            distinct(gene)
        }
      } else {
        gene_names <- distinct(tablefile, gene)
        list_data$x_plot_range <- c(1, num_bins)
        list_data$STATE[3] <- my_landl
      }
      # sets master gene list name
      my_name <- paste("Compleat\nn =", n_distinct(gene_names$gene))
      if (file_count > 0) {
        list_data$gene_info <- list_data$gene_info %>% 
          dplyr::mutate(gene_list = if_else(str_detect(gene_list,"^Compleat\nn =") , my_name, gene_list))
        names(list_data$gene_file)[1] <- my_name
      }
      # shiny progress bar
      setProgress(5, detail = "building data and adding to tools")
      if (list_data$STATE[2] == 0 &
          n_distinct(list_data$table_file$set) < 3) {
        oo <- legend_nickname[x]
      } else {
        oo <- "0"
      }
      # saves data in list of lists
      list_data$table_file <- distinct(bind_rows(list_data$table_file, tablefile))
      list_data$gene_file[[my_name]]$use <- gene_names
      list_data$gene_file[[my_name]]$info <- paste("full gene list",
                                                   Sys.Date())
      list_data$gene_info <- distinct(bind_rows(list_data$gene_info,tibble(
        gene_list = my_name,
        set = legend_nickname[x],
        mycol = my_color[x],
        onoff = oo,
        sub = " ",
        plot_set = " "
      )))
      
      file_count <- 1
      setProgress(6, detail = "done with this file")
    }
    return(list_data)
  }

# reads in 1 column gene list file(s), tests, fills out info and returns list_data
LoadGeneFile <-
  function(file_path,
           file_name,
           list_data) {
    my_remote_file <- NULL
    legend_nickname <- NULL
    if(length(list_data$table_file) < 1){
      print("needed  or does shiny handle? add message?")
      return(list_data)
    }
    # shiny progress bar
    setProgress(1, detail = "start gathering information on file(s)")
    # tests if loading a file with a list of address to remote files, requirs .url.txt in file name
    for (i in seq_along(file_name)) {
      if (str_detect(file_name[i], ".url.txt")) {
        meta_data <- suppressMessages(read_delim(file_path[i],delim = " ",col_names = FALSE))
        my_remote_file <- c(my_remote_file, meta_data$file)
        if(i > 1){
          file_path[i] <- NULL
        } else {
          file_path <- NULL
        }
      }
    }
    if (!is.null(my_remote_file)) {
      file_path <- c(file_path, my_remote_file)
    }
    # shiny progress bar
    setProgress(2, detail = "downloading/reading in file")
    # loop thourgh each item in file_path
    for (x in seq_along(file_path)) {
      # createing nickname
      legend_nickname[x] <- last(str_split(file_name[x],"/",simplify = T)) %>% 
        str_remove(., ".txt") %>% str_replace("\\.","_")
      # checks if file with same nickname has already been loaded
      if (legend_nickname[x] %in% names(list_data$gene_file)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "has already been loaded"),
          size = "s",
          easyClose = TRUE
        ))
        next()
      }
      # gets number of columns in file, used to guess how to deal with file
      #  and checks if file exits
      num_col <-
        try(count_fields(file_path[x],
                         n_max = 1,
                         skip = 1,
                         tokenizer = tokenizer_tsv()),silent = T)
      if ("try-error" %in% class(num_col)) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "cant find file to load"),
          size = "s",
          easyClose = TRUE
        ))
        next
      }
      if(num_col == 1){
        #normal gene list
        col_names <- c("gene")
      } else {
        showModal(
          modalDialog(
            title = "I dont know how to load this file",
            "I expected a 1 colunm file",
            size = "s",
            easyClose = TRUE
          )
        )
        next()
      }
      # reads in file
      tablefile <-
        suppressMessages(read_tsv(file_path[x],
                                  comment = "#",
                                  col_names = col_names)) %>%
        distinct(gene)
      # shiny progress bar
      setProgress(3, detail = "Checking form problems")
      # checks gene list is a subset of what has been loaded
      gene_names <-
        semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
      # test data is compatible with already loaded data
      if (n_distinct(gene_names$gene) == 0) {
        showModal(
          modalDialog(
            title = "Information message",
            " couldn't find a exact match, checking for partial match ",
            size = "s",
            easyClose = TRUE
          )
        )
        # tries to grep lists and find matches
        # shiny progress bar
        setProgress(4, detail = "looking for gene name matches")
        tablefile <-
          distinct(tibble(gene = str_subset(
            list_data$gene_file[[1]]$use$gene, tablefile$gene
          )))
        if (n_distinct(tablefile$gene) == 0) {
          showModal(
            modalDialog(
              title = "Information message",
              " No genes found after pattern matching search",
              size = "s",
              easyClose = TRUE
            )
          )
          return()
        }
        showModal(
          modalDialog(
            title = "Information message",
            " Don't forget to save the gene list for future use",
            size = "s",
            easyClose = TRUE
          )
        )
      }
      # adds full n count to nickname
      my_name <- paste0(legend_nickname[x], "\nn = ", n_distinct(tablefile$gene))
      # preps meta data
      gene_info <- list_data$gene_info %>% 
        dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
        dplyr::mutate(gene_list = my_name, sub = " ", onoff = "0",
                      plot_set = " ")
      # shiny progress bar
      setProgress(5, detail = "building data and adding list")
      # saves data in list of lists
      list_data$gene_file[[my_name]]$use <- distinct(tablefile, gene)
      list_data$gene_file[[my_name]]$full <- tablefile
      list_data$gene_file[[my_name]]$info <-
        paste("Loaded gene list from file",
              legend_nickname[x],
              Sys.Date())
      list_data$gene_info <- bind_rows(list_data$gene_info, gene_info)
      
      setProgress(6, detail = "done with this file")
    }
    return(list_data)
  }

# server ----
server <- function(input, output, session) {
  output$user <- renderUser({
    dashboardUser(
      name = "BenTools V9.a", 
      image = "https://scontent-den4-1.xx.fbcdn.net/v/t1.6435-1/45791244_10217878289252358_1224105918908596224_n.jpg?_nc_cat=104&ccb=1-3&_nc_sid=7206a8&_nc_ohc=xBluQRu0PkAAX8FYD0X&_nc_ht=scontent-den4-1.xx&oh=d3403edbc08de6d92ca15e6d508ac958&oe=6091B67A", 
      title = "Benjamin Erickson",
      subtitle = "BBErickson@gmail.com"
    )
  })
  # disables tabs on start
  addCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
  addCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
  
  # sidebar observe
  
  # loads data file(s) ----
  observeEvent(input$filetable, {
    print("load file")
    shinyjs::disable("startoff")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
    LD <- LoadTableFile(input$filetable$datapath,
                                   input$filetable$name,
                                   LIST_DATA)
    })

    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
    } else {
      showModal(modalDialog(
        title = "Information message",
        paste("No files loaded"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_file),
      selected = names(LIST_DATA$gene_file)[1]
    )
    ff <- distinct(LIST_DATA$table_file, set)$set
    updateSelectInput(session,
                      "selectdataoption",
                      choices = ff)
    # first time starting
    if (LIST_DATA$STATE[1] == 0) {
      LIST_DATA$STATE[1] <<- 1
      shinyjs::show("startoff")
    }
    # enables tabs after loading file
    shinyjs::enable("startoff")
    shinyjs::reset("filetable")
    removeCssClass(selector = "a[data-value='qcOptions']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='mainplot']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='filenorm']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='genelists']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='sorttool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='ratiotool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='clustertool']", class = "inactiveLink")
    removeCssClass(selector = "a[data-value='cdftool']", class = "inactiveLink")
  })
  
  
  # loads gene list file ----
  observeEvent(input$filegene1, {
    print("load gene file")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   # load info, update select boxes, switching works and chaning info and ploting
                   LD <- LoadGeneFile(
                     input$filegene1$datapath,
                     input$filegene1$name,
                     LIST_DATA
                   )
                 })
    if (!is.null(LD)) {
      LIST_DATA <<- LD
    }
    shinyjs::reset("filegene1")
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_file),
      selected = last(names(LIST_DATA$gene_file))
    )
  })
  
}

# UI -----
ui <- dashboardPage(
  skin = "purple-light",
  options = list(sidebarExpandOnHover = F),
  header = dashboardHeader(userOutput("user")),
  sidebar = dashboardSidebar(id = "sidebar", minified = TRUE, collapsed = TRUE,
                             tags$head(tags$style(".inactiveLink {
                            pointer-events: none;
                           cursor: default;
                           }")), # disables tabs on start
                             sidebarMenu(
                               id = "leftSideTabs",
                               menuItem("Load Data", tabName = "loaddata", icon = icon("file-import")),
                               menuItem("QC/Options", tabName = "qcOptions", icon = icon("clipboard-check")),
                               menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
                               menuItem("Norm data", tabName = "filenorm", icon = icon("files-o")),
                               menuItem("Compare Lists", tabName = "genelists", icon = icon("gears")),
                               menuItem("Filter Tool", tabName = "sorttool", icon = icon("gears")),
                               menuItem("Ratio Tool", tabName = "ratiotool", icon = icon("gears")),
                               menuItem("Cluster Tools", tabName = "clustertool", icon = icon("gears")),
                               menuItem("CDF Tools", tabName = "cdftool", icon = icon("gears"))
                             )
  ),
  body = dashboardBody(
    useShinyjs(),
    tabItems(
      # load data tab ----
      tabItem(tabName = "loaddata",
              box(status = "navy",
                  solidHeader = TRUE,
                  title = "Load .table/URL.txt file",
                  width = 5,
                  style = "height: 200px;",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ),
                  helpText("load windowed bedGraph file(s)")
              ),
              hidden(div(
                id = "startoff",
              box(
                title = "Select Gene list",
                width = 7,
                style = "height: 200px;" ,
                solidHeader = TRUE,
                status = "navy",
                sidebar = boxSidebar(
                  id = "sidebarGenelistColor",
                  icon = icon("palette"),
                  width = 50,
                  pickerInput(inputId = "kbrewer",
                              label = "color brewer theme", 
                              choices = c("select", kBrewerList),
                              selected = "select"
                  ),
                  actionButton("BttnNewColor", "Set color same as Compleat")
                ),
                fileInput("filegene1",
                          label = "Load gene list",
                          accept = c('.txt'),
                fluidRow(align="center",
                         pickerInput("selectgenelistoptions", "", 
                            width = 300, choices = "Compleat"),
                actionButton("actionremovegene", "Remove Gene list"))
              )
              ))
      ),
      tabItem(tabName = "qcOptions",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "mainplot",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "filenorm",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "genelists",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "sorttool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "ratiotool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "clustertool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      ),
      tabItem(tabName = "cdftool",
              box(status = "purple",
                  solidHeader = TRUE,
                  title = "QC Options",
                  fileInput(
                    "filetable",
                    label = "",
                    accept = c('.table'),
                    multiple = TRUE
                  ))
      )
    )
  ),
  controlbar = dashboardControlbar(),
  title = "DashboardPage"
)

# exicute ----
shinyApp(ui = ui, server = server)