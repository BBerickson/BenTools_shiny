# Created by Benjamin Erickson BBErickson@gmail.com

# 
# 
# Go through each tab update 100%
# 1 load table file screen, load table file function with progress bar
# 2 load gene list, function, remove load color, save ... functions 
# 3 brewer list colored, and working function, remove line and dot options (place some where else?)
# 4 remove files and gene lists
# 
# find and fix TODO
# look into shinydashboardplus for better box's

# bonus: filter tool for calling real signal using a mix of peak and percentile/ or machine learning project
# 


print("Bentools V8m Shiny")

# program for loading packages ----
my_packages <- function(x) {
  for (i in x) {
    #  require returns TRUE invisibly if it was able to load package
    if (!require(i , character.only = TRUE)) {
      #  If package was not able to be loaded then re-install
      if(i == "patchwork"){
        devtools::install_github("thomasp85/patchwork")
        }else{
        install.packages(i , dependencies = TRUE)
        print(paste("installing ", i, " : please wait"))
      }
      #  Load package after installing
      require(i , character.only = TRUE)
    }
  }
}

# run load needed packages using my_packages(x) ----
suppressPackageStartupMessages(my_packages(
  c("devtools",
    "RCurl",
    "shiny",
    "DT",
    "shinydashboard",
    "shinyWidgets",
    "shinycssloaders",
    "tidyverse",
    "tidyselect",
    "fastcluster",
    "shinyjs",
    "colourpicker",
    "RColorBrewer",
    "patchwork",
    "zip",
    "ggpubr")
))

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 500MB. ----
options(shiny.maxRequestSize = 500 * 1024 ^ 2)

# LIST_DATA ----
LIST_DATA <<- list(
  table_file = NULL,
  # gene bin score set
  gene_file = NULL,
  # holds $Compleat genes from files and $gene file(s)
  gene_info = NULL,
  # for holding metaa data gene file(s) [c("gene_list", "set", "dot", "line", "color", plot?, nrom)]
  ttest = NULL,
  # t.test results $use is for numbers $gene_info for holding plotting options
  clust = NULL,
  # Cluster holder
  x_plot_range = c(0, 0),
  STATE = c(0, 0) # flow control,
  # [1] 1 = at least one file has been loaded and lets reactives fill in info
  #
  # [2] 0 = first time switching tab auto ploting
  #     1 = hiden plot button, reatives for plot enabled
  #     2 = on/off reactive picker changed, shows plot button, reatives for plot disabled
)

# types lines to be plotted ----
kLineOptions <-
  c(
    "solid line",
    "dashed line",
    "dotted line",
    "dot dash line",
    "long dash line",
    "two dash line",
    "No line"
  )

# types of dots to be used in plotting ----
kDotOptions <-
  c(
    "none",
    "square",
    "circle",
    "triangle point up",
    "diamond",
    "Biger circle",
    "smaller circle"
  )

# Brewer color sets to be avalible ----
kBrewerList <-
  c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral")

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

# color Brewer set that is active to use in plot ----
kListColorSet <- brewer.pal(11, kBrewerList[9])[-c(4:7)]

# math options avalible ----
kMathOptions <- c("mean", "sum", "median", "var")

# functions ----

# basic test for valid rgb color and type and
# switches hex <--> rgb, able to apply tint to rgb output from hex
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
    # tests if loading a file with a list of address to remote files, requirs .url.txt in file name
    for (i in seq_along(file_name)) {
      if (str_detect(file_name[i], ".url.txt")) {
        meta_data <- suppressMessages(read_delim(file_path[i],delim = " ", 
                                                 col_names = c("file","type","nick","color"),
                                                 col_types = cols(file=col_character(),
                                                                  type=col_character(),
                                                                  nick=col_character(),
                                                                  color=col_character())))
        my_remote_file <- c(my_remote_file, meta_data$file)
        my_color <- c(my_color, meta_data$color) 
        legend_nickname <- c(legend_nickname, str_replace(meta_data$nick,"\\.","_")) 
        my_landl <- last(meta_data$type)
        if(i > 1){
          file_path[i] <- NULL
        } else {
          file_path <- NULL
        }
      }else {
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
    # loop thourgh each item in file_path
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
      # checking/createing meta data
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
        my_color[x] <- sample(brewer.pal(11, sample(kBrewerList, size=1))[-c(4:7)],size = 1)
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
      # checks the number of bins and gene nameing skeem is consistant
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
        gene_names <-
          semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
        # test data is compatible with already loaded data
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
        mydot = kDotOptions[1],
        myline = kLineOptions[1],
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

# read in and match up names and change colors
LoadColorFile <- function(file_path, list_data, gene_list) {
  # check number oc columns
  num_bins <-
    count_fields(file_path,
                 n_max = 1,
                 tokenizer = tokenizer_delim(" "))
  # guessing list of colors
  if (num_bins == 1) {
    color_file <-
      suppressMessages(read_table(file_path, col_names = c("mycol", "set")))
  } else if(num_bins == 2) {
    color_file <-
      suppressMessages(read_table(file_path, col_names = c("set", "mycol")))
  } else {
    showModal(modalDialog(
      title = "Information message",
      paste(
        "I can't use this color file"
      ),
      size = "s",
      easyClose = TRUE
    ))
    return(list_data)
  }
  # check if color and convert to hex
  for (i in seq_along(color_file$mycol)) {
    if (is.na(color_file$mycol[i])) {
      color_file$mycol[i] <- sample(brewer.pal(11, sample(kBrewerList,size=1))[-c(4:7)],size = 1)
    } else {
      color_file$mycol[i] <- RgbToHex(color_file$mycol[i], convert = "hex")
    }
    if(is.na(color_file$set[i])){
      color_file$set[i] <- distinct(list_data$gene_info, set)$set[i] 
    }
  }
  color_file <- color_file %>% 
    dplyr::filter(!is.na(set))
  list_data <- full_join(list_data$gene_info,color_file,by="set") %>% 
    dplyr::filter(!is.na(gene_list)) %>%
    dplyr::mutate(mycol=if_else(!is.na(mycol.y),mycol.y,mycol.x)) %>% 
    dplyr::select(-mycol.y,-mycol.x)   
  return(list_data)
}

# records check box on/off
CheckBoxOnOff <- function(check_box, list_data) {
  if(!all(is.na(names(check_box)))){
    list_data <- full_join(list_data,check_box,by=c("set","gene_list")) %>%
      dplyr::filter(!is.na(set)) %>% 
      dplyr::mutate(onoff=if_else(is.na(onoff.y),"0",set)) %>% 
      dplyr::select(-onoff.y,-onoff.x) %>%
      distinct()
  }
  #### TODO
  # else {
  #   list_data <- mutate(list_data, onoff = "0")
  # }
  list_data
}

# make a new normalized file by deviding one file by the other
MakeNormFile <-
  function(list_data,
           nom,
           dnom,
           gbyg,
           divzerofix,
           addfiles,
           nickname) {
    # check 2 files have been selected
    if (nom == "" | dnom == "") {
      return(NULL)
    }
    # set up tool info and progress bar
    myname <- "bin_by_bin"
    # get data files
    nd <- list_data$table_file %>% dplyr::filter(set == nom | set== dnom)
    if(nchar(nickname)<1){
      nickname <- paste(nom, addfiles, dnom,sep = " ")
    }
    # applies custom norm factor(s)
    my_sel <- list_data$gene_info %>% 
      dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
      dplyr::filter(set == nom | set== dnom)
    new_gene_list <- full_join(dplyr::filter(nd,set == nom), 
                               dplyr::filter(nd,set == dnom), by = c("gene", "bin")) %>% 
      replace_na(., list(score = 0))
    legend_nickname <- nickname
    if (addfiles == "+") {
      new_gene_list <- transmute(
        new_gene_list,
        gene = gene,
        bin = bin,
        set = legend_nickname,
        score = score.x + score.y
      )
    } else {
      # if min/2 find Na's and 0's, and replace
      if (divzerofix) {
        new_gene_list <- new_gene_list %>% 
          dplyr::mutate(score=if_else(set == dnom & score == 0,NA,score))
        myname <- paste0(myname, "_0->min/2")
        new_min_for_dom <-
          min(new_gene_list$score, na.rm = TRUE) / 2
        new_gene_list <-
          replace_na(new_gene_list, list(score = new_min_for_dom))
      }
      # files numbers are replaced with mean of bins if applied
      if (gbyg != "bin by bin") {
        myname <- "mean_of_bins"
        if (divzerofix) {
          myname <- paste0(myname, "_0->min/2")
        }
        new_gene_list <- new_gene_list %>% 
          group_by(bin, set) %>%
          dplyr::mutate(score = mean(score, na.rm = TRUE)) %>% ungroup()
      }
      # make gene list and do math
      legend_nickname <- paste0(nickname, ": ", myname)
      new_gene_list <- transmute(
        new_gene_list,
        gene = gene,
        bin = bin,
        set = legend_nickname,
        score = score.x / score.y
      ) %>% na_if(Inf)
    }
    # output test
    if (n_distinct(new_gene_list$gene) < 1) {
      showModal(
        modalDialog(
          title = "Information message",
          " No genes left, try replacing Inf and/or bin by bin",
          size = "s",
          easyClose = TRUE
        )
      )
      return(NULL)
    }
    # adds meta data 
    list_data$table_file <- bind_rows(list_data$table_file, new_gene_list)
    list_data$gene_info <- bind_rows(list_data$gene_info,
                                     tibble(
                                       gene_list = names(list_data$gene_file)[1],
                                       set = legend_nickname,
                                       mydot = kDotOptions[1],
                                       myline = kLineOptions[1],
                                       mycol = sample(brewer.pal(11, sample(kBrewerList,size=1))[-c(4:7)],size = 1),
                                       onoff = "0",
                                       sub = " ",
                                       plot_set = " "
                                     ))
    return(list_data)
  }

# removes gene list
RemoveGeneList <-
  function(list_data, list_name) {
    list_data$gene_file[[list_name]] <- NULL
    list_data$gene_info <- list_data$gene_info %>% 
      dplyr::filter(gene_list != list_name)
    if (list_data$STATE[2] != 0) {
      list_data$STATE[2] <- 2
    }
    list_data
  }

# removes data file
RemoveFile <- function(list_data, file_name, remove_all) {
  if (n_distinct(list_data$table_file$set) > 1 & !remove_all) {
    list_data$gene_info <- list_data$gene_info %>% 
      dplyr::filter(set != file_name)
    list_data$table_file <- list_data$table_file %>% 
      dplyr::filter(set != file_name)
    # builds new main in Compleat gene lists
    list_data$gene_file[[1]]$use <-
      distinct(list_data$table_file, gene)
    my_name <-
      paste("Compleat\nn =", n_distinct(list_data$gene_file[[1]]$use$gene))
    names(list_data$gene_file)[1] <- my_name
    if (list_data$STATE[2] != 0) {
      list_data$STATE[2] <- 2
    }
  } else {
    # builds empty list of lists
    list_data <- list(
      table_file = NULL,
      # gene bin score set
      gene_file = NULL,
      # holds $Compleat genes from files and $gene file(s)
      gene_info = tibble(
        gene_list = "none",
        set = "none",
        mydot = kDotOptions[1],
        myline = kLineOptions[1],
        mycol = sample(brewer.pal(11, sample(kBrewerList,size=1))[-c(4:7)],size = 1),
        onoff = "0",
        sub = " ",
        plot_set = " "
      ),
      # for holding metaa data gene file(s) [c("gene_list", "set", "dot", "line", "color", plot?, nrom)]
      ttest = NULL,
      # t.test results $use is for numbers $gene_info for holding plotting options
      clust = NULL,
      # Cluster holder
      x_plot_range = c(0, 0),
      STATE = c(0, 0) # flow control,
      # [1] 1 = at least one file has been loaded and lets reactives fill in info
      #
      # [2] 0 = first time switching tab auto ploting
      #     1 = hiden plot button, reatives for plot enabled
      #     2 = on/off reactive picker changed, shows plot button, reatives for plot disabled
    )
  }
  list_data
}

# Total, exclusive and intersected gene lists
IntersectGeneLists <-
  function(list_data, list_name) {
    if (is.null(list_name)) {
      return(NULL)
    }
    setProgress(1, detail = paste("building list"))
    outlist <- NULL
    # grab selected gene list(s)
    lapply(list_name, function(j) {
      outlist[[j]] <<- list_data$gene_file[[j]]$use %>% dplyr::select(gene)
    })
    # calapes into one list
    outlist <- bind_rows(outlist)
    # remove any pre used data
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Gene_List_")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Gene_List_"))
    
    # recored for info
    setProgress(2, detail = paste("building Total list"))
    if (n_distinct(outlist$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_Total\nn =", n_distinct(outlist$gene))
      # recored for info
      list_data$gene_file[[nick_name1]]$full <- distinct(outlist)
      list_data$gene_file[[nick_name1]]$use <- distinct(dplyr::select(outlist, gene))
      list_data$gene_file[[nick_name1]]$info <-
        paste("Gene_List_Total",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date())
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_Total"), 
                                           onoff = "0",
                                           plot_set = " ")))
    }
    setProgress(3, detail = paste("building intersect list"))
    intersected <- dplyr::filter(outlist, duplicated(gene))
    if (n_distinct(intersected$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_intersect\nn =", n_distinct(intersected$gene))
      # recored for info
      list_data$gene_file[[nick_name1]]$full <- intersected
      list_data$gene_file[[nick_name1]]$use <- dplyr::select(intersected, gene)
      list_data$gene_file[[nick_name1]]$info <-
        paste("Gene_List_intersect",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date())
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_intersect"), 
                                           onoff = "0",
                                           plot_set = " ")))
      setProgress(4, detail = paste("building exclusive list"))
      exclusive <- anti_join(outlist, intersected, by = "gene")
      if (n_distinct(exclusive$gene) == 0) {
        exclusive <- distinct(outlist)
      }
      nick_name1 <-
        paste("Gene_List_exclusive\nn =", n_distinct(exclusive$gene))
      # recored for info
      list_data$gene_file[[nick_name1]]$full <- exclusive
      list_data$gene_file[[nick_name1]]$use <- dplyr::select(exclusive, gene)
      list_data$gene_file[[nick_name1]]$info <-
        paste("Gene_List_exclusive",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date())
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste("Gene_List_exclusive"), 
                                           onoff = "0",
                                           plot_set = " ")))
      
    }
    list_data
  }

# sorts active gene list contain top % signal based on selected bins and file
FilterTop <-
  function(list_data,
           list_name,
           file_names,
           start_bin,
           end_bin,
           num,
           topbottom,
           my_filter_all) {
    if (is.null(file_names)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return(NULL)
    }
    lc <- 0
    outlist <- NULL
    lapply(file_names, function(j) {
      setProgress(lc + 1, detail = paste("sorting", j))
      apply_bins <-
        semi_join(dplyr::filter(list_data$table_file, set == j), 
                  list_data$gene_file[[list_name]]$use, by = 'gene')  
      apply_bins <- group_by(apply_bins, gene) %>%
        dplyr::filter(bin %in% start_bin:end_bin) %>%
        summarise(mysums = sum(score, na.rm = TRUE),.groups="drop") %>%
        mutate(myper = as.numeric(strtrim(cume_dist(mysums), 5))) %>%
        arrange(desc(mysums))
      gene_count <- nrow(apply_bins)
      if (topbottom == "Top%") {
        num2 <- c(1, ceiling(gene_count * (num / 100)))
      } else if (topbottom == "Middle%") {
        if (num == 100) {
          num2 <- c(1, gene_count)
        } else {
          med <- median(apply_bins$myper)
          num2 <-
            c(count(apply_bins, myper >= max(med, num / 100))[[2]][2],
              count(apply_bins, myper <= min(med, (100 - num) / 100))[[2]][1])
        }
      } else {
        num2 <-
          c(ceiling((gene_count + 1) - (gene_count * (num / 100))), gene_count)
      }
      if (any(is.na(num2))) {
        num2 <-
          c(ceiling((gene_count) - (gene_count * max(.5, num / 100))), ceiling(gene_count * max(.5, num / 100)))
      }
      outlist2 <- dplyr::mutate(apply_bins,!!j := myper) %>%
        dplyr::select(gene,!!j) %>%
        slice(num2[1]:num2[2])
      if (lc > 0) {
        if(my_filter_all){
          outlist <<- inner_join(outlist, outlist2, by = 'gene')
        } else{
          outlist <<- full_join(outlist, outlist2, by = 'gene')
        }
      } else {
        outlist <<- outlist2
      }
      lc <<- lc + 1
    })
    if (length(outlist$gene) == 0) {
      return(list_data)
    }
    old_name <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_name) > 3) {
      # remove old sort gene list keepint 4
      list_data$gene_file[[first(old_name)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_name))
    }
    setProgress(lc + 2, detail = "building list")
    topbottom2 <- paste(str_remove(topbottom,"%"), paste0(num, "%"))
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter ",topbottom2, "\nn = ", n_distinct(outlist$gene))), 33)
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- dplyr::select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <-
      paste(
        "Filter",
        topbottom2,
        "bins",
        start_bin,
        "to",
        end_bin,
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date(),
        list_data$gene_file[[list_name]]$info
      )
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter",
                                                      topbottom2,
                                                      "bins",
                                                      start_bin,
                                                      "to",
                                                      end_bin), 
                                         onoff = "0",
                                         plot_set = " ")))
    list_data
  }

# sort my percent
FilterPer <-
  function(list_data,
           list_name,
           file_names,
           start_end_bin,
           my_filter_all,
           my_per,
           my_type,
           anyall) {
    if (is.null(file_names)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return(NULL)
    }
    p_funs <- map(my_per/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
      set_names(paste0("my_p_",seq_along(my_per)))
    gene_list <- list_data$gene_file[[list_name]]$use
    out_list <- list_data$table_file %>% 
      dplyr::filter(set == file_names) %>% 
      semi_join(.,gene_list,by="gene") %>% 
      dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2]) 
    
    out_per <- out_list %>%
      group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() 
    
    if(my_type == "min%"){
      if(anyall){
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score >= my_p_1)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene,set) %>% dplyr::filter(all(score >= my_p_1)) %>% 
          ungroup() %>% distinct(gene)
      }
      topbottom2 <- paste(str_remove(my_type,"%"), paste0(my_per[1], "%"))
    } else if(my_type == "max%"){
      if(anyall){
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene,set) %>% dplyr::filter(all(score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      }
      topbottom2 <- paste(str_remove(my_type,"%"), paste0(my_per[2], "%"))
    } else {
      if(anyall){
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene) %>% dplyr::filter(all(score >= my_p_1 & score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      } else {
        out_list <- full_join(out_list,out_per,by=c("bin","set")) %>%
          group_by(gene,set) %>% dplyr::filter(all(score >= my_p_1 & score <= my_p_2)) %>% 
          ungroup() %>% distinct(gene)
      }
      topbottom2 <- paste(paste(str_remove(my_type,"%"), paste0(my_per[1], "%")),paste0(my_per[2], "%"),collapse = " and ")
    }
    if (length(out_list$gene) == 0) {
      return(list_data)
    }
    old_names <- grep("^Filter", names(list_data$gene_file), value = T)
    if (length(old_names) > 3) {
      # remove old sort gene list keeping 4
      list_data$gene_file[[first(old_names)]] <- NULL
      list_data$gene_info <- dplyr::filter(list_data$gene_info,
                                           gene_list != first(old_names))
    }
    nick_name <-
      strtrim(gsub("(.{30})",
                   "\\1... ",
                   paste0("Filter Prob ",topbottom2, "\nn = ", n_distinct(out_list$gene))), 33)
    if(length(file_names) == 1){
      my_per <- seq(my_per[1],my_per[2]/2.1,length.out = 5)
      p_funs <- map(my_per/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per)
      my_list <- list_data$table_file %>% 
        dplyr::filter(set == file_names) %>% 
        semi_join(.,gene_list,by="gene") %>% 
        dplyr::filter(bin %in% start_end_bin[1]:start_end_bin[2]) 
      out_per1 <- my_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        gather(.,key = set,value = "my_p_1",-bin,-set)
      
      my_per <- seq(my_per[2]/2,my_per[2],length.out = 5)
      p_funs <- map(my_per/100, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
        set_names(my_per)
      
      out_per2 <- my_list %>%
        group_by(set,bin) %>% summarize_at(vars(score), p_funs) %>% ungroup() %>% 
        gather(.,key = set,value = "my_p_2",-bin,-set)
      out_list1 <- full_join(out_per1,out_per2,by=c("bin","set")) %>% 
        dplyr::mutate(set=paste0(round(as.numeric(set),2),"%"))
    } else {
      out_list1 <- list_data$table_file %>% 
        dplyr::filter(set == file_names) %>% 
        semi_join(.,gene_list,by="gene") %>% 
        full_join(.,out_per,by=c("bin","set")) %>% 
        replace_na(list(my_p_1 = 0, my_p_2 = 0)) %>% 
        dplyr::select(-gene,-score) 
    }
    
    list_data$gene_file[[nick_name]]$full <- out_list %>% dplyr::mutate(min=my_per[1],max=my_per[2])
    list_data$gene_file[[nick_name]]$use <- out_list
    list_data$sortplot <- out_list1
    list_data$gene_file[[nick_name]]$info <-
      paste(
        "Filter Prob:",
        topbottom2,
        "bins",
        start_end_bin[1],
        "to",
        start_end_bin[2],
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date(),
        list_data$gene_file[[list_name]]$info
      )
    list_data$gene_info <-
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>%
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>%
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste("Filter Prob:",
                                                      topbottom2,
                                                      "bins",
                                                      start_end_bin[1],
                                                      "to",
                                                      start_end_bin[2]),
                                         onoff = "0",
                                         plot_set = " ")))
    
    list_data
  }

# a[1]/b[2] or (a[1]/a[2])/(b[1]/b[2]) make gene list
CompareRatios <-
  function(list_data,
           list_name,
           ratio1file,
           ratio2file,
           start1_bin,
           end1_bin,
           start2_bin,
           end2_bin,
           num,
           divzerofix,
           normbin = 0) {
    if (ratio1file == "") {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    outlist <- NULL
    if (ratio2file == "None" | ratio2file == "") {
      if(sum(start2_bin,end2_bin) == 0){
        showModal(modalDialog(
          title = "Information message",
          paste("no file or bins to compare to"),
          size = "s",
          easyClose = TRUE
        ))
        return()
      }
      ratiofile <- ratio1file
    } else {
      ratiofile <- c(ratio1file, ratio2file)
    }
    lc <- 0
    lapply(ratiofile, function(j) {
      df <-
        semi_join(dplyr::filter(list_data$table_file, set == j), 
                  list_data$gene_file[[list_name]]$use, by = 'gene') 
      if (normbin > 0) {
        df <- group_by(df, gene) %>%
          arrange(bin) %>% 
          dplyr::mutate(score = score / nth(score, normbin))
      }
      df <- group_by(df, gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T),.groups="drop") %>%
        ungroup()
      # if min/2 find Na's and 0's, and replace
      if(sum(start2_bin,end2_bin) == 0){
        df$sum2 <- 1
      }
      if (divzerofix) {
        df$sum2 <- na_if(df$sum2, 0)
        new_min <-
          min(df$sum2, na.rm = TRUE) / 2
        df <-
          replace_na(df, list(sum2 = new_min))
      }
      lc <<- lc + 1
      outlist[[lc]] <<-
        transmute(df, gene = gene, Ratio = sum1 / sum2) %>%
        na_if(Inf) %>% dplyr::select(gene, Ratio)
      
      if (lc > 1) {
        if (divzerofix) {
          outlist[[2]]$Ratio <- na_if(outlist[[2]]$Ratio, 0)
          outlist[[2]] <-
            replace_na(outlist[[2]], list(Ratio = new_min))
          outlist[[1]] <<-
            inner_join(outlist[[1]], outlist[[2]], by = 'gene') %>%
            transmute(gene = gene, Ratio = Ratio.x / Ratio.y) %>%
            na_if(Inf)  %>% dplyr::select(gene, Ratio)
        } else {
          outlist[[1]] <<-
            inner_join(outlist[[1]], outlist[[2]], by = 'gene') %>%
            transmute(gene = gene, Ratio = Ratio.x / Ratio.y) %>%
            na_if(Inf)  %>% dplyr::select(gene, Ratio)
        }
      }
    })
    #remove old info
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Ratio_")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Ratio_"))
    if(num < 0){
      num <- 1/num
    }
    nick_name <- NULL
    if(num != 0){
      upratio <- dplyr::filter(outlist[[1]], Ratio > num)
    } else {
      upratio <- NULL
    }
    
    if (n_distinct(upratio$gene) > 0) {
      nick_name1 <-
        paste("Ratio_Up_file1\nn =", n_distinct(upratio$gene))
      nick_name <- c(nick_name, nick_name1)
      list_data$gene_file[[nick_name1]]$full <- upratio %>% dplyr::mutate(set=nick_name1)
      list_data$gene_file[[nick_name1]]$use <- dplyr::select(upratio, gene)
      list_data$gene_file[[nick_name1]]$info <-
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "]/",
          ratio1file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "] ",
          "fold change cut off",
          num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name1,
                                           sub =  paste(
                                             "\nRatio_Up_file1",
                                             "fold change cut off",
                                             num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           plot_set = " ")))
    }
    if(num != 0){
      upratio <- dplyr::filter(outlist[[1]], Ratio < 1 / num & Ratio != 0)
    } else {
      upratio <- NULL
      } 
    if (n_distinct(upratio$gene) > 0) {
      nick_name2 <-
        paste("Ratio_Down_file1\nn =", n_distinct(upratio$gene))
      nick_name <- c(nick_name, nick_name2)
      list_data$gene_file[[nick_name2]]$full <- upratio %>% dplyr::mutate(set=nick_name2)
      list_data$gene_file[[nick_name2]]$use <- dplyr::select(upratio, gene)
      list_data$gene_file[[nick_name2]]$info <-
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "]/",
          ratio1file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "] ",
          "fold change cut off",
          num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name2,
                                           sub =  paste(
                                             "Ratio_Down_file1",
                                             "fold change cut off",
                                             num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           plot_set = " ")))
    }
    if(num != 0){
      upratio <-
        dplyr::filter(outlist[[1]], Ratio <= num &
                        Ratio >= 1 / num | Ratio == 0)
    } else {
      upratio <- outlist[[1]]
    }
    
    if (n_distinct(upratio$gene) > 0) {
      nick_name3 <-
        paste("Ratio_No_Diff\nn =", n_distinct(upratio$gene))
      nick_name <- c(nick_name, nick_name3)
      list_data$gene_file[[nick_name3]]$full <- upratio %>% dplyr::mutate(set=nick_name3)
      list_data$gene_file[[nick_name3]]$use <- dplyr::select(upratio, gene)
      list_data$gene_file[[nick_name3]]$info <-
        paste(
          "Ratio_Up_file1",
          ratio2file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "]/",
          ratio1file,
          "[",
          start1_bin,
          "to",
          end1_bin,
          "]/[",
          start2_bin,
          "to",
          end2_bin,
          "] ",
          "fold change cut off",
          num,
          divzerofix,
          "from",
          list_name,
          "gene list",
          Sys.Date()
        )
      list_data$gene_info <- 
        distinct(bind_rows(list_data$gene_info,
                           list_data$gene_info %>% 
                             dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                             dplyr::mutate(gene_list = nick_name3,
                                           sub =  paste(
                                             "Ratio_No_Diff",
                                             "fold change cut off",
                                             num,
                                             divzerofix,
                                             "from",
                                             list_name
                                           ), 
                                           onoff = "0",
                                           plot_set = " ")))
    }
    list_data$boxRatio <- NULL
    for (nn in nick_name) {
      list_data$boxRatio <- bind_rows(list_data$boxRatio,list_data$gene_file[[nn]]$full)
    }
    list_data
  }

# Change the number of clusters
ClusterNumList <- function(list_data,
                           list_name,
                           clusterfile,
                           start_bin,
                           end_bin,
                           num,
                           myname) {
  if (is_empty(list_data$clust) | clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  setProgress(3, detail = "spliting into clusters")
  list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^Group_|^Cluster_")]
  list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^Group_|^Cluster_"))
  if (myname == "Cluster_") {
    gene_list <-
      dplyr::mutate(list_data$clust$use, cm = cutree(list_data$clust$cm, num))
  } else {
    gene_list <-
      dplyr::mutate(list_data$clust$full, cm = ntile(cm, as.numeric(num)))
  }
  if (is_empty(grep("Cluster_", list_name)) |
      is_empty(grep("Group_", list_name))) {
    list_name <- 1
  }
  for (nn in 1:num) {
    outlist <- dplyr::filter(gene_list, cm == nn)
    nick_name <-
      paste(paste0(myname, nn, "\nn ="), n_distinct(outlist$gene))
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- dplyr::select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <-
      paste(
        nick_name,
        "bins",
        start_bin,
        "to",
        end_bin,
        "from",
        list_name,
        clusterfile,
        num,
        myname,
        "total",
        Sys.Date()
      )
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(list_data$gene_file)[1]) %>% 
                           dplyr::mutate(gene_list = nick_name,
                                         sub =  paste(
                                           myname,
                                           "by",
                                           num,
                                           "from",
                                           list_name
                                         ), 
                                         onoff = "0",
                                         plot_set = " ")))
    setProgress(4, detail = paste("finishing cluster", nn))
  }
  list_data
}

# finds 2 - 4 clusers from the one active file, plotting the patterns and displaying the gene lists
FindClusters <- function(list_data,
                         list_name,
                         clusterfile,
                         start_bin,
                         end_bin) {
  if (clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  setProgress(1, detail = paste("gathering data"))
  df <-
    semi_join(dplyr::filter(list_data$table_file, set == clusterfile), 
              list_data$gene_file[[list_name]]$use, by = 'gene') 
  setProgress(2, detail = "hierarchical clustering using ward method")
  list_data$clust <- list()
  list_data$clust$cm <-
    hclust.vector(as.data.frame(spread(df, bin, score))[, c((start_bin:end_bin) + 2)], method = "ward")
  list_data$clust$use <- distinct(df, gene)
  list_data
}

# finds 2 - 4 groups from the one active file, plotting the patterns and displaying the gene lists
FindGroups <- function(list_data,
                       list_name,
                       clusterfile,
                       start_bin,
                       end_bin) {
  if (clusterfile == "") {
    showModal(modalDialog(
      title = "Information message",
      paste("No file selected to work on"),
      size = "s",
      easyClose = TRUE
    ))
    return(NULL)
  }
  setProgress(1, detail = paste("gathering data"))
  df <-
    semi_join(dplyr::filter(list_data$table_file, set == clusterfile), 
              list_data$gene_file[[list_name]]$use, by = 'gene') 
  setProgress(2, detail = "finding groups")
  list_data$clust <- list()
  list_data$clust$full <- group_by(df, gene) %>%
    dplyr::filter(bin %in% start_bin:end_bin) %>%
    summarise(cm = sum(score, na.rm = TRUE),.groups="drop")
  list_data$clust$use <- distinct(df, gene)
  list_data
}

# Cumulative Distribution data prep "PI / EI"
CumulativeDistribution <-
  function(list_data,
           onoffs,
           start1_bin,
           end1_bin,
           start2_bin,
           end2_bin) {
    print("cdf function")
    if (is.null(onoffs)) {
      showModal(modalDialog(
        title = "Information message",
        paste("No file selected to work on"),
        size = "s",
        easyClose = TRUE
      ))
      return()
    }
    # remove old data sets
    list_data$gene_file <- list_data$gene_file[!str_detect(names(list_data$gene_file),"^CDF")]
    list_data$gene_info <- list_data$gene_info %>% dplyr::filter(!str_detect(gene_list,"^CDF"))
    outlist <- NULL
    for (list_name in names(onoffs)) {
      setProgress(1, detail = paste("dividing one by the other"))
        # Complete within gene list and sum regions
      outlist[[list_name]] <-
          semi_join(dplyr::filter(list_data$table_file, set == onoffs[[list_name]]), 
                    list_data$gene_file[[list_name]]$use, by = 'gene') %>% 
          group_by(gene,set) %>%
          summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                    sum2 = sum(score[start2_bin:end2_bin],	na.rm = T),.groups="drop") %>% 
          ungroup() %>% 
          dplyr::mutate(., value = sum1 / sum2) %>%
          na_if(Inf) %>% na_if(0) %>% 
          group_by(., gene,set) %>%
          dplyr::mutate(test = sum(value)) %>%
          ungroup() %>% 
          dplyr::filter(!is.na(test)) %>%
          group_by(., set) %>%
          arrange(value) %>%
          dplyr::transmute(
            gene=gene,
            bin = row_number(),
            set = set,
            plot_set = paste(list_name, "-", set),
            value = value
          ) %>%
          ungroup()
    }
    
    # unlist and binds all together
    outlist <- bind_rows(outlist) %>% distinct()
    
    
    setProgress(2, detail = paste("building list"))
    # removes top and bottom %
    if (sum(start1_bin, end1_bin) > sum(start2_bin, end2_bin)) {
      use_header <- "Log2 EI Cumulative plot"
    } else {
      use_header <- "Log2 PI Cumulative plot"
    }
    if (n_distinct(outlist$gene) > 0) {
      nick_name1 <- paste0("CDF ", use_header)
      list_data$gene_file[[nick_name1]]$full <- outlist
      list_data$gene_file[[nick_name1]]$use <- outlist %>% select(gene)
      list_data$gene_file[[nick_name1]]$info <-
        paste(
          use_header,
          "CDF",
          "bins",
          start1_bin,
          "to",
          end1_bin,
          "/",
          start2_bin,
          "to",
          end2_bin,
          "from",
          names(onoffs),
          "gene list(s)",
          paste(distinct(outlist, plot_set), collapse = " "),
          Sys.Date()
        )
    } else {
      nick_name1 <- paste("CDF n = 0")
    }
    for (list_name in names(onoffs)) {
    list_data$gene_info <- 
      distinct(bind_rows(list_data$gene_info,
                         list_data$gene_info %>% 
                           dplyr::filter(gene_list == names(LIST_DATA$gene_file)[1] &
                                           set %in% onoffs[[list_name]]) %>% 
                           dplyr::mutate(gene_list = nick_name1,
                                         sub =  paste(
                                           "CDF from",
                                           list_name
                                         ), 
                                         onoff = "0",
                                         plot_set = paste(list_name, "-", set),
                                         myheader = use_header)))
    setProgress(5, detail = "finishing up")
    }
    list_data
  }

# makes sure t test wont crash on an error
try_t_test <- function(db,my_set,my_math ="none",my_test="t.test",padjust="fdr",
                       alternative="two.sided", exact=FALSE, paired=FALSE){
  print("try t.test")
  exact <- if_else(exact=="TRUE",TRUE,FALSE)
  paired <- if_else(paired=="TRUE",TRUE,FALSE)
  combn(unique(db$set),2) -> my_comparisons
  my_comparisons2 <- list()
  db_out <- list()
  for(cc in 1:ncol(my_comparisons)){
    my_comparisons2[[cc]] <- (c(my_comparisons[1,cc],my_comparisons[2,cc]))
  }
  db <- spread(db,set,score) 
  for(i in my_comparisons2){
    db2 <- dplyr::select(db, gene,bin,all_of(i)) %>% 
      rename(score.x=all_of(names(.)[3]), score.y=all_of(names(.)[4])) 
    
    myTtest <- tibble(bin=NA,p.value=NA)
    
    for(t in unique(db2$bin)){
      x.score <- dplyr::filter(db2, bin ==t)
      y.score <- dplyr::filter(db2, bin ==t)
      kk <- try(get(my_test)(x.score$score.x,y.score$score.y,
                             alternative = alternative, 
                             exact=exact, 
                             paired=paired)$p.value)
      if("try-error" %in% class(kk) | !is.numeric(kk)){
        kk <-1
      }
      myTtest <- myTtest %>% add_row(bin = t, p.value=kk)
    }
    myTtest <- myTtest %>% dplyr::filter(!is.na(bin))
    if(padjust != "NO"){
      myTtest <- myTtest %>% dplyr::mutate(p.value=p.adjust(p.value,method = padjust))
    }
    if(my_math =="-log"){
      myTtest <- myTtest %>% dplyr::mutate(p.value=if_else(p.value==0,2.2e-16,p.value)) %>% dplyr::mutate(p.value=-log(p.value))
    } else if(my_math =="-log10"){
      myTtest <- myTtest %>% dplyr::mutate(p.value=if_else(p.value==0,2.2e-16,p.value)) %>% dplyr::mutate(p.value=-log10(p.value))
    }
    db_out[[str_c(i,collapse = "-")]] <- myTtest%>% 
      dplyr::mutate(., set = paste(
        paste0(gsub("(.{20})", "\\1\n", 
                    str_split_fixed(str_c(i,collapse = "-"), "\n",n=2)[,1]),
               gsub("(.{20})", "\\1\n", 
                    str_split_fixed(str_c(i,collapse = "-"), "\n",n=2)[,2])),
        gsub("(.{20})", "\\1\n", 
             str_split_fixed(my_set, "\n", n=2)[,1]),
        str_split_fixed(my_set, "\n", n=2)[,2],
        sep = '\n'
      ))
    
  }
  
  db_out
}

# makes tibble with all active data
Active_list_data <-
  function(list_data) {
    table_file <- list_data$table_file
    gene_file <- list_data$gene_file
    gene_info <- list_data$gene_info
    list_data_out <- NULL
    print("active data function")
    for ( i in names(gene_file) ){
      # checks to see if at least one file in list is acitve
      if (gene_info %>% dplyr::filter(gene_list == i & onoff != 0) %>% nrow() == 0) {
        next()
      } else {
        my_sel <- gene_info %>% dplyr::filter(gene_list == i & onoff != 0)
        list_data_out[[i]] <-
          table_file %>% 
          dplyr::filter(set %in% my_sel$onoff) %>%
          semi_join(., gene_file[[i]]$use, by = "gene") %>% 
          dplyr::mutate(., gene_list = i)
        my_sel2 <- my_sel %>% dplyr::mutate(.,plot_set = paste(
          gsub("(.{20})", "\\1\n", 
               str_split_fixed(i, "\nn = ", n=2)[,1]),
          paste0("n = ", n_distinct(list_data_out[[i]]$gene)),
          gsub("(.{20})", "\\1\n", set),
          sep = '\n'
        )) %>% select(set,plot_set)
        list_data_out[[i]] <- list_data_out[[i]] %>% inner_join(.,my_sel2,by="set")
      }
    }
    return(bind_rows(list_data_out))
  }

# Applys math to active data in each gene list
ApplyMath <-
  function(list_data,
           use_math,
           relative_frequency,
           normbin) {
    print("apply math fun")
    setProgress(1, detail = paste("Gathering info"))
    # applys math to data file
    if (relative_frequency == "rel gene frequency") {
      list_data <- list_data %>% group_by(plot_set, gene) %>%
        dplyr::mutate(score = score / sum(score, na.rm = TRUE)) %>%
        ungroup()
    }
    list_data <- list_data %>% group_by(plot_set, bin) %>%
      summarise(value = get(use_math)(score, na.rm = T), .groups="drop")
    
    if (normbin > 0) {
      list_data <- list_data %>% 
        group_by(plot_set) %>%
        arrange(bin) %>%
        dplyr::mutate(value = value / nth(value, normbin)) %>%
        ungroup()
    } else if (relative_frequency == "relative frequency") {
      list_data <- list_data %>%
        group_by(plot_set) %>%
        dplyr::mutate(value = value / sum(value)) %>%
        ungroup()
    }
    return(list_data %>% dplyr::mutate(set=plot_set))
  }

# applys t.test to active data
ApplyTtest <-
  function(list_data,
           switchttest,
           use_tmath,
           switchttesttype,
           padjust,my_alt, my_exact, my_paired) {
    print("apply t test")
    # t.test comparing files in same gene list
    ttest <- NULL
    if(switchttest == "by lists" & n_distinct(list_data$gene_list) > 1){
      list_data <- list_data %>% 
        rename(set=gene_list,gene_list=set)
    }
    n_test <- list_data %>% group_by(gene_list) %>% 
      summarise(n=n_distinct(set), .groups= "drop")
    if(switchttest != "none"){
      for(i in distinct(list_data, gene_list)$gene_list){
        if(n_test %>% dplyr::filter(gene_list == i) %>% dplyr::select(n) > 1){
          ttest[[i]] <- bind_rows(try_t_test(list_data %>% dplyr::filter(gene_list == i),i,
                                             use_tmath,switchttesttype,padjust,
                                             my_alt, noquote(my_exact), noquote(my_paired))) %>% 
            dplyr::mutate(mydot = kDotOptions[1],
                          myline = kLineOptions[1],
                          mycol = "#000000" )
        } 
      }
    } 
    bind_rows(ttest)
  }

# gather relavent plot option data
MakePlotOptionFrame <- function(list_data) {
  print("plot options fun")
  gene_info <- list_data$gene_info
  list_data_frame <- NULL
  options_main <- NULL
  # checks to see if at least one file in list is acitve
  if (gene_info %>% dplyr::filter(onoff != 0) %>% nrow() == 0) {
    return(NULL)
  } else {
    gene_info <- gene_info %>% dplyr::filter(onoff != 0)
    my_lines <-
      match(gene_info$myline, kLineOptions)
    my_dots <-
      match(gene_info$mydot, kDotOptions)
    gene_info <- gene_info %>%
      dplyr::mutate(
        myline = if_else(my_lines > 6, 0, as.double(my_lines)),
        mydot = if_else(my_dots == 1, 0, my_dots + 13),
        mysizedot = if_else(my_dots == 1, 0.01, 4.5),
        set = plot_set
      )
  }
  # tint if same color is used more then once
  ldf <- duplicated(gene_info$mycol)
  for (i in seq_along(gene_info$mycol)) {
    if (ldf[i]) {
      # mycol = RgbToHex(color_file$color[i], convert = "hex", tint = mytint)
      gene_info$mycol[i] <- RgbToHex(gene_info$mycol[i], convert = "hex", tint = log(i,10))
    }
  }
  return(gene_info)
}

# gather relavent plot option data
MakePlotOptionttest <- function(list_data, Y_Axis_TT,my_ttest_log,hlineTT,pajust,ttype) {
  if(is_empty(list_data)){
    return(NULL)
  }
  print("plot options ttest fun")
  out_options <- group_by(list_data,set) %>% distinct(.,myline,mydot,mycol) %>% 
    dplyr::mutate(myline=match(myline,kLineOptions)) %>% 
    dplyr::mutate(mydot=match(mydot, kDotOptions)) %>% 
    dplyr::mutate(if_else(mydot == 1, 0.01, 4.5)) %>% 
    dplyr::mutate(myline=if_else(myline > 6, 0, as.double(myline))) %>% 
    dplyr::mutate(mydot=if_else(mydot == 1, 0, mydot+13))
  list_data_frame <- NULL
  
  ldf <- duplicated(out_options$mycol)
  for (i in seq_along(out_options$mycol)) {
    if (ldf[i]) {
      out_options$mycol[i] <-
        RgbToHex(out_options$mycol[i], convert = "hex", tint = log(i,10))
    }
  }
  
  list_data_frame$options_main_tt <- out_options
  list_data_frame$ylimTT <- Y_Axis_TT
  if(pajust != "none"){
    pp <- paste0("p.adjust:", pajust)
  } else{
    pp <- "p.value"
  }
  if(my_ttest_log == "-log"){
    list_data_frame$hlineTT <- -log(hlineTT)
    list_data_frame$ylabTT <- paste0(my_ttest_log,"(",pp,")"," ",ttype)
  } else if(my_ttest_log == "-log10"){
    list_data_frame$hlineTT <- -log10(hlineTT)
    list_data_frame$ylabTT <- paste0(my_ttest_log,"(",pp,")"," ",ttype)
  } else{
    list_data_frame$hlineTT <- hlineTT
    list_data_frame$ylabTT <- paste0(pp," ",ttype)
  }
  
  return(list_data_frame)
}

# Sets y lable fix
YAxisLable <-
  function(use_math = "mean",
           relative_frequency = "none",
           norm_bin = 0,
           smoothed = F,
           log_2 = F) {
    use_y_label <- paste(use_math, "of bin counts")
    if (relative_frequency == "rel gene frequency") {
      use_y_label <- paste("RF per gene :", use_y_label)
    } else if (relative_frequency == "relative frequency") {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : RF")
    }
    if (norm_bin > 0) {
      if (relative_frequency == "rel gene frequency") {
        use_y_label <- paste(use_y_label, " : Norm bin ", norm_bin)
      } else {
        use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                             "bins : Normalize to bin ",
                             norm_bin)
      }
    }
    if (log_2) {
      use_y_label <- paste0("log2(", use_y_label, ")")
    }
    if (smoothed) {
      use_y_label <- paste0("smoothed(", use_y_label, ")")
    }
    use_y_label
  }

# Sets lines and labels
LinesLabelsListset <- function(body1bin = 20,
                               body2bin = 40,
                               tssbin = 15,
                               tesbin = 45,
                               binbp = 100,
                               totbins = 80,
                               everybin = 5,
                               tssname = "TSS",
                               tesname = "pA") {
  print("lines and labels fun")
  # I creat this in steps bin 1 to next land mark (TSS TES) then go from there to next land mark until end of bins
  everybp <- everybin * binbp
  if (everybp > 0) {
    my_5prim <- NULL
    my_3prim <- NULL
    if (tssbin > 0) {
      # set up 1 to TSS'
      TSSname <- seq(-tssbin * binbp, 0, by = everybp)
      TSSloc <- seq(1,  by = everybin, length.out = length(TSSname))
      # make sure TSS is included
      if (any(TSSname == 0)) {
        TSSloc[TSSname == 0] <- tssbin + .5
        TSSname[TSSname == 0] <- tssname
      } else if(any(TSSloc == tssbin)){
        TSSname[TSSloc == tssbin] <- tssname
        TSSloc[TSSloc == tssbin] <- tssbin + .5
      } else {
        TSSloc <- sort(c(TSSloc, tssbin + .5))
        TSSname <- append(TSSname, tssname)
      }
      my_5prim <-
        tibble(lloc = TSSloc, lname = as.character(TSSname))
      # keep distance between TSS and numbers
      nextStart <- tssbin + everybin
      if (body1bin > 0 & tesbin > 0 & body2bin > 0) {
        nextEnd <- body1bin
      } else if (tesbin == 0) {
        nextEnd <- totbins
      } else {
        nextEnd <- 0
      }
      if (nextStart <= nextEnd) {
        TSSname1 <- seq(everybp, (nextEnd - tssbin) * binbp, by = everybp)
        TSSloc1 <-
          seq(nextStart,
              by = everybin,
              length.out = length(TSSname1))
        # make sure body brake is included
        if (!any(TSSloc1 == nextEnd)) {
          TSSname1 <- append(TSSname1, (nextEnd - tssbin) * binbp)
          TSSloc1 <- c(TSSloc1, nextEnd)
        }
        my_5prim2 <-
          tibble(lloc = TSSloc1, lname = as.character(TSSname1))
        my_5prim <-
          full_join(my_5prim, my_5prim2, by = c("lloc", "lname"))
      } else if (nextEnd > 0) {
        my_5prim2 <-
          tibble(lloc = nextEnd, lname = as.character((nextEnd - tssbin) * binbp))
        my_5prim <-
          full_join(my_5prim, my_5prim2, by = c("lloc", "lname"))
      }
    }
    if (tesbin > 0) {
      # next to TES'
      if (body1bin > 0 & tssbin > 0 & body2bin > 0) {
        nextStart <- body2bin
      } else if (tssbin == 0) {
        nextStart <- 1
      } else {
        nextStart <- totbins
      }
      if (nextStart <= tesbin) {
        TESname <-  abs(seq((nextStart - tesbin) * binbp, 0, by = everybp))
        if(nextStart == 1){
          TESname <- TESname + binbp
        }
        TESloc <-
          seq(nextStart,
              by = everybin,
              length.out = length(TESname))
        # make sure TES is included
        if (any(TESname == 0)) {
          TESloc[TESname == 0] <- tesbin + .5
          TESname[TESname == 0] <- tesname
        } else if(any(TESloc == tesbin)){
          TESname[TESloc == tesbin] <- tesname
          TESloc[TESloc == tesbin] <- tesbin + .5
        }else {
          TESloc <- sort(c(TESloc, tesbin + .5))
          TESname <- append(TESname, tesname)
        }
        my_3prim <-
          tibble(lloc = TESloc, lname = as.character(TESname))
      } else {
        my_3prim <- tibble(lloc = tesbin + .5, lname = tesname)
      }
      # TES to end
      nextStart <- tesbin + everybin
      if (nextStart < totbins) {
        TESname1 <- seq(everybp, (totbins - tesbin) * binbp, by = everybp)
        TESloc1 <-
          seq(nextStart,
              by = everybin,
              length.out = length(TESname1))
        if (!any(TESloc1 == totbins)) {
          TESname1 <- append(TESname1, (totbins - tesbin) * binbp)
          TESloc1 <- c(TESloc1, totbins)
        }
        my_3prim2 <-
          tibble(lloc = TESloc1, lname = as.character(TESname1))
        my_3prim <-
          full_join(my_3prim, my_3prim2, by = c("lloc", "lname"))
      }
    }
    if (!is.null(my_5prim) & !is.null(my_3prim)) {
      my_53prim <-
        arrange(full_join(my_5prim, my_3prim, by = c("lloc", "lname")), lloc)
      use_plot_breaks <- my_53prim$lloc
      use_plot_breaks_labels <- my_53prim$lname
      use_plot_breaks_labels <-
        use_plot_breaks_labels[seq_along(use_plot_breaks)]
    } else if (!is.null(my_5prim)) {
      use_plot_breaks <- my_5prim$lloc
      use_plot_breaks_labels <- my_5prim$lname
      use_plot_breaks_labels <-
        use_plot_breaks_labels[seq_along(use_plot_breaks)]
    } else if (!is.null(my_3prim)) {
      use_plot_breaks <- my_3prim$lloc
      use_plot_breaks_labels <- my_3prim$lname
      use_plot_breaks_labels <-
        use_plot_breaks_labels[seq_along(use_plot_breaks)]
    } else {
      # just print bin numbers
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
      use_plot_breaks_labels <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
    }
    # no bp bin labels
  } else {
    if (everybin > 0) {
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
      use_plot_breaks_labels <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
    } else {
      use_plot_breaks <- .5
      use_plot_breaks_labels <- "none"
    }
  }
  # virtical line set up
  use_plot_breaks <- na_if(use_plot_breaks, 0.5)
  use_plot_breaks_labels <-
    use_plot_breaks_labels[!is.na(use_plot_breaks)]
  use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
  list(mybrakes = use_plot_breaks,
       mylabels = use_plot_breaks_labels)
}

# Sets plot lines and labels colors
LinesLabelsListPlot <-
  function(body1bin,
           body1color,
           body1line,
           body2bin,
           body2color,
           body2line,
           tssbin,
           tsscolor,
           tssline,
           tesbin,
           tescolor,
           tesline,
           use_plot_breaks_labels,
           use_plot_breaks,
           vlinesize,
           linesize,
           fontsizex,
           fontsizey,
           legendsize,
           ttestlinesize) {
    print("lines and labels plot fun")
    if (length(use_plot_breaks_labels) > 0) {
      mycolors <- rep("black", length(use_plot_breaks))
      use_virtical_line <- c(NA, NA, NA, NA)
      if (tssbin > 0) {
        mycolors[which(use_plot_breaks == tssbin  + .5)] <- tsscolor
        use_virtical_line[1] <- tssbin  + .5
        if (tssbin < body1bin &
            body1bin < body2bin &
            body2bin < tesbin & tesbin <= last(use_plot_breaks)) {
          use_virtical_line[3:4] <- c(body1bin, body2bin)
        }
      }
      if (tesbin > 0) {
        mycolors[which(use_plot_breaks == tesbin  + .5)] <- tescolor
        use_virtical_line[2] <- tesbin + .5
      }
    } else {
      use_plot_breaks <- .5
      use_plot_breaks_labels <- "none"
      use_virtical_line <- c(NA, NA, NA, NA)
    }
    # virtical line set up
    use_virtical_line_color <-
      c(tsscolor, tescolor, body1color, body2color)
    use_virtical_line_type <-
      c(tssline, tesline, body1line, body2line)
    use_plot_breaks <- na_if(use_plot_breaks, 0.5)
    use_virtical_line <- na_if(use_virtical_line, 0.5)
    use_plot_breaks_labels <-
      use_plot_breaks_labels[!is.na(use_plot_breaks)]
    use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
    use_virtical_line_type <-
      use_virtical_line_type[!is.na(use_virtical_line)]
    use_virtical_line_color <-
      use_virtical_line_color[!is.na(use_virtical_line)]
    use_virtical_line <-
      use_virtical_line[!is.na(use_virtical_line)]
    list(
      myline = virtical_line_data_frame <- data.frame(
        use_virtical_line,
        use_virtical_line_type,
        use_virtical_line_color,
        stringsAsFactors = FALSE
      ),
      mycolors = mycolors,
      mybrakes = use_plot_breaks,
      mylabels = use_plot_breaks_labels,
      mysize = c(vlinesize, linesize, fontsizex, fontsizey, legendsize, ttestlinesize)
    )
  }

# lines and labels preset helper
LinesLabelsPreSet <- function(mytype) {
  # 5|4, 4|3, tss, pA, bp/bin, max bins, every bin
  if (mytype == kLinesandlabels[1]) {
    tt <- c(20, 40, 15, 45, 100, LIST_DATA$x_plot_range[2], 10)
  } else if (mytype == kLinesandlabels[2]) {
    tt <- c(0, 0, 40, 0, 25, LIST_DATA$x_plot_range[2], 20)
  } else if (mytype == kLinesandlabels[3]) {
    tt <- c(0, 0, 0, 10, 100, LIST_DATA$x_plot_range[2], 10)
  } else if (mytype == kLinesandlabels[4]) {
    tt <- c(0, 0, 5, 0, 50, LIST_DATA$x_plot_range[2], 10)
  } else if (mytype == kLinesandlabels[5]) {
    tt <- c(0,
            0,
            floor(LIST_DATA$x_plot_range[2] * .25),
            0,
            100,
            LIST_DATA$x_plot_range[2],
            10)
  } else if (mytype == kLinesandlabels[6]) {
    tt <- c(
      0,
      0,
      floor(LIST_DATA$x_plot_range[2] * .25),
      ceiling(LIST_DATA$x_plot_range[2] * .75),
      100,
      LIST_DATA$x_plot_range[2],
      10
    )
  } else if (mytype == kLinesandlabels[7]) {
    tt <- c(0,
            0,
            0,
            ceiling(LIST_DATA$x_plot_range[2] * .75),
            100,
            LIST_DATA$x_plot_range[2],
            10)
  } else {
    tt <-
      c(
        ceiling(LIST_DATA$x_plot_range[2] * .33),
        floor(LIST_DATA$x_plot_range[2] * .66),
        floor(LIST_DATA$x_plot_range[2] * .25),
        ceiling(LIST_DATA$x_plot_range[2] * .75),
        100,
        LIST_DATA$x_plot_range[2],
        ceiling(LIST_DATA$x_plot_range[2] * .1)
      )
  }
  tt
}

# help get min and max from apply math data set
MyXSetValues <-
  function(apply_math,
           xBinRange,
           yBinRange = c(0, 100),
           log_2 = F) {
    tt <- group_by(apply_math, set) %>%
      dplyr::filter(bin %in% xBinRange[1]:xBinRange[2]) %>%
      ungroup() %>%
      summarise(min(value, na.rm = T), max(value, na.rm = T),.groups="drop") %>%
      unlist(., use.names = FALSE)
    tt <-
      c(tt[1] + (tt[1] * (yBinRange[1] / 100)), tt[2] + (tt[2] * ((yBinRange[2] -
                                                                     100) / 100)))
    if (log_2) {
      tt <- log2((abs(tt))^(sign(tt)))
    }
    tt
  }

# main ggplot function
GGplotLineDot <-
  function(list_long_data_frame,
           xBinRange,
           plot_options,
           yBinRange,
           line_list,
           use_smooth,
           plot_ttest,
           use_log2,
           use_y_label,
           plot_occupancy) {
    print("ggplot")
    legend_space <- lengths(strsplit(sort(plot_options$set), "\n")) / 1.1
    if (use_log2) {
      gp <-
        ggplot(
          list_long_data_frame,
          aes(
            x = as.numeric(bin),
            y = log2(value),
            color = set,
            shape = set,
            size = set,
            linetype = set
          )
        )
    } else {
      gp <-
        ggplot(
          list_long_data_frame,
          aes(
            x = as.numeric(bin),
            y = value,
            color = set,
            shape = set,
            size = set,
            linetype = set
          )
        )
    }
    if (use_smooth) {
      gp <- gp +
        geom_smooth(se = FALSE,
                    size = line_list$mysize[2],
                    span = .2) 
    } else{
      gp <- gp +
        geom_line(size = line_list$mysize[2],alpha=0.8)
    }
    gp <- gp +
      geom_point(stroke = .001) +
      scale_size_manual(values = plot_options$mysizedot) +
      scale_color_manual(values = plot_options$mycol) +
      scale_shape_manual(values = plot_options$mydot) +
      scale_linetype_manual(values = plot_options$myline) +
      ylab(use_y_label) +
      geom_vline(
        data = line_list$myline,
        aes(xintercept = use_virtical_line),
        size = line_list$mysize[1],
        linetype = line_list$myline$use_virtical_line_type,
        color = line_list$myline$use_virtical_line_color
      ) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      theme(axis.title.y = element_text(size =  line_list$mysize[4] + 4, margin = margin(2, 10, 2, 2))) +
      theme(axis.text.y = element_text(size = line_list$mysize[4],
                                       face = 'bold')) 
    if(!is_empty(LIST_DATA$ttest)){
      use_col_tt <- plot_ttest$options_main_tt$mycol
      use_line_tt <- plot_ttest$options_main_tt$myline
      names(use_col_tt) <- plot_ttest$options_main_tt$set
      names(use_line_tt) <- plot_ttest$options_main_tt$set
      gp2 <- ggplot(LIST_DATA$ttest, aes(y=p.value,x=bin,
                                         color=set,
                                         shape = set,
                                         size = set,
                                         linetype = set)) + 
        geom_line(size = line_list$mysize[6],alpha=0.8) +
        scale_color_manual(values = use_col_tt) +
        scale_linetype_manual(values = use_line_tt)+
        geom_hline(yintercept = plot_ttest$hlineTT,color="blue") + 
        theme_bw() +
        geom_vline(
          data = line_list$myline,
          aes(xintercept = use_virtical_line),
          size = line_list$mysize[1],
          linetype = line_list$myline$use_virtical_line_type,
          color = line_list$myline$use_virtical_line_color
        ) +
        xlab(paste(Sys.Date(), paste(unique(
          plot_options$sub
        ), collapse = ", "), collapse = ", ")) +
        ylab(plot_ttest$ylabTT)+
        scale_x_continuous(breaks = line_list$mybrakes[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
                           labels = line_list$mylabels[between(line_list$mybrakes, xBinRange[1], xBinRange[2])]) +
        theme(axis.title.x = element_text(size =  line_list$mysize[3], vjust = .5)) +
        theme(axis.title.y = element_text(size =  line_list$mysize[4], margin = margin(2, 10, 2, 2))) +
        theme(axis.text.y = element_text(size = line_list$mysize[4],
                                         face = 'bold'))+
        theme(
          axis.text.x = element_text(
            color = line_list$mycolors[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
            size = line_list$mysize[3],
            angle = -45,
            hjust = .1,
            vjust = .9,
            face = 'bold'
          )
        ) +
        theme(
          legend.title = element_blank(),
          legend.key = element_rect(size = line_list$mysize[5] / 2, color = 'white'),
          legend.key.height = unit(legend_space/1.2, "line"),
          legend.text = element_text(size = line_list$mysize[5]/1.2, face = 'bold')
        ) +
        coord_cartesian(xlim = xBinRange, ylim = plot_ttest$ylimTT)
      my_occupancy <- c(6 - plot_occupancy, plot_occupancy + .5)
      gp <- gp + coord_cartesian(xlim = xBinRange, ylim = unlist(yBinRange))
      suppressMessages(print(gp + gp2 + plot_layout(ncol = 1, heights = my_occupancy)))
      return(suppressMessages(gp+ gp2 + plot_layout(ncol = 1, heights = my_occupancy)))
    } else{
      gp <- gp + 
        xlab(paste(Sys.Date(), paste(unique(
          plot_options$sub
        ), collapse = ", "), collapse = ", ")) +
        scale_x_continuous(breaks = line_list$mybrakes[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
                           labels = line_list$mylabels[between(line_list$mybrakes, xBinRange[1], xBinRange[2])]) +
        theme(axis.title.x = element_text(size =  line_list$mysize[3], vjust = .5)) +
        theme(
          axis.text.x = element_text(
            #fix for coord_cartesian [between(line_list$mybrakes, xBinRange[1], xBinRange[2])]
            color = line_list$mycolors[between(line_list$mybrakes, xBinRange[1], xBinRange[2])],
            size = line_list$mysize[3],
            angle = -45,
            hjust = .1,
            vjust = .9,
            face = 'bold'
          )
        ) +
        theme(
          legend.title = element_blank(),
          legend.key = element_rect(size = line_list$mysize[5] / 2, color = 'white'),
          legend.key.height = unit(legend_space, "line"),
          legend.text = element_text(size = line_list$mysize[5], face = 'bold')
        )  +
        coord_cartesian(xlim = xBinRange, ylim = unlist(yBinRange))
      
      suppressMessages(print(gp))
      return(suppressMessages(gp))
    }
  }

# CDG ggplot function
GGplotC <-
  function(df2,
           plot_options,
           use_header) {
    print("ggplot CDF")
    use_col <- plot_options$mycol
    names(use_col) <- plot_options$set
    legend_space <- max(1, (lengths(strsplit(
      plot_options$set, "\n"
    ))))
    gp <- ggplot(df2, aes(log2(value), color = set)) +
      stat_ecdf(show.legend = TRUE, size = 1.8) +
      scale_color_manual(name = "Sample", values = use_col) +
      ylab("Fraction of genes") +
      ggtitle(use_header) +
      theme_bw() +
      theme(legend.title = element_blank()) +
      theme(axis.title.y = element_text(size =  15)) +
      theme(axis.title.x = element_text(size =  13, vjust = .5)) +
      theme(axis.text.x = element_text(
        size = 12,
        angle = -45,
        hjust = .1,
        vjust = .9,
        face = 'bold'
      )) +
      theme(
        legend.title = element_blank(),
        legend.key = element_rect(size = 5, color = 'white'),
        legend.key.height = unit(legend_space, "line"),
        legend.text = element_text(size = 10)
      )
    suppressMessages(print(gp))
    return(suppressMessages(gp))
  }

# server ----
server <- function(input, output, session) {
  # remove on non-local deployment
  session$onSessionEnded(stopApp)
  
  # reacive values ----
  reactive_values <- reactiveValues(
    pickerfile_controler = "",
    Y_Axis_Lable = NULL,
    Y_Axis_numbers = NULL,
    Lines_Labels_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL,
    Plot_controler_cluster = NULL,
    Plot_controler_ratio = NULL,
    Picker_controler = NULL,
    makeplot = 0,
    onoff = list(),
    binset = FALSE
  )
  
  # change tab controls ----
  observeEvent(input$tabs, ignoreInit = TRUE, {
    print("switching tabs")
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
    toggle("showcdftoolpicker",
           condition = (input$tabs == "cdftool" &
                          LIST_DATA$STATE[1] != 0))
    toggle("showgenelistspicker",
           condition = (input$tabs == "genelists" &
                          length(LIST_DATA$gene_file) > 0))
    # load tab
    if (input$tabs == "loaddata" & LIST_DATA$STATE[1] != 0) {
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_file),
        selected = names(LIST_DATA$gene_file)[1]
      )
    }
    # file norm tab
    if (input$tabs == "filenorm" & LIST_DATA$STATE[1] != 0) {
      updatePickerInput(
        session,
        "pickernumerator",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == names(LIST_DATA$gene_file)[1]),
          mycol)$mycol, sep = ":"))
      )
      updatePickerInput(
        session,
        "pickerdenominator",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == names(LIST_DATA$gene_file)[1]),
          mycol)$mycol, sep = ":"))
      )
      output$valueboxnormfile <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("cogs"),
                 color = "yellow")
      })
    }
    if (input$tabs == "genelists" &
        length(LIST_DATA$gene_file) != 0) {
      shinyjs::hide('actiongenelistsdatatable')
      if (length(LIST_DATA$gene_file) >= 1) {
        og <- input$pickergenelists
        if (!all(NULL %in% names(LIST_DATA$gene_file))) {
          og <- NULL
        }
        updatePickerInput(
          session,
          "pickergenelists",
          choices = names(LIST_DATA$gene_file),
          selected = og,
          choicesOpt = list(style = rep("color:black", length(
            names(LIST_DATA$gene_file)
          )))
        )
      }
      
    }
    if (input$tabs == "sorttool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectsortfile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- names(LIST_DATA$gene_file)[1]
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
        choices = distinct(LIST_DATA$gene_info, set)$set,
        selected = NULL,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == input$selectsortfile),
          mycol)$mycol, sep = ":"))
      )
      if (sum(grepl("^Filter", names(LIST_DATA$gene_file))) > 0) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$use$gene),
            "Gene List Filter",
            icon = icon("list"),
            color = "green"
          )
        })
      } else {
        shinyjs::hide('actionsortdatatable')
        output$valueboxsort <- renderValueBox({
          valueBox(0,
                   "Gene List Filter",
                   icon = icon("list"),
                   color = "green")
        })
      }
    }
    
    if (input$tabs == "ratiotool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectratiofile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectratiofile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      updatePickerInput(
        session,
        "pickerratio1file",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == input$selectratiofile),
          mycol)$mycol, sep = ":"))
      )
      updatePickerInput(
        session,
        "pickerratio2file",
        choices = c("None", distinct(LIST_DATA$gene_info, set)$set),
        selected = "None",
        choicesOpt = list(style = paste("color", c(
          "#000000",
          dplyr::select(
            dplyr::filter(LIST_DATA$gene_info,
                          gene_list == input$selectratiofile),
            mycol)$mycol), sep = ":"))
      )
      if (sum(grepl("Ratio_", names(LIST_DATA$gene_file))) == 0) {
        shinyjs::hide('actionratiodatatable')
      }
      
      if (any(grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio1 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file))]]$use$gene),
            "Ratio Up file1",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxratio1 <- renderValueBox({
          valueBox(0,
                   "Ratio Up file1",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (any(grep("Ratio_Down_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio2 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_Down_file1\nn =",
                                                 names(LIST_DATA$gene_file))]]$use$gene),
            "Ratio Up file2",
            icon = icon("list"),
            color = "blue"
          )
        })
      } else{
        output$valueboxratio2 <- renderValueBox({
          valueBox(0,
                   "Ratio Up file2",
                   icon = icon("list"),
                   color = "blue")
        })
      }
      if (any(grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio3 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file))]]$use$gene),
            "Ratio No Diff",
            icon = icon("list"),
            color = "yellow"
          )
        })
      } else{
        output$valueboxratio3 <- renderValueBox({
          valueBox(0,
                   "Ratio No Diff",
                   icon = icon("list"),
                   color = "yellow")
        })
      }
    }
    
    if (input$tabs == "clustertool" & LIST_DATA$STATE[1] != 0) {
      if(LIST_DATA$STATE[[2]] == 0){
        
      }
      ol <- input$selectclusterfile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- names(LIST_DATA$gene_file)[1]
      }
      og <- input$pickerclusterfile
      if (!og %in% distinct(LIST_DATA$gene_info, set)$set) {
        og <- NULL
      }
      updateSelectInput(
        session,
        "selectclusterfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      updatePickerInput(
        session,
        "pickerclusterfile",
        choices = distinct(LIST_DATA$gene_info, set)$set,
        selected = og,
        choicesOpt = list(style = paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == input$selectclusterfile),
          mycol)$mycol, sep = ":"))
      )
      if (!is.null(reactive_values$clustergroups) & any(grep(
        paste0(reactive_values$clustergroups, "1\nn ="),
        names(LIST_DATA$gene_file)
      ) > 0)) {
        output$valueboxcluster1 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep(
              paste0(reactive_values$clustergroups, "1\nn ="),
              names(LIST_DATA$gene_file)
            )]]$use$gene),
            "Cluster 1",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        shinyjs::hide("plotcluster")
        shinyjs::hide('actionclusterdatatable')
        output$valueboxcluster1 <- renderValueBox({
          valueBox(0,
                   "Cluster 1",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (!is.null(reactive_values$clustergroups) & any(grep(
        paste0(reactive_values$clustergroups, "2\nn ="),
        names(LIST_DATA$gene_file)
      ) > 0)) {
        output$valueboxcluster2 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep(
              paste0(reactive_values$clustergroups, "2\nn ="),
              names(LIST_DATA$gene_file)
            )]]$use$gene),
            "Cluster 2",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxcluster2 <- renderValueBox({
          valueBox(0,
                   "Cluster 2",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (!is.null(reactive_values$clustergroups) & any(grep(
        paste0(reactive_values$clustergroups, "3\nn ="),
        names(LIST_DATA$gene_file)
      ) > 0)) {
        output$valueboxcluster3 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep(
              paste0(reactive_values$clustergroups, "3\nn ="),
              names(LIST_DATA$gene_file)
            )]]$use$gene),
            "Gene List 3",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxcluster3 <- renderValueBox({
          valueBox(0,
                   "Cluster 3",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (!is.null(reactive_values$clustergroups) & any(grep(
        paste0(reactive_values$clustergroups, "4\nn ="),
        names(LIST_DATA$gene_file)
      ) > 0)) {
        output$valueboxcluster4 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep(
              paste0(reactive_values$clustergroups, "4\nn ="),
              names(LIST_DATA$gene_file)
            )]]$use$gene),
            "Cluster 4",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxcluster4 <- renderValueBox({
          valueBox(0,
                   "Gene List 4",
                   icon = icon("list"),
                   color = "green")
        })
      }
    }
    #CDF
    if (input$tabs == "cdftool" & LIST_DATA$STATE[1] != 0) {
      #update cdf dynamic picker
      pickercdf <- list()
      for (i in names(LIST_DATA$gene_file)[grep("CDF ", names(LIST_DATA$gene_file), invert = T)]) {
        pickercdf[[i]] <-
          list(div(
            style = "margin-bottom: -20px;",
            pickerInput(
              inputId = gsub(" ", "-cdfspace2-", gsub("\n", "-cdfspace1-", i)),
              label = i,
              width = "99%",
              choices = distinct(LIST_DATA$gene_info, set)$set,
              multiple = T,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 0"
              ),
              choicesOpt = list(style = paste("color", dplyr::select(
                dplyr::filter(LIST_DATA$gene_info,
                              gene_list == names(LIST_DATA$gene_file)[1]),
                mycol)$mycol, sep = ":"))
            )
          ))
      }
      output$DynamicCDFPicker_main <- renderUI({
        pickercdf[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Group_|^Cluster_")]
      })
      if(any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        output$DynamicCDFPicker_filter <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Filter")]
        })
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        output$DynamicCDFPicker_comparisons <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Gene_List_")]
        })
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        output$DynamicCDFPicker_ratio <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Ratio_")]
        })
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_"))){
        output$DynamicCDFPicker_clusters <- renderUI({
          pickercdf[str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_")]
        })
      }
      if (sum(grepl("CDF ", names(LIST_DATA$gene_file))) == 0) {
        output$plotcdf <- renderPlot({
          NULL
        })
        shinyjs::hide('plotcdf')
        shinyjs::hide('actioncdfdatatable')
        my_count <- 0
      } else {
        my_count <-
          n_distinct(LIST_DATA$gene_file[[grep("CDF ", names(LIST_DATA$gene_file))]]$use$gene)
      }
      output$valueboxcdf <- renderValueBox({
        valueBox(my_count,
                 "Gene List",
                 icon = icon("list"),
                 color = "green")
      })
    }
    # first time switch tab auto plot
    if (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0) {
      reactive_values$Picker_controler <- 
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$table_file, set)$set)
      if (LIST_DATA$STATE[2] < 1) {
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...',
                     value = 0,
                     {
                       list_data_frame <- Active_list_data(LIST_DATA)
                       if (!is_empty(list_data_frame)) {
                         LIST_DATA$gene_info <<- list_data_frame %>% 
                           distinct(set,gene_list,plot_set) %>%
                           full_join(LIST_DATA$gene_info,.,by=c("set","gene_list")) %>% 
                           dplyr::filter(!is.na(set)) %>% 
                           dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                           dplyr::select(-plot_set.y,-plot_set.x) %>%
                           distinct()
                         reactive_values$Apply_Math <-
                           ApplyMath(
                             list_data_frame,
                             input$myMath,
                             input$selectplotnrom,
                             as.numeric(input$selectplotBinNorm)
                           )
                         
                         reactive_values$Y_Axis_Lable <- YAxisLable()
                         reactive_values$Plot_Options <-
                           MakePlotOptionFrame(LIST_DATA)
                       }
                     })
      } else if (LIST_DATA$STATE[2] == 2) {
        shinyjs::show("actionmyplotshow")
        reactive_values$Apply_Math <- NULL
      }
    }
  })
  
  # sets initual states ----
  observeEvent(reactive_values$binset, {
    #Plot
    updateSliderInput(
      session,
      "sliderplotBinRange",
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = LIST_DATA$x_plot_range
    )
    updateSelectInput(
      session,
      "selectplotBinNorm",
      choices = c(0:LIST_DATA$x_plot_range[2]),
      selected = 0
    )
    
    #Filter
    updateSliderInput(
      session,
      "slidersortbinrange",
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = LIST_DATA$x_plot_range
    )
    
    #Ratio
    updateSliderInput(
      session,
      "sliderbinratio1",
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = c(
        floor(LIST_DATA$x_plot_range[2] / 5.5),
        floor(LIST_DATA$x_plot_range[2] / 4.4)
      )
    )
    updateSliderInput(
      session,
      "sliderbinratio2",
      min = 0,
      max = LIST_DATA$x_plot_range[2],
      value = c(
        floor(LIST_DATA$x_plot_range[2] / 4.4) + 1,
        floor(LIST_DATA$x_plot_range[2] / 1.77)
      )
    )
    updateSliderInput(
      session,
      "sliderRatioBinNorm",
      min = 0,
      max = LIST_DATA$x_plot_range[2],
      value = 0
    )
    
    #Cluster
    updateSliderInput(
      session,
      "sliderbincluster",
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = LIST_DATA$x_plot_range
    )
    updateSelectInput(
      session,
      "selectplotBinNormcluster",
      choices = c(0:LIST_DATA$x_plot_range[2]),
      selected = 0
    )
    #CDF
    updateSliderInput(
      session,
      "sliderbincdf1",
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = c(
        floor(LIST_DATA$x_plot_range[2] / 5.5),
        floor(LIST_DATA$x_plot_range[2] / 4.4)
      )
    )
    updateSliderInput(
      session,
      "sliderbincdf2",
      min = LIST_DATA$x_plot_range[1],
      max = LIST_DATA$x_plot_range[2],
      value = c(
        floor(LIST_DATA$x_plot_range[2] / 4.4) + 1,
        floor(LIST_DATA$x_plot_range[2] / 1.77)
      )
    )
  })
  
  # loads data file(s) ----
  observeEvent(input$filetable, {
    print("load file")
    shinyjs::disable("startoff")
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
    
    if (LIST_DATA$STATE[1] == 0) {
      shinyjs::show("filegene1")
      shinyjs::show("downloadGeneList")
      shinyjs::show("filecolor")
      shinyjs::show("hiddensave")
      shinyjs::show("startoff")
      print("1st slider and plot lines Ylable")
      reactive_values$binset <- TRUE
      shinyjs::removeClass(selector = "body", class = "sidebar-collapse")
      # tries to guess lines and labels type
      num_bins <- LIST_DATA$x_plot_range[2]
      if (num_bins == 80 & LIST_DATA$STATE[3] == '543') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[1])
      } else if (num_bins == 80 & LIST_DATA$STATE[3] == '5') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[2])
      } else if (num_bins <= 60 & LIST_DATA$STATE[3] == '543') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[3])
      } else if (num_bins == 205 & LIST_DATA$STATE[3] == '5') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[4])
      } else if (LIST_DATA$STATE[3] == '5') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[5])
      } else if (LIST_DATA$STATE[3] == '4') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[6])
      } else if (LIST_DATA$STATE[3] == '3') {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[7])
      } else {
        updateSelectInput(session, "selectlineslabels", selected = kLinesandlabels[8])
      }
      LIST_DATA$STATE[1] <<- 1
    }
    
    shinyjs::enable("startoff")
    shinyjs::reset("filetable")
    ff <- distinct(LIST_DATA$table_file, set)$set
    updateSelectInput(session,
                      "selectdataoption",
                      choices = ff)
    
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
  
  # loads color file ----
  observeEvent(input$filecolor, {
    my_sel <- input$selectdataoption
    my_list <- input$selectgenelistoptions
    print("load color file")
    # load info, update select boxes, switching works and chaning info and ploting
    LIST_DATA <<- LoadColorFile(input$filecolor$datapath,
                                LIST_DATA, my_list)
    updateColourInput(session, "colourhex", value = paste(dplyr::filter(LIST_DATA$gene_info,gene_list == my_list & set == my_sel) %>% 
                                                            dplyr::select(mycol)))
    if (LIST_DATA$STATE[1] != 0 &
        !is.null(reactive_values$Apply_Math) &
        LIST_DATA$STATE[2] != 2) {
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options,
          reactive_values$Y_Axis_numbers,
          reactive_values$Lines_Labels_List,
          input$checkboxsmooth, reactive_values$Plot_Options_ttest,
          input$checkboxlog2,
          reactive_values$Y_Axis_Lable,
          input$sliderplotOccupancy
        )
    }
    shinyjs::reset("filecolor")
  })
  
  # update desplay selected item info ----
  observeEvent(c(input$selectdataoption, input$selectgenelistoptions),
               ignoreInit = TRUE,
               {
                 if (LIST_DATA$STATE[1] == 0) {
                   return()
                 }
                 my_sel <- LIST_DATA$gene_info %>% 
                   dplyr::filter(gene_list == input$selectgenelistoptions & 
                                   set == input$selectdataoption)
                 print("options update")
                 updateColourInput(session, "colourhex", value = paste(my_sel$mycol))
                 updateTextInput(session,
                                 "textnickname",
                                 value = paste(my_sel$set))
                 updateSelectInput(session,
                                   "selectdot",
                                   selected = paste(my_sel$mydot))
                 updateSelectInput(session,
                                   "selectline",
                                   selected = paste(my_sel$myline))
                 if (input$selectgenelistoptions == names(LIST_DATA$gene_file)[1]) {
                   shinyjs::disable("actionremovegene")
                   updateActionButton(session,"BttnNewColor",label = "Set all lists same color")
                   shinyjs::enable("textnickname")
                   shinyjs::enable("actionoptions")
                 } else {
                   shinyjs::enable("actionremovegene")
                   updateActionButton(session,"BttnNewColor",label = "Set color same as Compleat")
                   shinyjs::disable("textnickname")
                   shinyjs::disable("actionoptions")
                 }
               })
  
  # save functions, gene list, color list ----
  output$downloadGeneList <- downloadHandler(
    filename = function() {
      if (input$radiogroupsave == "Save Compleat color - file pair") {
        paste(Sys.Date(), ".color.txt", sep = "")
      } else if (input$radiogroupsave == "Save full Table file"){
        paste(input$selectdataoption, ".table", sep = "")
      } else if(input$radiogroupsave == "Save Gene list as bed"){
        paste(gsub("\nn = ", " n = ", input$selectgenelistoptions),
              Sys.Date(),
              ".bed",
              sep = "_")
      } else {
        paste(gsub("\nn = ", " n = ", input$selectgenelistoptions),
              Sys.Date(),
              ".txt",
              sep = "_")
      }
    },
    content = function(file) {
      if ((input$radiogroupsave == "Save Compleat color - file pair")) {
        new_comments <- NULL
        for (i in distinct(LIST_DATA$gene_info, set)$set) {
          new_comments <-
            c(new_comments,
              paste(i, LIST_DATA$gene_info %>% 
                      dplyr::filter(gene_list == input$selectgenelistoptions &
                                      set == i) %>% dplyr::select(mycol),
                    sep = ","))
        }
        write_lines(new_comments, file)
      } else if (input$radiogroupsave == "Save full Table file"){
        new_comments <- LIST_DATA$table_file %>% dplyr::filter(set == input$selectdataoption) %>% 
          semi_join(., LIST_DATA$gene_file[[input$selectgenelistoptions]]$use)
        write_tsv(new_comments, file)
      } else if (input$radiogroupsave == "Save Gene list as bed") {
        new_comments <-
          gsub("(:|;)", "*",
               LIST_DATA$gene_file[[input$selectgenelistoptions]]$use$gene) %>% 
          sub("-","*",.) %>%  
          sub("(?=[-+])","*",.,perl=TRUE) %>% 
          str_split_fixed(.,"\\*",5) %>% as_tibble() %>% 
          dplyr::mutate(score=0) %>% dplyr::select(V1,V2,V3,V5,score,V4)
        write_tsv(new_comments,file,col_names = FALSE)
      } else {
        new_comments <- paste("#", Sys.Date(), "\n# File(s) used:")
        new_comments <-
          c(new_comments, paste("#", distinct(LIST_DATA$gene_info, set)$set))
        new_comments <-
          c(new_comments,  paste("\n#", gsub(
            "\nn = ", " n = ",  input$selectgenelistoptions
          )))
        new_comments <-
          c(new_comments, paste("#", gsub(
            "\nn = ", " n = ",
            paste(LIST_DATA$gene_file[[input$selectgenelistoptions]]$info)
          )))
        new_comments <-
          c(new_comments, LIST_DATA$gene_file[[input$selectgenelistoptions]]$use$gene)
        write_lines(new_comments, file)
      }
    }
  )
  
  # record new nickname  ----
  observeEvent(input$actionoptions, ignoreInit = TRUE, {
    # sets/resets nickname
    if (nchar(input$textnickname) == 0) {
      updateTextInput(session,
                      "textnickname",
                      value = paste(input$selectdataoption))
    } else {
      if (input$textnickname != input$selectdataoption) {
        print("new nickname")
        if (any(input$textnickname == distinct(LIST_DATA$gene_info, set)$set)) {
          updateTextInput(session,
                          "textnickname",
                          value = paste0(input$selectdataoption,"-",input$textnickname,
                                         "-dup"))
        }
        LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
          dplyr::mutate(set = if_else(set == input$selectdataoption,
                                      input$textnickname, set)) %>% 
          dplyr::mutate(onoff = if_else(onoff == input$selectdataoption,
                                        input$textnickname, onoff))
        LIST_DATA$table_file <<- LIST_DATA$table_file %>%
          dplyr::mutate(set = if_else(set == input$selectdataoption,
                                      input$textnickname, set))
        if (LIST_DATA$STATE[2] != 0) {
          LIST_DATA$STATE[2] <<- 2
        }
      }
    }
    ff <- distinct(LIST_DATA$table_file, set)$set
    updateSelectInput(session,
                      "selectdataoption",
                      choices = ff)
  })
  
  # records new dot options ----
  observeEvent(input$selectdot, {
    if (!is.null(names(LIST_DATA$gene_file))) {
      print(input$selectgenelistoptions)
      print(input$selectdataoption)
      if (input$selectdot != dplyr::filter(LIST_DATA$gene_info,
                                           gene_list == input$selectgenelistoptions & set == input$selectdataoption) %>% 
          dplyr::select(mydot)){
        print("new dot")
        LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
          dplyr::mutate(mydot = if_else(gene_list == input$selectgenelistoptions & set == input$selectdataoption,
                                        input$selectdot, mydot))
        if (LIST_DATA$STATE[1] != 0 &
            !is.null(reactive_values$Apply_Math) &
            LIST_DATA$STATE[2] != 2) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
              reactive_values$Plot_Options,
              reactive_values$Y_Axis_numbers,
              reactive_values$Lines_Labels_List,
              input$checkboxsmooth, reactive_values$Plot_Options_ttest,
              input$checkboxlog2,
              reactive_values$Y_Axis_Lable,
              input$sliderplotOccupancy
            )
        }
      }
    }
  })
  
  # records new line options ----
  observeEvent(input$selectline, {
    if (!is.null(names(LIST_DATA$gene_file))) {
      if (input$selectdot != dplyr::filter(LIST_DATA$gene_info,
                                           gene_list == input$selectgenelistoptions & set == input$selectdataoption) %>% 
          dplyr::select(myline)){
        print("new line")
        LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
          dplyr::mutate(myline = if_else(gene_list == input$selectgenelistoptions & set == input$selectdataoption,
                                         input$selectline, myline))
        if (LIST_DATA$STATE[1] != 0 &
            !is.null(reactive_values$Apply_Math) &
            LIST_DATA$STATE[2] != 2) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
              reactive_values$Plot_Options,
              reactive_values$Y_Axis_numbers,
              reactive_values$Lines_Labels_List,
              input$checkboxsmooth, reactive_values$Plot_Options_ttest,
              input$checkboxlog2,
              reactive_values$Y_Axis_Lable,
              input$sliderplotOccupancy
            )
        }
      }
    }
  })
  
  # update color based on rgb text input ----
  observeEvent(input$actionmyrgb, {
    print("color rgb")
    updateColourInput(session, "colourhex", value = RgbToHex(input$textrgbtohex, convert = "hex"))
  })
  
  # update and save color selected ----
  observeEvent(input$colourhex, ignoreInit = TRUE, {
    print("update text color")
    updateTextInput(session,
                    "textrgbtohex",
                    value = RgbToHex(x = input$colourhex, convert = "rgb"))
    if (!is.null(names(LIST_DATA$gene_file))) {
      my_sel <- LIST_DATA$gene_info %>% 
        dplyr::filter(gene_list == input$selectgenelistoptions & 
                        set == input$selectdataoption)
      if (input$colourhex != my_sel$mycol) {
        print("color new")
        LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
          dplyr::mutate(mycol=if_else(gene_list == input$selectgenelistoptions & 
                                        set == input$selectdataoption,
                                      input$colourhex, mycol))
        
        if (LIST_DATA$STATE[1] != 0 &
            !is.null(reactive_values$Apply_Math) &
            LIST_DATA$STATE[2] != 2) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          reactive_values$Plot_controler <-
            GGplotLineDot(
              reactive_values$Apply_Math,
              input$sliderplotBinRange,
              reactive_values$Plot_Options,
              reactive_values$Y_Axis_numbers,
              reactive_values$Lines_Labels_List,
              input$checkboxsmooth, reactive_values$Plot_Options_ttest,
              input$checkboxlog2,
              reactive_values$Y_Axis_Lable,
              input$sliderplotOccupancy
            )
        }
        my_sel <- LIST_DATA$gene_info %>% 
          dplyr::filter(gene_list == input$selectgenelistoptions) %>% 
          dplyr::select(mycol)
        reactive_values$Picker_controler <-
          paste("color", unlist(my_sel), sep = ":")
      }
    }
  })
  
  # reactive picker watcher (watches everything but graps only gene lists with white space fix) ----
  observeEvent(reactiveValuesToList(input)[gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", names(LIST_DATA$gene_file)))],
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 reactive_values$picker <-
                   reactiveValuesToList(input)[gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", names(LIST_DATA$gene_file)))]
               })
  
  # records check box on/off ----
  observeEvent(reactive_values$picker,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 print("checkbox on/off")
                 # needed for controling flow of first time auto plot
                 if (LIST_DATA$STATE[1] != 0) {
                   ttt <- reactive_values$picker
                   checkboxonoff <- list()
                   for(i in names(ttt)){
                     onoff_name <-
                       gsub("-bensspace2-", " ", gsub("-bensspace1-", "\n", i))
                     if(!is_empty(ttt[[i]])){
                       checkboxonoff <- bind_rows(checkboxonoff, tibble(gene_list = onoff_name,
                                                                        onoff = ttt[[i]], set = ttt[[i]]))
                     } else if(!all(is.na(names(ttt)))){
                       checkboxonoff <- bind_rows(checkboxonoff, tibble(gene_list = onoff_name,
                                                                        onoff = NA, set = NA))
                     }
                   }  
                   reactive_values$onoff <- checkboxonoff
                 }
               })
  
  
  # sets and resets plot button on/off ----
  observeEvent(reactive_values$onoff,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 print("toggle on/off")
                 print(LIST_DATA$STATE[2])
                 LIST_DATA$gene_info <<-
                   CheckBoxOnOff(reactive_values$onoff,
                                 LIST_DATA$gene_info)
                 if (LIST_DATA$STATE[2] > 0) {
                   shinyjs::show("actionmyplotshow")
                   shinyjs::disable("numericYRangeHigh")
                   shinyjs::disable("numericYRangeLow")
                   shinyjs::disable("numericYRangeHighpval")
                   shinyjs::disable("numericYRangeLowpval")
                   LIST_DATA$STATE[2] <<- 2
                 } else if(LIST_DATA$STATE[2] == -2){
                   LIST_DATA$STATE[2] <<- -1
                 } else {
                   LIST_DATA$STATE[2] <<- 1
                 }
               })
  
  # plots when action button is pressed ----
  observeEvent(input$actionmyplot, ignoreInit = TRUE, {
    print("plot button")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   list_data_frame <- Active_list_data(LIST_DATA)
                   if (!is_empty(list_data_frame)) {
                     LIST_DATA$gene_info <- list_data_frame %>% 
                       distinct(set,gene_list,plot_set) %>%
                       full_join(LIST_DATA$gene_info,.,by=c("set","gene_list")) %>% 
                       dplyr::filter(!is.na(set)) %>% 
                       dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                       dplyr::select(-plot_set.y,-plot_set.x) %>%
                       distinct()
                     reactive_values$Apply_Math <-
                       ApplyMath(
                         list_data_frame,
                         input$myMath,
                         input$selectplotnrom,
                         as.numeric(input$selectplotBinNorm)
                       )
                     mm <- 0
                     if(!is.null(LIST_DATA$ttest$set)){
                       mm <- round(extendrange(range(LIST_DATA$ttest$p.value,na.rm = T,finite=T),f = .1),digits = 2)
                       p_cutoff <- input$hlinettest
                       if(input$selectttestlog =="-log"){
                         p_cutoff <-  -log(input$hlinettest)
                       } else if(input$selectttestlog =="-log10"){
                         p_cutoff <-  -log10(input$hlinettest)
                       }
                       if(mm[1]>0){
                         mm[1] <- 0
                       }
                       if(mm[2]<p_cutoff){
                         mm[2] <- p_cutoff
                       }
                       updateNumericInput(session,
                                          "numericYRangeHighpval",
                                          value = round(max(mm), 4))
                       updateNumericInput(session,
                                          "numericYRangeLowpval",
                                          value = round(min(mm), 4))
                     }
                     if(input$switchttest!="none"){
                       updateSelectInput(session,"selectttestitem", choices = distinct(LIST_DATA$ttest,set)$set)
                     }else{
                       updateSelectInput(session,"selectttestitem", choices = "none",
                                         selected="none")
                     }
                     reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
                     LIST_DATA$STATE[2] <<- 1
                     if (!is.null(reactive_values$Apply_Math) & is.null(LIST_DATA$ttest$set)){
                       updateNumericInput(session,
                                          "numericYRangeHighpval",
                                          value = 0)
                       updateNumericInput(session,
                                          "numericYRangeLowpval",
                                          value = 0)
                     }
                   } else {
                     LIST_DATA$STATE[2] <<- 2
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
                 })
    shinyjs::hide("actionmyplotshow")
    shinyjs::enable("numericYRangeHigh")
    shinyjs::enable("numericYRangeLow")
    shinyjs::enable("numericYRangeHighpval")
    shinyjs::enable("numericYRangeLowpval")
  })
  
  # updates Apply_Math ----
  observeEvent(reactive_values$Apply_Math, {
    print("updates reactive_values$Apply_Math")
    reactive_values$Y_Axis_numbers <-
      MyXSetValues(
        reactive_values$Apply_Math,
        input$sliderplotBinRange,
        c(0,100),
        input$checkboxlog2
      )
    my_step <-
      (max(reactive_values$Y_Axis_numbers) - min(reactive_values$Y_Axis_numbers)) /
      20
    updateNumericInput(session,
                       "numericYRangeHigh",
                       value = round(max(reactive_values$Y_Axis_numbers), 4),
                       step = my_step)
    updateNumericInput(session,
                       "numericYRangeLow",
                       value = round(min(reactive_values$Y_Axis_numbers), 4),
                       step = my_step)
    reactive_values$makeplot <- reactive_values$makeplot + 1
  })
  
  # renders plot ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  output$plot1sort <- renderPlot({
    reactive_values$Plot_controler_sort_min
  })
  output$plot2sort <- renderPlot({
    reactive_values$Plot_controler_sort_max
  })
  output$plotcluster <- renderPlot({
    reactive_values$Plot_controler_cluster
  })
  output$plotratio <- renderPlot({
    reactive_values$Plot_controler_ratio
  })
  # updates norm applymath ----
  observeEvent(c(input$myMath,
                 input$selectplotBinNorm,
                 input$selectplotnrom),
               ignoreInit = TRUE,
               {
                 reactive_values$Y_Axis_Lable <-
                   YAxisLable(
                     input$myMath,
                     input$selectplotnrom,
                     as.numeric(input$selectplotBinNorm),
                     input$checkboxsmooth,
                     input$checkboxlog2
                   )
                 if (LIST_DATA$STATE[1] != 0 &
                     LIST_DATA$STATE[2] != 2) {
                   print("apply math")
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  list_data_frame <- Active_list_data(LIST_DATA)
                                  if (!is_empty(list_data_frame)) {
                                    LIST_DATA$gene_info <- list_data_frame %>% 
                                      distinct(set,gene_list,plot_set) %>%
                                      full_join(LIST_DATA$gene_info,.,by=c("set","gene_list")) %>% 
                                      dplyr::filter(!is.na(set)) %>% 
                                      dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                      dplyr::select(-plot_set.y,-plot_set.x) %>%
                                      distinct()
                                    reactive_values$Apply_Math <-
                                      ApplyMath(
                                        list_data_frame,
                                        input$myMath,
                                        input$selectplotnrom,
                                        as.numeric(input$selectplotBinNorm)
                                      )
                                  }
                                })
                 }
               })
  
  # occupancy slider t.test is trigger ----
  observeEvent(input$sliderplotOccupancy, ignoreInit = T, {
    print("y slider")
    if (!is.null(reactive_values$Apply_Math)& input$switchttest != "none") {
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      if(input$switchttest!="none"){
        updateSelectInput(session,"selectttestitem", choices = distinct(LIST_DATA$ttest,set)$set)
      }else{
        updateSelectInput(session,"selectttestitem", choices = "none",
                          selected="none")
      }
    }
    if (!is.null(reactive_values$Apply_Math)) {
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options,
          reactive_values$Y_Axis_numbers,
          reactive_values$Lines_Labels_List,
          input$checkboxsmooth, reactive_values$Plot_Options_ttest,
          input$checkboxlog2,
          reactive_values$Y_Axis_Lable,
          input$sliderplotOccupancy
        )
    }
  })
  
  # t.test select file ----
  observeEvent(input$selectttestitem, ignoreInit = T, {
    print("t.test select")
    if (input$selectttestitem != "none" & !is_empty(LIST_DATA$ttest)) {
      myline <- LIST_DATA$ttest %>% dplyr::filter(set == input$selectttestitem) %>% distinct(myline)
      mycol <- LIST_DATA$ttest %>% dplyr::filter(set == input$selectttestitem) %>% distinct(mycol)
      updateSelectInput(session, "selectlinettest", selected = myline)
      updateColourInput(session, "selectcolorttest", value = paste(mycol))
    }
  })
  
  # t.test select plot options change ----
  observeEvent(input$actionttest, ignoreInit = T, {
    print("t.test select")
    if (LIST_DATA$STATE[2] !=2 & 
        !is.null(LIST_DATA$ttest$set) &
        input$selectttestitem != "none") {
      LIST_DATA$ttest <<- LIST_DATA$ttest %>% 
        dplyr::mutate(myline=ifelse(set == input$selectttestitem,input$selectlinettest,myline))
      LIST_DATA$ttest <<- LIST_DATA$ttest %>% 
        dplyr::mutate(mycol=ifelse(set == input$selectttestitem,input$selectcolorttest,mycol))
      reactive_values$Plot_Options_ttest <- MakePlotOptionttest(LIST_DATA$ttest,
                                                                c(round(min(input$numericYRangeLowpval), 4),round(max(input$numericYRangeHighpval), 4)),
                                                                input$selectttestlog,input$hlinettest,input$padjust,input$switchttesttype)
      if(reactive_values$Lines_Labels_List$mysize[6] == as.numeric(input$selectttestlinesize)){
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options,
            reactive_values$Y_Axis_numbers,
            reactive_values$Lines_Labels_List,
            input$checkboxsmooth, reactive_values$Plot_Options_ttest,
            input$checkboxlog2,
            reactive_values$Y_Axis_Lable,
            input$sliderplotOccupancy
          )
      }else{
        reactive_values$Lines_Labels_List$mysize[6] <- as.numeric(input$selectttestlinesize)
      }
    }
  })
  
  # t.test my_math ----
  observeEvent(c(reactive_values$Apply_Math,
                 input$selectttestlog,
                 input$switchttesttype,input$padjust,
                 input$selectttestalt,
                 input$selectttestexact,
                 input$selectttestpaired,
                 input$sliderplotOccupancy,
                 input$switchttest), ignoreInit = T, {
                   if (LIST_DATA$STATE[2] != 2 & 
                       input$switchttest != "none"){
                     print("t.test select")
                     withProgress(message = 'Calculation in progress',
                                  detail = 'This may take a while...',
                                  value = 0,
                                  {
                                    list_data_frame <- Active_list_data(LIST_DATA)
                                    if (!is_empty(list_data_frame)) {
                                      LIST_DATA$gene_info <- list_data_frame %>% 
                                        distinct(set,gene_list,plot_set) %>%
                                        full_join(LIST_DATA$gene_info,.,by=c("set","gene_list")) %>% 
                                        dplyr::filter(!is.na(set)) %>% 
                                        dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                        dplyr::select(-plot_set.y,-plot_set.x) %>%
                                        distinct()
                                      ttest <-
                                        ApplyTtest(list_data_frame,
                                                   input$switchttest,
                                                   input$selectttestlog,
                                                   input$switchttesttype,
                                                   input$padjust,
                                                   input$selectttestalt,
                                                   input$selectttestexact,
                                                   input$selectttestpaired)
                                    }
                                  })
                     LIST_DATA$ttest <<- ttest
                     mm <- 0
                     if (!is_empty(list_data_frame)) {
                       mm <- round(extendrange(range(ttest$p.value, na.rm = T,finite=T),f = .1),digits = 2)
                       p_cutoff <- input$hlinettest
                       if(input$selectttestlog == "-log"){
                         p_cutoff <- -log(input$hlinettest)
                       } else if(input$selectttestlog == "-log10"){
                         p_cutoff <- -log10(input$hlinettest)
                       }
                       if(mm[1] > 0){
                         mm[1] <- 0
                       }
                       if(mm[2] < p_cutoff){
                         mm[2] <- p_cutoff
                       }
                       updateNumericInput(session,
                                          "numericYRangeHighpval",
                                          value = round(max(mm), 4))
                       updateNumericInput(session,
                                          "numericYRangeLowpval",
                                          value = round(min(mm), 4))
                       
                       if(input$switchttest!="none"){
                         updateSelectInput(session,"selectttestitem", choices = distinct(ttest,set)$set)
                       }else{
                         updateSelectInput(session,"selectttestitem", choices = "none",
                                           selected="none")
                       }
                       if (is.null(ttest$set)){
                         updateNumericInput(session,
                                            "numericYRangeHighpval",
                                            value = 0)
                         updateNumericInput(session,
                                            "numericYRangeLowpval",
                                            value = 0)
                         reactive_values$Plot_Options_ttest <- NULL
                       }
                     } 
                     if(!is_empty(list_data_frame) & 
                        reactive_values$Lines_Labels_List$mysize[6] == as.numeric(input$selectttestlinesize)){
                       reactive_values$Plot_Options_ttest <<- MakePlotOptionttest(ttest,c(round(min(mm), 4),round(max(mm), 4)),
                                                                                  input$selectttestlog,input$hlinettest,input$padjust,input$switchttesttype)
                       reactive_values$Plot_controler <-
                         GGplotLineDot(
                           reactive_values$Apply_Math,
                           input$sliderplotBinRange,
                           reactive_values$Plot_Options,
                           reactive_values$Y_Axis_numbers,
                           reactive_values$Lines_Labels_List,
                           input$checkboxsmooth, reactive_values$Plot_Options_ttest,
                           input$checkboxlog2,
                           reactive_values$Y_Axis_Lable,
                           input$sliderplotOccupancy
                         )
                     }else{
                       reactive_values$Lines_Labels_List$mysize[6] <- as.numeric(input$selectttestlinesize)
                     }
                   }
                 })
  
  # updates plot actionButtonYXrange ----
  observeEvent(input$actionButtonYXrange, ignoreInit = T, {
    print("actionButtonYXrange")
    if (!is.null(reactive_values$Apply_Math) &
        input$actionButtonYXrange & LIST_DATA$STATE[2] != 2) {
      reactive_values$Y_Axis_numbers <-
        c(input$numericYRangeLow,input$numericYRangeHigh)
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options,
          reactive_values$Y_Axis_numbers,
          reactive_values$Lines_Labels_List,
          input$checkboxsmooth, reactive_values$Plot_Options_ttest,
          input$checkboxlog2,
          reactive_values$Y_Axis_Lable,
          input$sliderplotOccupancy
        )
    }
    if (LIST_DATA$STATE[2] == 2) {
      my_step <-
        (max(reactive_values$Y_Axis_numbers) - min(reactive_values$Y_Axis_numbers)) /
        20
      updateNumericInput(session,
                         "numericYRangeHigh",
                         value = round(max(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
      updateNumericInput(session,
                         "numericYRangeLow",
                         value = round(min(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
    }
  })
  
  # plots when bin slider or other triggers is triggered ----
  observeEvent(reactive_values$makeplot,
    ignoreInit = TRUE,
    {
      if (!is.null(reactive_values$Apply_Math) &
          LIST_DATA$STATE[2] != 2) {
        print("reactive_values$makeplot plots")
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options,
            reactive_values$Y_Axis_numbers,
            reactive_values$Lines_Labels_List,
            input$checkboxsmooth, reactive_values$Plot_Options_ttest,
            input$checkboxlog2,
            reactive_values$Y_Axis_Lable,
            input$sliderplotOccupancy
          )
      }
    }
  )
  
  # quick lines and labels preset change ----
  observeEvent(input$selectlineslabels, ignoreInit = TRUE, {
    if (input$selectlineslabels == "") {
      return()
    }
    print("quick Lines & Labels")
    myset <- LinesLabelsPreSet(input$selectlineslabels)
    updateNumericInput(session, "numericbody1", value = myset[1])
    updateNumericInput(session, "numericbody2", value = myset[2])
    updateNumericInput(session, "numerictss", value = myset[3])
    updateNumericInput(session, "numerictes", value = myset[4])
    updateNumericInput(session, "numericbinsize", value = myset[5])
    updateNumericInput(session, "numericlabelspaceing", value = myset[7])
    
  })
  
  # keep sizes real numbers ---
  observeEvent(
    c(
      input$selectvlinesize,
      input$selectlinesize,
      input$selectfontsizex,
      input$selectfontsizey,
      input$selectlegendsize
    ),
    ignoreInit = TRUE,
    {
      mynum <- c(2, 2.5, 13, 13, 10)
      myset <- c(
        input$selectvlinesize,
        input$selectlinesize,
        input$selectfontsizex,
        input$selectfontsizey,
        input$selectlegendsize
      )
      # keep bin positions in bounds > 0
      for (i in seq_along(myset)) {
        if (is.na(myset[i]) | myset[i] < 0) {
          myset[i] <- mynum[i]
          updateNumericInput(session, "selectvlinesize", value = myset[1])
          updateNumericInput(session, "selectlinesize", value = myset[2])
          updateNumericInput(session, "selectfontsizex", value = myset[3])
          updateNumericInput(session, "selectfontsizey", value = myset[4])
          updateNumericInput(session, "selectlegendsize", value = myset[5])
        }
      }
    }
  )
  
  # Update lines and labels ----
  observeEvent(
    c(
      input$numericbody1,
      input$numericbody2,
      input$numerictss,
      input$numerictssname,
      input$numerictes,
      input$numerictesname,
      input$numericbinsize,
      input$numericlabelspaceing
    ),
    ignoreInit = TRUE,
    {
      print("observe line and labels")
      myset <- c(
        input$numericbody1,
        input$numericbody2,
        input$numerictss,
        input$numerictes,
        input$numericbinsize,
        input$numericlabelspaceing
      )
      # keep bin positions in bounds > 0
      for (i in seq_along(myset)) {
        if (is.na(myset[i]) | myset[i] < 0) {
          myset[i] <- 0
          updateNumericInput(session, "numericbody1", value = myset[1])
          updateNumericInput(session, "numericbody2", value = myset[2])
          updateNumericInput(session, "numerictss", value = myset[3])
          updateNumericInput(session, "numerictes", value = myset[4])
          updateNumericInput(session, "numericbinsize", value = myset[5])
          updateNumericInput(session, "numericlabelspaceing", value = myset[6])
        }
      }
      Lines_Labels_List <- LinesLabelsListset(myset[1],
                                              myset[2],
                                              myset[3],
                                              myset[4],
                                              myset[5],
                                              LIST_DATA$x_plot_range[2],
                                              myset[6],
                                              input$numerictssname,
                                              input$numerictesname)
      if (input$selectlineslabels != "" & LIST_DATA$STATE[2] != 0) {
        updateSelectInput(session, "selectlineslabels", selected = "")
        my_pos <-
          suppressWarnings(as.numeric((Lines_Labels_List$mybrakes)))
        my_label <- Lines_Labels_List$mylabels
      }
      # set lable and posistion numbers
      updateTextInput(session,
                      "landlnames",
                      value = paste(Lines_Labels_List$mylabels, collapse = " "))
      updateTextInput(session,
                      "landlposition",
                      value = paste(Lines_Labels_List$mybrakes , collapse = " "))
    }
  )
  
  # checks that number of names == position ----
  observeEvent(c(input$landlnames, input$landlposition), ignoreInit = TRUE, {
    my_pos <-
      suppressWarnings(as.numeric(unlist(
        strsplit(input$landlposition, split = "\\s+")
      )))
    my_label <- unlist(strsplit(input$landlnames, split = "\\s+"))
    if (any(is.na(my_pos))) {
      my_pos <- my_pos[is.na(my_pos)]
      updateTextInput(session, "landlposition", value = my_pos)
    }
    if (length(my_pos) == length(my_label)) {
      shinyjs::enable("actionlineslabels")
      updateActionButton(session, "actionlineslabels", label = "Plot with new Lines and Labels")
    } else {
      updateActionButton(session, "actionlineslabels", label = "Labels must equel # of positions")
      shinyjs::disable("actionlineslabels")
    }
    if (LIST_DATA$STATE[2] == 0) {
      if (length(my_pos) == 0) {
        my_label <- "none"
        my_pos <- LIST_DATA$x_plot_range[2] * 2
      }
      reactive_values$Lines_Labels_List <-
        LinesLabelsListPlot(
          input$numericbody1,
          input$selectbody1color,
          input$selectbody1line,
          input$numericbody2,
          input$selectbody2color,
          input$selectbody2line,
          input$numerictss,
          input$selecttsscolor,
          input$selecttssline,
          input$numerictes,
          input$selecttescolor,
          input$selecttesline,
          my_label,
          my_pos,
          input$selectvlinesize,
          input$selectlinesize,
          input$selectfontsizex,
          input$selectfontsizey,
          input$selectlegendsize,
          input$selectttestlinesize
        )
    }
  })
  
  # action button update lines and labels ----
  observeEvent(input$actionlineslabels, ignoreInit = TRUE, {
    print("action lines and labels")
    if (!is.null(reactive_values$Apply_Math)) {
      reactive_values$Y_Axis_numbers <-
        MyXSetValues(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          c(0,100),
          input$checkboxlog2
        )
      my_step <-
        (max(reactive_values$Y_Axis_numbers) - min(reactive_values$Y_Axis_numbers)) /
        20
      updateNumericInput(session,
                         "numericYRangeHigh",
                         value = round(max(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
      updateNumericInput(session,
                         "numericYRangeLow",
                         value = round(min(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
    }
    my_pos <-
      suppressWarnings(as.numeric(unlist(
        strsplit(input$landlposition, split = "\\s+")
      )))
    my_label <- unlist(strsplit(input$landlnames, split = "\\s+"))
    if (length(my_pos) == 0) {
      my_label <- "none"
      my_pos <- LIST_DATA$x_plot_range[2] * 2
    }
    
    # if tss or tes location make sure there is text
    if(nchar(trimws(input$numerictssname)) == 0 & input$numerictss > 0){
      updateTextInput(session, "numerictssname", value = "TSS")
    }
    if(nchar(trimws(input$numerictesname)) == 0 & input$numerictes > 0){
      updateTextInput(session, "numerictesname", value = "pA")
    }
    
    reactive_values$Lines_Labels_List <-
      LinesLabelsListPlot(
        input$numericbody1,
        input$selectbody1color,
        input$selectbody1line,
        input$numericbody2,
        input$selectbody2color,
        input$selectbody2line,
        input$numerictss,
        input$selecttsscolor,
        input$selecttssline,
        input$numerictes,
        input$selecttescolor,
        input$selecttesline,
        my_label,
        my_pos,
        input$selectvlinesize,
        input$selectlinesize,
        input$selectfontsizex,
        input$selectfontsizey,
        input$selectlegendsize,
        input$selectttestlinesize
      )
  })
  
  # replot with smooth update ----
  observeEvent(input$checkboxsmooth, ignoreInit = TRUE, {
    if (!is.null(reactive_values$Apply_Math) &
        LIST_DATA$STATE[2] != 2) {
      reactive_values$Y_Axis_Lable <-
        YAxisLable(
          input$myMath,
          input$selectplotnrom,
          as.numeric(input$selectplotBinNorm),
          input$checkboxsmooth,
          input$checkboxlog2
        )
      reactive_values$Y_Axis_numbers <-
        MyXSetValues(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          c(0,100),
          input$checkboxlog2
        )
      my_step <-
        (max(reactive_values$Y_Axis_numbers) - min(reactive_values$Y_Axis_numbers)) /
        20
      updateNumericInput(session,
                         "numericYRangeHigh",
                         value = round(max(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
      updateNumericInput(session,
                         "numericYRangeLow",
                         value = round(min(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options,
          reactive_values$Y_Axis_numbers,
          reactive_values$Lines_Labels_List,
          input$checkboxsmooth, reactive_values$Plot_Options_ttest,
          input$checkboxlog2,
          reactive_values$Y_Axis_Lable,
          input$sliderplotOccupancy
        )
    }
  })
  
  # replot with log2 update ----
  observeEvent(input$checkboxlog2, ignoreInit = TRUE, {
    if (!is.null(reactive_values$Apply_Math) &
        LIST_DATA$STATE[2] != 2) {
      reactive_values$Y_Axis_Lable <-
        YAxisLable(
          input$myMath,
          input$selectplotnrom,
          as.numeric(input$selectplotBinNorm),
          input$checkboxsmooth,
          input$checkboxlog2
        )
      reactive_values$Y_Axis_numbers <-
        MyXSetValues(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          c(0,100),
          input$checkboxlog2
        )
      my_step <-
        (max(reactive_values$Y_Axis_numbers) - min(reactive_values$Y_Axis_numbers)) /
        20
      updateNumericInput(session,
                         "numericYRangeHigh",
                         value = round(max(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
      updateNumericInput(session,
                         "numericYRangeLow",
                         value = round(min(reactive_values$Y_Axis_numbers), 4),
                         step = my_step)
      reactive_values$Plot_controler <-
        GGplotLineDot(
          reactive_values$Apply_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Options,
          reactive_values$Y_Axis_numbers,
          reactive_values$Lines_Labels_List,
          input$checkboxsmooth, reactive_values$Plot_Options_ttest,
          input$checkboxlog2,
          reactive_values$Y_Axis_Lable,
          input$sliderplotOccupancy
        )
    }
  })
  
  # quick color set change ----
  observeEvent(input$kbrewer, ignoreInit = TRUE, {
    print("kbrewer")
    if (!is.null(LIST_DATA$gene_info[[1]]) &
        input$kbrewer != "select") {
      kListColorSet <<- brewer.pal(11, input$kbrewer)[-c(4:7)][1:n_distinct(LIST_DATA$gene_info$set)]
      print("kbrewer update")
      LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
        group_by(gene_list) %>% 
        dplyr::mutate(.,mycol=if_else(gene_list == input$selectgenelistoptions, kListColorSet, mycol)) %>% 
        ungroup()
      
      updateColourInput(session, "colourhex", value =
                          paste(dplyr::filter(LIST_DATA$gene_info, 
                                              gene_list == input$selectgenelistoptions & 
                                                set == input$selectdataoption) %>% 
                                  dplyr::select(mycol)))
      if (!is.null(reactive_values$Apply_Math) &
          LIST_DATA$STATE[2] == 1) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options,
            reactive_values$Y_Axis_numbers,
            reactive_values$Lines_Labels_List,
            input$checkboxsmooth, reactive_values$Plot_Options_ttest,
            input$checkboxlog2,
            reactive_values$Y_Axis_Lable,
            input$sliderplotOccupancy
          )
      }
      updateSelectInput(session, "kbrewer", selected = "select")
      reactive_values$Picker_controler <-
        paste("color", dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == input$selectgenelistoptions),
          mycol)$mycol, sep = ":")
    }
  })
  
  
  # new gene list color set quick fix ----
  observeEvent(input$BttnNewColor, ignoreInit = TRUE,{
    print("set same color")
    if (length(names(LIST_DATA$gene_file)) > 1) {
      my_sel <- LIST_DATA$gene_info %>% 
        dplyr::filter(gene_list == names(LIST_DATA$gene_file)[1]) %>% 
        dplyr::select(mycol)
      LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% group_by(gene_list) %>% 
        dplyr::mutate(mycol = my_sel$mycol) %>% ungroup()
      
      updateColourInput(session, "colourhex", value =
                          dplyr::select(dplyr::filter(LIST_DATA$gene_info, 
                                                      gene_list == input$selectgenelistoptions &
                                                        set == input$selectdataoption),mycol)$mycol)
      if (!is.null(reactive_values$Apply_Math) &
          LIST_DATA$STATE[2] != 2) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options,
            reactive_values$Y_Axis_numbers,
            reactive_values$Lines_Labels_List,
            input$checkboxsmooth, reactive_values$Plot_Options_ttest,
            input$checkboxlog2,
            reactive_values$Y_Axis_Lable,
            input$sliderplotOccupancy
          )
      }
    }
  })
  
  # update plot picker ----
  observeEvent(
    reactive_values$Picker_controler,
    ignoreNULL = FALSE,
    ignoreInit = TRUE,
    {
      print("plot pickers update")
      
      pickerlist <- list()
      for (i in names(LIST_DATA$gene_file)) {
        pickerlist[[i]] <-
          list(div(
            style = "margin-bottom: -10px;",
            pickerInput(
              inputId = gsub(" ", "-bensspace2-", gsub("\n", "-bensspace1-", i)),
              label = i,
              width = "99%",
              choices = distinct(LIST_DATA$gene_info, set)$set,
              selected =  dplyr::select(dplyr::filter(LIST_DATA$gene_info, 
                                                      gene_list == i & onoff != 0), 
                                        onoff)$onoff,
              multiple = T,
              options = list(
                `actions-box` = TRUE,
                `selected-text-format` = "count > 0"),
              choicesOpt = list(style = paste("color", 
                                              dplyr::select(dplyr::filter(LIST_DATA$gene_info, 
                                                                          gene_list == i), mycol)$mycol,
                                              sep = ":"))
            )
          ))
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        output$DynamicGenePicker_sort <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Filter")]
        })
        shinyjs::show("showpickersort")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Filter"))){
        shinyjs::hide("showpickersort")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        output$DynamicGenePicker_comparisons <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Gene_List_")]
        })
        shinyjs::show("showpickercomparisons")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Gene_List_"))){
        shinyjs::hide("showpickercomparisons")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        output$DynamicGenePicker_ratio <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Ratio_")]
        })
        shinyjs::show("showpickerratio")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Ratio_"))){
        shinyjs::hide("showpickerratio")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_"))){
        output$DynamicGenePicker_clusters <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_")]
        })
        shinyjs::show("showpickercluster")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^Group_|^Cluster_"))){
        shinyjs::hide("showpickercluster")
      }
      if(any(str_detect(names(LIST_DATA$gene_file),"^CDF"))){
        output$DynamicGenePicker_cdf <- renderUI({
          pickerlist[str_detect(names(LIST_DATA$gene_file),"^CDF")]
        })
        shinyjs::show("showpickercdf")
      } else if(!any(str_detect(names(LIST_DATA$gene_file),"^CDF"))){
        shinyjs::hide("showpickercdf")
      }
      output$DynamicGenePicker_main <- renderUI({
        pickerlist[!str_detect(names(LIST_DATA$gene_file),"^Filter|^Gene_List_|^Ratio_|^Group_|^Cluster_|^CDF")]
      })
      
    }
  )
  
  # Remove data file ----
  observeEvent(input$actionremovefile, ignoreInit = TRUE, {
    print("remove file")
    LIST_DATA <<-
      RemoveFile(LIST_DATA,
                 input$selectdataoption,
                 input$checkboxremovefile)
    if (LIST_DATA$STATE[1] != 0) {
      ff <- distinct(LIST_DATA$gene_info, set)$set
      updateSelectInput(session,
                        "selectdataoption",
                        choices = ff)
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_file),
        selected = names(LIST_DATA$gene_file)[1]
      )
    } else {
      updateSelectInput(session,
                        "selectdataoption",
                        choices = '')
      updateSelectInput(session,
                        "selectgenelistoptions",
                        choices = 'Load Data File')
      shinyjs::hide("filegene1")
      shinyjs::hide("downloadGeneList")
      shinyjs::hide("checkboxsavesplit")
      shinyjs::hide("filecolor")
      shinyjs::hide("hiddensave")
      shinyjs::hide("startoff")
      shinyjs::hide('actiongenelistsdatatable')
      shinyjs::hide('genelists1table')
      shinyjs::hide('genelists2table')
      shinyjs::hide('genelists3table')
      shinyjs::hide('actionsortdatatable')
      shinyjs::hide('sorttable')
      shinyjs::hide('ratio1table')
      shinyjs::hide('ratio2table')
      shinyjs::hide('ratio3table')
      shinyjs::hide('plotcluster')
      shinyjs::hide("cluster1table")
      shinyjs::hide("cluster2table")
      shinyjs::hide("cluster3table")
      shinyjs::hide("cluster4table")
      shinyjs::hide('cdftable')
      shinyjs::hide('plotcdf')
      reactive_values$Apply_Math <- NULL
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
      updateAwesomeCheckbox(session, "checkboxremovefile", value = FALSE)
    }
    reactive_values$binset <- FALSE
    reactive_values$Picker_controler <- 
      c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
  })
  
  # Remove gene list ----
  observeEvent(input$actionremovegene, ignoreInit = TRUE, {
    print("remove gene list")
    shinyjs::hide('actiongenelistsdatatable')
    shinyjs::hide('genelists1table')
    shinyjs::hide('genelists2table')
    shinyjs::hide('genelists3table')
    shinyjs::hide('actionsortdatatable')
    shinyjs::hide('sorttable')
    shinyjs::hide('ratio1table')
    shinyjs::hide('ratio2table')
    shinyjs::hide('ratio3table')
    shinyjs::hide('plotcluster')
    shinyjs::hide("cluster1table")
    shinyjs::hide("cluster2table")
    shinyjs::hide("cluster3table")
    shinyjs::hide("cluster4table")
    shinyjs::hide('cdftable')
    shinyjs::hide('plotcdf')
    LIST_DATA <<-
      RemoveGeneList(LIST_DATA, input$selectgenelistoptions)
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_file),
      selected = names(LIST_DATA$gene_file)[1]
    )
  })
  
  observeEvent(c(input$pickernumerator, input$adddata,
                 input$pickerdenominator), {
                   if (input$pickernumerator != "") {
                     updateTextInput(session, "textnromname",value = paste(input$pickernumerator, input$adddata,
                                                                           input$pickerdenominator))
                     output$valueboxnormfile <- renderValueBox({
                       valueBox("0%",
                                "Done",
                                icon = icon("cogs"),
                                color = "yellow")
                     })
                   }
                 })
  
  # create norm file ----
  observeEvent(input$actionnorm, ignoreInit = TRUE, {
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- MakeNormFile(
                     LIST_DATA,
                     input$pickernumerator,
                     input$pickerdenominator,
                     input$radiogenebygene,
                     input$checkboxnormzero,
                     input$adddata,
                     input$textnromname
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      updatePickerInput(session,
                        "pickernumerator", selected = "",
                        choices = distinct(LIST_DATA$gene_info, set)$set,
                        choicesOpt = list(style = paste("color", 
                                                        dplyr::select(
                                                          dplyr::filter(LIST_DATA$gene_info,
                                                                        gene_list == names(
                                                                          LIST_DATA$gene_file)[1]), 
                                                          mycol)$mycol, 
                                                        sep = ":")))
      updatePickerInput(session,
                        "pickerdenominator", selected = "",
                        choices = distinct(LIST_DATA$gene_info, set)$set,
                        choicesOpt = list(style = paste("color", 
                                                        dplyr::select(
                                                          dplyr::filter(LIST_DATA$gene_info,
                                                                        gene_list == names(
                                                                          LIST_DATA$gene_file)[1]),
                                                          mycol)$mycol, 
                                                        sep = ":")))
      updateTextInput(session, "textnromname", value = "")
      output$valueboxnormfile <- renderValueBox({
        valueBox(
          "Done",
          paste("Compleat n =", n_distinct(LIST_DATA$gene_file[[1]]$use$gene)),
          icon = icon("thumbs-up", lib = "glyphicon"),
          color = "green"
        )
      })
      ff <- distinct(LIST_DATA$gene_info, set)$set
      updateSelectInput(session,
                        "selectdataoption",
                        choices = ff)
      # update main plot picker with new file name trigger
      reactive_values$Picker_controler <- 
        c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
    } else {
      #no new data file created
      output$valueboxnormfile <- renderValueBox({
        valueBox("0%",
                 "Done",
                 icon = icon("cogs"),
                 color = "red")
      })
    }
  })
  
  # Gene action ----
  observeEvent(input$actiongenelists, {
    print("gene lists action")
    shinyjs::hide('actiongenelistsdatatable')
    shinyjs::hide('genelists1table')
    shinyjs::hide('genelists2table')
    shinyjs::hide('genelists3table')
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- IntersectGeneLists(LIST_DATA,
                                            input$pickergenelists)
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
      ol <- input$pickergenelists
      if (!any(ol %in% names(LIST_DATA$gene_file))) {
        ol <- grep("Gene_List_", names(LIST_DATA$gene_file), value = TRUE)
      } else {
        
      }
      updateSelectInput(
        session,
        "selectsortfile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      shinyjs::show('actiongenelistsdatatable')
      if (any(grep("Gene_List_intersect\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxgene1 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Gene_List_intersect\nn =",
                                                 names(LIST_DATA$gene_file))]]$use$gene),
            "Gene List intersect",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxgene1 <- renderValueBox({
          valueBox(0,
                   "Gene List intersect",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (any(grep("Gene_List_Total\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxgene2 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Gene_List_Total\nn =", names(LIST_DATA$gene_file))]]$full$gene),
            "Gene List Total",
            icon = icon("list"),
            color = "yellow"
          )
        })
      } else{
        output$valueboxgene2 <- renderValueBox({
          valueBox(0,
                   "Gene List Total",
                   icon = icon("list"),
                   color = "yellow")
        })
      }
      if (any(grep("Gene_List_exclusive\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxgene3 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Gene_List_exclusive\nn =",
                                                 names(LIST_DATA$gene_file))]]$use$gene),
            "Gene List exclusive",
            icon = icon("list"),
            color = "red"
          )
        })
      } else{
        output$valueboxgene3 <- renderValueBox({
          valueBox(0,
                   "Gene List exclusive",
                   icon = icon("list"),
                   color = "red")
        })
      }
    } else {
      output$valueboxgene1 <- renderValueBox({
        valueBox(0,
                 "Gene List intersect",
                 icon = icon("list"),
                 color = "green")
      })
      output$valueboxgene2 <- renderValueBox({
        valueBox(0,
                 "Gene List Total",
                 icon = icon("list"),
                 color = "yellow")
      })
      output$valueboxgene3 <- renderValueBox({
        valueBox(0,
                 "Gene List exclusive",
                 icon = icon("list"),
                 color = "red")
      })
      return()
    }
  })
  
  # Gene lists show gene list ----
  observeEvent(input$actiongenelistsdatatable, ignoreInit = TRUE, {
    print("generiate gene lists table")
    shinyjs::hide('actiongenelistsdatatable')
    if (any(grep("Gene_List_intersect\nn =", names(LIST_DATA$gene_file)) >
            0)) {
      newnames1 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_intersect\nn =",
               names(LIST_DATA$gene_file),
               value = TRUE
             ))
      mytab <- "Intersected Gene Lists"
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists1table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_intersect\nn =",
                                                     names(LIST_DATA$gene_file))]]$use,
                           rownames = FALSE,
                           colnames = newnames1,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Gene_List_intersect\nn =",
                                                               names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show('genelists1table')
    } else {
      output$genelists1table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Gene_List_exclusive n = 0",
            options = list(searching = FALSE)
          )
        )
      mytab <- "Total Gene Lists"
    }
    if (any(grep("Gene_List_Total\nn =", names(LIST_DATA$gene_file)) >
            0)) {
      newnames2 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_Total\nn =",
               names(LIST_DATA$gene_file),
               value = TRUE
             ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists2table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_Total\nn =",
                                                     names(LIST_DATA$gene_file))]]$use,
                           rownames = FALSE,
                           colnames = newnames2,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Gene_List_Total\nn =",
                                                               names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show('genelists2table')
    } else {
      output$genelists2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Gene_List_exclusive n = 0",
            options = list(searching = FALSE)
          )
        )
      if (mytab == "Total Gene Lists") {
        mytab <- "Exclusive Gene Lists"
      }
    }
    if (any(grep("Gene_List_exclusive\nn =", names(LIST_DATA$gene_file)) >
            0)) {
      newnames3 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_exclusive\nn =",
               names(LIST_DATA$gene_file),
               value = T
             ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists3table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_exclusive\nn =",
                                                     names(LIST_DATA$gene_file))]]$use,
                           rownames = FALSE,
                           colnames = newnames3,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Gene_List_exclusive\nn =",
                                                               names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show('genelists3table')
    } else {
      output$genelists3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Gene_List_exclusive n = 0",
            options = list(searching = FALSE)
          )
        )
      if (mytab == "Exclusive Gene Lists") {
        mytab <- "Total Gene Lists"
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
      choices = distinct(LIST_DATA$gene_info, set)$set,
      selected = reactive_values$pickerfile_controler,
      choicesOpt = list(style = paste("color", 
                                      dplyr::select(
                                        dplyr::filter(LIST_DATA$gene_info,
                                                      gene_list == input$selectsortfile), 
                                        mycol)$mycol, 
                                      sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  # sort sum tool action ----
  observeEvent(input$actionsorttool, {
    print("sort tool")
    shinyjs::hide('actionsortdatatable')
    shinyjs::hide('sorttable')
    
    if (input$slidersortpercent < 50 &
        input$selectsorttop == "Middle%") {
      updateSliderInput(session, "slidersortpercent", value = 50)
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- FilterTop(
                     LIST_DATA,
                     input$selectsortfile,
                     input$pickersortfile,
                     input$slidersortbinrange[1],
                     input$slidersortbinrange[2],
                     input$slidersortpercent,
                     input$selectsorttop,
                     input$checkboxfilterall
                   )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
      LD <- LIST_DATA
      mylist <- last(grep("^Filter", names(LIST_DATA$gene_file)))
      LD$gene_info <- LD$gene_info %>%
        dplyr::mutate(onoff=if_else(gene_list == names(LD$gene_file)[mylist] &
                                      set %in% input$pickersortfile, set, "0"))
      list_data_frame <- Active_list_data(LD)
      if (!is_empty(list_data_frame)) {
        Apply_Cluster_Math <-
          ApplyMath(
            list_data_frame,
            input$myMathcluster,
            input$radioplotnromcluster,
            as.numeric(input$selectplotBinNormcluster)
          )
        reactive_values$Plot_controler_sort_min <- ggplot()
        reactive_values$Plot_controler_sort_max <- ggplot()
        gp1 <-
          ggplot(Apply_Cluster_Math ,aes(as.numeric(bin),value,color=set)) +
          geom_line() +
          ylab("Mean bin value") +
          theme(legend.position="bottom",
                legend.title = element_blank(),
                axis.title.x=element_blank())
        print(gp1)
        reactive_values$Plot_controler_sort_min <- gp1
      }
      shinyjs::show('actionsortdatatable')
      if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$use$gene),
            "Gene List Filter",
            icon = icon("list"),
            color = "green"
          )
        })
      } else {
        output$valueboxsort <- renderValueBox({
          valueBox(0,
                   "Gene List Filter",
                   icon = icon("list"),
                   color = "green")
        })
      }
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
    }
    # updating select and keeping track if sort on sort
    ol <- input$selectsortfile
    if (!ol %in% names(LIST_DATA$gene_file)) {
      ol <- grep("^Filter", names(LIST_DATA$gene_file), value = TRUE)
      reactive_values$pickerfile_controler <- input$pickersortfile
    } else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectsortfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
  })
  # sort % numiric controler ----
  observeEvent(c(input$numericsortmin,input$numericsortmax), ignoreInit = TRUE,{
    if (!is.numeric(input$numericsortmin)) {
      updateNumericInput(session, "numericsortmin", value = 1)
    }
    if (!is.numeric(input$numericsortmax)) {
      updateNumericInput(session, "numericsortmax", value = 99.5)
    }
    if (input$numericsortmin < 0 | input$numericsortmin > input$numericsortmax) {
      updateNumericInput(session, "numericsortmin", value = 1)
    }
    if (input$numericsortmax < input$numericsortmin | input$numericsortmax > 100) {
      updateNumericInput(session, "numericsortmax", value = 99.5)
    }
  })
  
  # sort min max between % tool action ----
  observeEvent(input$actionsortper, ignoreInit = TRUE, {
    print("sort % tool")
    
    sortmin <- FilterPer(LIST_DATA, 
                         input$selectsortfile,
                         input$pickersortfile,
                         input$slidersortbinrange,
                         input$checkboxfilterall,
                         c(input$numericsortmin,input$numericsortmax),
                         input$selectsortper,
                         input$checkboxfilterall)
    
    if(!is_empty(sortmin$sortplot)){
      LIST_DATA <<- sortmin
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      } 
      if(input$selectsortper != "max%" ){  
        gp1 <-
          ggplot(sortmin$sortplot %>% dplyr::filter(!is.na(my_p_1)) ,aes(as.numeric(bin),my_p_1,color=set)) + 
          geom_line() + 
          ylab("Min % value") +
          theme(legend.position="bottom", 
                legend.title = element_blank(),
                axis.title.x=element_blank())
      } else{
        gp1 <- ggplot()
      }
      if (input$selectsortper != "min%" ){
        gp2 <-
          ggplot(sortmin$sortplot %>% dplyr::filter(!is.na(my_p_2)),aes(as.numeric(bin),my_p_2,color=set)) + 
          geom_line() + 
          ylab("Max % value") +
          theme(legend.position="bottom", 
                legend.title = element_blank(),
                axis.title.x=element_blank())
      } else {
        gp2 <- ggplot()
      }
      if (input$selectsortper == "min%" ){
        print(gp1)
      } else if(input$selectsortper == "max%" ){
        print(gp2)
      } else{
        print(gp1 + gp2 + plot_layout(ncol = 1))
      }
      reactive_values$Plot_controler_sort_min <- gp1
      reactive_values$Plot_controler_sort_max <- gp2
      if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxsort <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$use$gene),
            "Gene List Filter",
            icon = icon("list"),
            color = "green"
          )
        })
      } else {
        output$valueboxsort <- renderValueBox({
          valueBox(0,
                   "Gene List Filter",
                   icon = icon("list"),
                   color = "green")
        })
      }
    } else {
      output$valueboxsort <- renderValueBox({
        valueBox(0,
                 "Gene List Filter",
                 icon = icon("list"),
                 color = "green")
      })
      reactive_values$Plot_controler_sort_min <-
        ggplot()
      reactive_values$Plot_controler_sort_max <-
        ggplot()
    }
    ol <- input$selectsortfile
    if (!ol %in% names(LIST_DATA$gene_file)) {
      ol <- grep("^Filter", names(LIST_DATA$gene_file), value = TRUE)
      reactive_values$pickerfile_controler <- input$pickersortfile
    } else {
      reactive_values$pickerfile_controler <- ""
    }
    updateSelectInput(
      session,
      "selectsortfile",
      choices = names(LIST_DATA$gene_file),
      selected = ol
    )
  })
  
  # Filter gene list show data table ----
  observeEvent(input$actionsortdatatable, ignoreInit = TRUE, {
    print("show data table")
    if (any(grep("^Filter", names(LIST_DATA$gene_file)) > 0)) {
      newnames <-
        gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full))
      dt <- datatable(
        rowid_to_column(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full),
        rownames = FALSE,
        colnames = c("ID", strtrim(newnames, 30)),
        class = 'cell-border stripe compact',
        filter = 'top',
        caption = LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$info,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = FALSE,
          width = 5,
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
        
      ) %>% formatPercentage(names(LIST_DATA$gene_file[[last(grep("^Filter", names(LIST_DATA$gene_file)))]]$full)[-1])
    } else {
      dt <- datatable(
        LIST_DATA$gene_file[[1]]$empty,
        rownames = FALSE,
        colnames = "none",
        options = list(searching = FALSE)
      )
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   output$sorttable <- DT::renderDataTable(dt)
                 })
    shinyjs::hide('actionsortdatatable')
    shinyjs::show('sorttable')
  })
  
  # sort tool gene list $use ----
  observeEvent(input$sorttable_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 oldname <-
                   last(grep("^Filter ", names(LIST_DATA$gene_file)))
                 oldname1 <- names(LIST_DATA$gene_file)[oldname] %>% str_split_fixed(., " n = ",n=2)
                 newname <-
                   paste(oldname1[1], "n =",
                         length(input$sorttable_rows_all))
                 if (newname != names(LIST_DATA$gene_file)[oldname]) {
                   print("sort filter $use")
                   names(LIST_DATA$gene_file)[oldname] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$sorttable_rows_all])
                   ol <- input$selectsortfile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                     reactive_values$pickerfile_controler <-
                       input$pickersortfile
                   } else {
                     reactive_values$pickerfile_controler <- ""
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   }
                   updateSelectInput(
                     session,
                     "selectsortfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxsort <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Gene List Filter",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  # ratio tool picker control ----
  observeEvent(input$selectratiofile, ignoreInit = TRUE, {
    print("ratio picker update")
    if (reactive_values$pickerfile_controler[1] == "") {
      reactive_values$pickerfile_controler <- c("", "")
    }
    updatePickerInput(
      session,
      "pickerratio1file",
      choices = distinct(LIST_DATA$gene_info, set)$set,
      selected = reactive_values$pickerfile_controler[1],
      choicesOpt = list(style = paste("color", dplyr::select(
        dplyr::filter(LIST_DATA$gene_info,
                      gene_list == input$selectratiofile),
        mycol)$mycol, sep = ":"))
    )
    updatePickerInput(
      session,
      "pickerratio2file",
      choices = c("None", distinct(LIST_DATA$gene_info, set)$set),
      selected = reactive_values$pickerfile_controler[2],
      choicesOpt = list(style = paste("color", c(
        "#000000",
        dplyr::select(
          dplyr::filter(LIST_DATA$gene_info,
                        gene_list == input$selectratiofile),
          mycol)$mycol), sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  # ratio tool gene lists $use ----
  observeEvent(input$ratio1table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <-
                   paste("Ratio_Up_file1\nn =",
                         length(input$ratio1table_rows_all))
                 oldname <-
                   grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldname]) {
                   print("ratio1 filter $use")
                   names(LIST_DATA$gene_file)[oldname] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$ratio1table_rows_all])
                   
                   ol <- input$selectratiofile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       c(input$pickerratio1file, input$pickerratio2file)
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   } else {
                     reactive_values$pickerfile_controler <- ""
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   }
                   updateSelectInput(
                     session,
                     "selectratiofile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxratio1 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Ratio Up file1",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  observeEvent(input$ratio2table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <-
                   paste("Ratio_Down_file1\nn =",
                         length(input$ratio2table_rows_all))
                 oldname <-
                   grep("Ratio_Down_file1\nn =", names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldname]) {
                   print("ratio2 filter $use")
                   names(LIST_DATA$gene_file)[oldname] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$ratio2table_rows_all])
                   
                   ol <- input$selectratiofile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       c(input$pickerratio1file, input$pickerratio2file)
                   } else {
                     reactive_values$pickerfile_controler <- ""
                   }
                   updateSelectInput(
                     session,
                     "selectratiofile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxratio2 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use$gene),
                       "Ratio Up file2",
                       icon = icon("list"),
                       color = "blue"
                     )
                   })
                 }
               })
  
  observeEvent(input$ratio3table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <-
                   paste("Ratio_No_Diff\nn =", length(input$ratio3table_rows_all))
                 oldname <-
                   grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldname]) {
                   print("no ratio filter $use")
                   names(LIST_DATA$gene_file)[oldname] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$ratio3table_rows_all])
                   
                   ol <- input$selectratiofile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       c(input$pickerratio1file, input$pickerratio2file)
                   } else {
                     reactive_values$pickerfile_controler <- ""
                   }
                   updateSelectInput(
                     session,
                     "selectratiofile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxratio3 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use$gene),
                       "Ratio Up file3",
                       icon = icon("list"),
                       color = "yellow"
                     )
                   })
                 }
               })
  
  # Ratio tool action ----
  observeEvent(input$actionratiotool, ignoreInit = TRUE, {
    print("ratio tool action")
    shinyjs::hide('ratio1table')
    shinyjs::hide('ratio2table')
    shinyjs::hide('ratio3table')
    if (is.numeric(input$numericratio)) {
      if (input$numericratio < 0) {
        updateNumericInput(session, "numericratio", value = 2)
      }
    } else {
      updateNumericInput(session, "numericratio", value = 2)
    }
    if (input$sliderbinratio2[1] == 0 & input$sliderbinratio2[2] > 0) {
      # showModal(modalDialog(
      #   title = "Information message",
      #   paste("Bins regions should not overlab, \nBins set to 1/4 3/4"),
      #   size = "s",
      #   easyClose = TRUE
      # ))
      updateSliderInput(session,
                        "sliderbinratio2",
                        value = c(0,0))
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
                       input$checkratiozero,
                       input$sliderRatioBinNorm
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
      shinyjs::show('actionratiodatatable')
      
      ol <- input$selectratiofile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <-
          grep(strsplit(ol, "\nn = ")[[1]][1],
               names(LIST_DATA$gene_file),
               value = TRUE)
        reactive_values$pickerfile_controler <-
          c(input$pickerratio1file, input$pickerratio2file)
      } else {
        reactive_values$pickerfile_controler <- ""
      }
      updateSelectInput(
        session,
        "selectratiofile",
        choices = names(LIST_DATA$gene_file),
        selected = ol
      )
      if (any(grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio1 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file))]]$use$gene),
            "Ratio Up file1",
            icon = icon("list"),
            color = "green"
          )
        })
      } else{
        output$valueboxratio1 <- renderValueBox({
          valueBox(0,
                   "Ratio Up file1",
                   icon = icon("list"),
                   color = "green")
        })
      }
      if (any(grep("Ratio_Down_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio2 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_Down_file1\nn =",
                                                 names(LIST_DATA$gene_file))]]$use$gene),
            "Ratio Down file1",
            icon = icon("list"),
            color = "blue"
          )
        })
      } else{
        output$valueboxratio2 <- renderValueBox({
          valueBox(0,
                   "Ratio Down file1",
                   icon = icon("list"),
                   color = "blue")
        })
      }
      if (any(grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file)) > 0)) {
        output$valueboxratio3 <- renderValueBox({
          valueBox(
            n_distinct(LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file))]]$use$gene),
            "Ratio No Diff",
            icon = icon("list"),
            color = "yellow"
          )
        })
      } else{
        output$valueboxratio3 <- renderValueBox({
          valueBox(0,
                   "Ratio No Diff",
                   icon = icon("list"),
                   color = "yellow")
        })
      }
      if(!is.null(LIST_DATA$boxRatio)){
        my_range <- range(LIST_DATA$boxRatio$Ratio,na.rm = T) 
        updateNumericInput(session, "textboxmaxratio",
                           value = my_range[2])
        updateNumericInput(session, "textboxminratio",
                           value = my_range[1])
      } else {
        updateNumericInput(session, "textboxmaxratio",
                           value = 0)
        updateNumericInput(session, "textboxminratio",
                           value = 0)
      }
    } else {
      output$valueboxratio1 <- renderValueBox({
        valueBox(0,
                 "Ratio Up file1",
                 icon = icon("list"),
                 color = "green")
      })
      output$valueboxratio2 <- renderValueBox({
        valueBox(0,
                 "Ratio Down file1",
                 icon = icon("list"),
                 color = "blue")
      })
      output$valueboxratio3 <- renderValueBox({
        valueBox(0,
                 "Ratio No Diff",
                 icon = icon("list"),
                 color = "yellow")
      })
      return()
    }
  })
  
  # Ratio show gene list ----
  observeEvent(input$actionratiodatatable, ignoreInit = TRUE, {
    print("generiate ratio table")
    shinyjs::hide('actionratiodatatable')
    if (any(grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
      newnames1 <-
        gsub("\n", " ", grep(
          "Ratio_Up_file1\nn =",
          names(LIST_DATA$gene_file),
          value = TRUE
        ))
      mytab <- "Up Fold Change file 1"
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$ratio1table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn =",
                                                     names(LIST_DATA$gene_file))]]$full,
                           rownames = FALSE,
                           colnames = newnames1,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn", names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show('ratio1table')
    } else {
      output$ratio1table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = strtrim(newnames1, 24),
            options = list(searching = FALSE)
          )
        )
      mytab <- "Up Fold Change file 2"
    }
    if (any(grep("Ratio_Down_file1\nn =", names(LIST_DATA$gene_file)) > 0)) {
      newnames2 <-
        gsub("\n",
             " ",
             grep(
               "Ratio_Down_file1\nn =",
               names(LIST_DATA$gene_file),
               value = TRUE
             ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$ratio2table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Ratio_Down_file1\nn =",
                                                     names(LIST_DATA$gene_file))]]$full,
                           rownames = FALSE,
                           colnames = newnames2,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Ratio_Down_file1\nn", names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show('ratio2table')
    } else {
      output$ratio2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = strtrim(newnames2, 24),
            options = list(searching = FALSE)
          )
        )
      if (mytab == "Up Fold Change file 2") {
        mytab <- "No Fold Change"
      }
    }
    if (any(grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_file)) > 0)) {
      newnames3 <-
        gsub("\n", " ", grep(
          "Ratio_No_Diff\nn =",
          names(LIST_DATA$gene_file),
          value = T
        ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$ratio3table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn =",
                                                     names(LIST_DATA$gene_file))]]$full,
                           rownames = FALSE,
                           colnames = newnames3,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn", names(LIST_DATA$gene_file))]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show('ratio3table')
    } else {
      output$ratio3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = "Ratio_No_Diff n = 0",
            options = list(searching = FALSE)
          )
        )
      if (mytab == "No Fold Change") {
        mytab <- "Up Fold Change file 1"
      }
    }
    updateTabItems(session, "ratiotooltab", mytab)
  })
  
  # update Ratio boxplot ----
  observeEvent(c(input$textboxmaxratio, input$textboxminratio,
                 input$checkboxoutlierratio),ignoreInit = TRUE, ignoreNULL = TRUE,{
                   if(!is.null(LIST_DATA$boxRatio)){
                     my_range <- c(ceiling(input$textboxmaxratio), floor(input$textboxminratio)) 
                     if(input$checkboxoutlierratio){
                       gb <- ggboxplot(LIST_DATA$boxRatio, x= "set", y = "Ratio", 
                                       color="set",short.panel.labs = FALSE, notch = T,
                                       outlier.shape = NA)
                     } else{
                       gb <- ggboxplot(LIST_DATA$boxRatio, x= "set", y = "Ratio", 
                                       color="set",short.panel.labs = FALSE, notch = T)
                     }
                     if(n_distinct(LIST_DATA$boxRatio$set)>1){
                       # add remove outlier and set range
                       combn(unique(LIST_DATA$boxRatio$set),2) -> tt
                       my_comparisons2 <- list()
                       for(i in 1:ncol(tt)){
                         my_comparisons2[[i]] <- (c(tt[1,i],tt[2,i]))
                       }
                       gb <- gb +
                         stat_compare_means(comparisons = my_comparisons2, method = "t.test",
                                            label.y = my_range[1])
                       
                       
                     } 
                     print(gb + coord_cartesian(ylim = my_range))
                     reactive_values$Plot_controler_ratio <- gb +
                       coord_cartesian(ylim = my_range)
                   } 
                   
                 })
  
  # cluster tool picker control ----
  observeEvent(input$selectclusterfile, ignoreInit = TRUE, {
    print("cluster picker update")
    shinyjs::hide('plotcluster')
    shinyjs::hide("cluster1table")
    shinyjs::hide("cluster2table")
    shinyjs::hide("cluster3table")
    shinyjs::hide("cluster4table")
    updatePickerInput(
      session,
      "pickerclusterfile",
      choices = distinct(LIST_DATA$gene_info, set)$set,
      selected = reactive_values$pickerfile_controler,
      choicesOpt = list(style = paste("color", dplyr::select(
        dplyr::filter(LIST_DATA$gene_info,
                      gene_list == input$selectclusterfile),
        mycol)$mycol, sep = ":"))
    )
    reactive_values$pickerfile_controler <- ""
  })
  
  # cluster tool gene lists $use ----
  observeEvent(input$cluster1table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 print("cluster tool gene lists $use 1")
                 oldnamenum <-
                   grep(paste0(reactive_values$clustergroups, "1\nn ="),
                        names(LIST_DATA$gene_file))
                 newname <-
                   paste0(reactive_values$clustergroups,
                          "1\nn = ",
                          length(input$cluster1table_rows_all))
                 print(newname)
                 if (newname != names(LIST_DATA$gene_file)[oldnamenum]) {
                   oldname <- names(LIST_DATA$gene_file)[oldnamenum]
                   print("cluster1 filter $use")
                   names(LIST_DATA$gene_file)[oldnamenum] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$cluster1table_rows_all])
                   LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
                     dplyr::mutate(gene_list=if_else(gene_list == oldname,newname,gene_list))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  LD <- LIST_DATA
                                  LD$gene_info <- LD$gene_info %>%
                                    dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Group_|^Cluster_") &
                                                                  set == input$pickerclusterfile, set, "0"))
                                  list_data_frame <- Active_list_data(LD)
                                  if (!is_empty(list_data_frame)) {
                                    ListColorSet <- brewer.pal(4,"Dark2")[1:n_distinct(distinct(list_data_frame,set,gene_list,plot_set))]
                                    LD$gene_info <- list_data_frame %>% 
                                      distinct(set,gene_list,plot_set) %>% 
                                      dplyr::mutate(mycol=ListColorSet)%>%
                                      full_join(LD$gene_info,.,by=c("set","gene_list")) %>% 
                                      dplyr::filter(!is.na(set)) %>% 
                                      dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                      dplyr::mutate(mycol=if_else(is.na(plot_set.y),mycol.x,mycol.y)) %>% 
                                      dplyr::select(-plot_set.y,-plot_set.x,-mycol.x,-mycol.y) %>%
                                      distinct()
                                    Apply_Cluster_Math <-
                                      ApplyMath(
                                        list_data_frame,
                                        input$myMathcluster,
                                        input$radioplotnromcluster,
                                        as.numeric(input$selectplotBinNormcluster)
                                      )
                                  }
                                })
                   if (!is_empty(list_data_frame)) {
                     Plot_Cluster_Options <-
                       MakePlotOptionFrame(LD)
                     Y_Axis_Cluster_numbers <-
                       MyXSetValues(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         log_2 = input$checkboxlog2cluster
                       )
                     reactive_values$Plot_controler_cluster <-
                       GGplotLineDot(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         Plot_Cluster_Options,
                         Y_Axis_Cluster_numbers,
                         reactive_values$Lines_Labels_List,
                         input$checkboxsmoothcluster,reactive_values$Plot_Options_ttest,
                         input$checkboxlog2cluster,
                         isolate(
                           YAxisLable(
                             input$myMathcluster,
                             input$radioplotnromcluster,
                             as.numeric(input$selectplotBinNorm),
                             input$checkboxsmoothcluster,
                             input$checkboxlog2cluster
                           )
                         ),
                         input$sliderplotOccupancy
                       )
                   }
                   
                   ol <- input$selectclusterfile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       input$pickerclusterfile
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   } else {
                     reactive_values$pickerfile_controler <- ""
                     reactive_values$Picker_controler <- 
                       c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   }
                   updateSelectInput(
                     session,
                     "selectclusterfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxcluster1 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Cluster 1",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  observeEvent(input$cluster2table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <-
                   paste0(reactive_values$clustergroups,
                          "2\nn = ",
                          length(input$cluster2table_rows_all))
                 oldnamenum <-
                   grep(paste0(reactive_values$clustergroups, "2\nn ="),
                        names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldnamenum]) {
                   print("cluster2 filter $use")
                   oldname <- names(LIST_DATA$gene_file)[oldnamenum]
                   names(LIST_DATA$gene_file)[oldnamenum] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$cluster2table_rows_all])
                   LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
                     dplyr::mutate(gene_list=if_else(gene_list == oldname,newname,gene_list))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  LD <- LIST_DATA
                                  LD$gene_info <- LD$gene_info %>%
                                    dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Group_|^Cluster_") &
                                                                  set == input$pickerclusterfile, set, "0"))
                                  list_data_frame <- Active_list_data(LD)
                                  if (!is_empty(list_data_frame)) {
                                    ListColorSet <- brewer.pal(4,"Dark2")[1:n_distinct(distinct(list_data_frame,set,gene_list,plot_set))]
                                    LD$gene_info <- list_data_frame %>% 
                                      distinct(set,gene_list,plot_set) %>% 
                                      dplyr::mutate(mycol=ListColorSet)%>%
                                      full_join(LD$gene_info,.,by=c("set","gene_list")) %>% 
                                      dplyr::filter(!is.na(set)) %>% 
                                      dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                      dplyr::mutate(mycol=if_else(is.na(plot_set.y),mycol.x,mycol.y)) %>% 
                                      dplyr::select(-plot_set.y,-plot_set.x,-mycol.x,-mycol.y) %>%
                                      distinct() 
                                    Apply_Cluster_Math <-
                                      ApplyMath(
                                        list_data_frame,
                                        input$myMathcluster,
                                        input$radioplotnromcluster,
                                        as.numeric(input$selectplotBinNormcluster)
                                      )
                                  }
                                })
                   if (!is_empty(list_data_frame)) {
                     Plot_Cluster_Options <-
                       MakePlotOptionFrame(LD)
                     Y_Axis_Cluster_numbers <-
                       MyXSetValues(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         log_2 = input$checkboxlog2cluster
                       )
                     reactive_values$Plot_controler_cluster <-
                       GGplotLineDot(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         Plot_Cluster_Options,
                         Y_Axis_Cluster_numbers,
                         reactive_values$Lines_Labels_List,
                         input$checkboxsmoothcluster,reactive_values$Plot_Options_ttest,
                         input$checkboxlog2cluster,
                         isolate(
                           YAxisLable(
                             input$myMathcluster,
                             input$radioplotnromcluster,
                             as.numeric(input$selectplotBinNorm),
                             input$checkboxsmoothcluster,
                             input$checkboxlog2cluster
                           )
                         ),
                         input$sliderplotOccupancy
                       )
                   }
                   
                   ol <- input$selectclusterfile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       input$pickerclusterfile
                   } else {
                     reactive_values$pickerfile_controler <- ""
                   }
                   updateSelectInput(
                     session,
                     "selectclusterfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxcluster2 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Cluster 2",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  observeEvent(input$cluster3table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <-
                   paste0(reactive_values$clustergroups,
                          "3\nn = ",
                          length(input$cluster3table_rows_all))
                 oldnamenum <-
                   grep(paste0(reactive_values$clustergroups, "3\nn ="),
                        names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldnamenum]) {
                   print("cluster3 filter $use")
                   oldname <- names(LIST_DATA$gene_file)[oldnamenum]
                   names(LIST_DATA$gene_file)[oldnamenum] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$cluster3table_rows_all])
                   LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
                     dplyr::mutate(gene_list=if_else(gene_list == oldname,newname,gene_list))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  LD <- LIST_DATA
                                  LD$gene_info <- LD$gene_info %>%
                                    dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Group_|^Cluster_") &
                                                                  set == input$pickerclusterfile, set, "0"))
                                  list_data_frame <- Active_list_data(LD)
                                  if (!is_empty(list_data_frame)) {
                                    ListColorSet <- brewer.pal(4,"Dark2")[1:n_distinct(distinct(list_data_frame,set,gene_list,plot_set))]
                                    LD$gene_info <- list_data_frame %>% 
                                      distinct(set,gene_list,plot_set) %>% 
                                      dplyr::mutate(mycol=ListColorSet)%>%
                                      full_join(LD$gene_info,.,by=c("set","gene_list")) %>% 
                                      dplyr::filter(!is.na(set)) %>% 
                                      dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                      dplyr::mutate(mycol=if_else(is.na(plot_set.y),mycol.x,mycol.y)) %>% 
                                      dplyr::select(-plot_set.y,-plot_set.x,-mycol.x,-mycol.y) %>%
                                      distinct() 
                                    Apply_Cluster_Math <-
                                      ApplyMath(
                                        list_data_frame,
                                        input$myMathcluster,
                                        input$radioplotnromcluster,
                                        as.numeric(input$selectplotBinNormcluster)
                                      )
                                  }
                                })
                   if (!is_empty(list_data_frame)) {
                     Plot_Cluster_Options <-
                       MakePlotOptionFrame(LD)
                     Y_Axis_Cluster_numbers <-
                       MyXSetValues(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         log_2 = input$checkboxlog2cluster
                       )
                     reactive_values$Plot_controler_cluster <-
                       GGplotLineDot(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         Plot_Cluster_Options,
                         Y_Axis_Cluster_numbers,
                         reactive_values$Lines_Labels_List,
                         input$checkboxsmoothcluster,reactive_values$Plot_Options_ttest,
                         input$checkboxlog2cluster,
                         isolate(
                           YAxisLable(
                             input$myMathcluster,
                             input$radioplotnromcluster,
                             as.numeric(input$selectplotBinNorm),
                             input$checkboxsmoothcluster,
                             input$checkboxlog2cluster
                           )
                         ),
                         input$sliderplotOccupancy
                       )
                   }
                   
                   ol <- input$selectclusterfile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       input$pickerclusterfile
                   } else {
                     reactive_values$pickerfile_controler <- ""
                   }
                   updateSelectInput(
                     session,
                     "selectclusterfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxcluster3 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Cluster 3",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  observeEvent(input$cluster4table_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <-
                   paste0(reactive_values$clustergroups,
                          "4\nn = ",
                          length(input$cluster4table_rows_all))
                 oldnamenum <-
                   grep(paste0(reactive_values$clustergroups, "4\nn ="),
                        names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldnamenum]) {
                   print("cluster4 filter $use")
                   oldname <- names(LIST_DATA$gene_file)[oldnamenum]
                   names(LIST_DATA$gene_file)[oldnamenum] <<- newname
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = LIST_DATA$gene_file[[newname]]$full$gene[input$cluster4table_rows_all])
                   LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
                     dplyr::mutate(gene_list=if_else(gene_list == oldname,newname,gene_list))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  LD <- LIST_DATA
                                  LD$gene_info <- LD$gene_info %>%
                                    dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Group_|^Cluster_") &
                                                                  set == input$pickerclusterfile, set, "0"))
                                  list_data_frame <- Active_list_data(LD)
                                  if (!is_empty(list_data_frame)) {
                                    ListColorSet <- brewer.pal(4,"Dark2")[1:n_distinct(distinct(list_data_frame,set,gene_list,plot_set))]
                                    LD$gene_info <- list_data_frame %>% 
                                      distinct(set,gene_list,plot_set) %>% 
                                      dplyr::mutate(mycol=ListColorSet)%>%
                                      full_join(LD$gene_info,.,by=c("set","gene_list")) %>% 
                                      dplyr::filter(!is.na(set)) %>% 
                                      dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                      dplyr::mutate(mycol=if_else(is.na(plot_set.y),mycol.x,mycol.y)) %>% 
                                      dplyr::select(-plot_set.y,-plot_set.x,-mycol.x,-mycol.y) %>%
                                      distinct() 
                                    Apply_Cluster_Math <-
                                      ApplyMath(
                                        list_data_frame,
                                        input$myMathcluster,
                                        input$radioplotnromcluster,
                                        as.numeric(input$selectplotBinNormcluster)
                                      )
                                  }
                                })
                   if (!is_empty(list_data_frame)) {
                     reactive_values$Plot_Cluster_Options <-
                       MakePlotOptionFrame(LD)
                     Y_Axis_Cluster_numbers <-
                       MyXSetValues(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         log_2 = input$checkboxlog2cluster
                       )
                     reactive_values$Plot_controler_cluster <-
                       GGplotLineDot(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         reactive_values$Plot_Cluster_Options,
                         Y_Axis_Cluster_numbers,
                         reactive_values$Lines_Labels_List,
                         input$checkboxsmoothcluster,reactive_values$Plot_Options_ttest,
                         input$checkboxlog2cluster,
                         isolate(
                           YAxisLable(
                             input$myMathcluster,
                             input$radioplotnromcluster,
                             as.numeric(input$selectplotBinNorm),
                             input$checkboxsmoothcluster,
                             input$checkboxlog2cluster
                           )
                         ),
                         input$sliderplotOccupancy
                       )
                   }
                   
                   ol <- input$selectclusterfile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <- newname
                     reactive_values$pickerfile_controler <-
                       input$pickerclusterfile
                   } else {
                     reactive_values$pickerfile_controler <- ""
                   }
                   updateSelectInput(
                     session,
                     "selectclusterfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   output$valueboxcluster4 <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Cluster 4",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  # Cluster tool action ----
  observeEvent(input$actionclustertool, ignoreInit = TRUE, {
    print("cluster tool action")
    shinyjs::hide('plotcluster')
    shinyjs::hide("cluster1table")
    shinyjs::hide("cluster2table")
    shinyjs::hide("cluster3table")
    shinyjs::hide("cluster4table")
    reactive_values$clustergroups <- NULL
    if (n_distinct(LIST_DATA$gene_file[[input$selectclusterfile]]$use) < 4) {
      reactive_values$clustergroups <- "Cluster_"
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
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
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
      reactive_values$clustergroups <- "Cluster_"
    } else {
      return()
    }
  })
  
  # Group tool action ----
  observeEvent(input$actiongroupstool, ignoreInit = TRUE, {
    print("group tool action")
    shinyjs::hide('plotcluster')
    shinyjs::hide("cluster1table")
    shinyjs::hide("cluster2table")
    shinyjs::hide("cluster3table")
    shinyjs::hide("cluster4table")
    reactive_values$clustergroups <- NULL
    if (n_distinct(LIST_DATA$gene_file[[input$selectclusterfile]]$use) < 4) {
      reactive_values$clustergroups <- "Group_"
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     FindGroups(
                       LIST_DATA,
                       input$selectclusterfile,
                       input$pickerclusterfile,
                       input$sliderbincluster[1],
                       input$sliderbincluster[2]
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
      reactive_values$clustergroups <- "Group_"
    } else {
      return()
    }
  })
  
  # Cluster tool numbers ----
  observeEvent(c(input$selectclusternumber, reactive_values$clustergroups),
               ignoreInit = TRUE,
               {
                 print("cluster tool number")
                 if (is.null(reactive_values$clustergroups)) {
                   return()
                 }
                 shinyjs::hide('actionclusterdatatable')
                 shinyjs::hide('actionclusterplot')
                 withProgress(message = 'Calculation in progress',
                              detail = 'This may take a while...',
                              value = 0,
                              {
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
                   shinyjs::show('actionclusterdatatable')
                   shinyjs::show('actionclusterplot')
                   shinyjs::hide("cluster1table")
                   shinyjs::hide("cluster2table")
                   shinyjs::hide("cluster3table")
                   shinyjs::hide("cluster4table")
                   ol <- input$selectclusterfile
                   if (!ol %in% names(LIST_DATA$gene_file)) {
                     ol <-
                       grep(strsplit(ol, "\nn")[[1]][1],
                            names(LIST_DATA$gene_file),
                            value = TRUE)
                     reactive_values$pickerfile_controler <-
                       input$pickerclusterfile
                   } else {
                     reactive_values$pickerfile_controler <- ""
                   }
                   updateSelectInput(
                     session,
                     "selectclusterfile",
                     choices = names(LIST_DATA$gene_file),
                     selected = ol
                   )
                   if (!is.null(reactive_values$clustergroups) &
                       any(grep(
                         paste0(reactive_values$clustergroups, "1\nn ="),
                         names(LIST_DATA$gene_file)
                       ) > 0)) {
                     output$valueboxcluster1 <- renderValueBox({
                       valueBox(
                         n_distinct(LIST_DATA$gene_file[[grep(
                           paste0(reactive_values$clustergroups, "1\nn ="),
                           names(LIST_DATA$gene_file)
                         )]]$use),
                         "Group 1",
                         icon = icon("list"),
                         color = "green"
                       )
                     })
                   } else{
                     output$valueboxcluster1 <- renderValueBox({
                       valueBox(0,
                                "Group 1",
                                icon = icon("list"),
                                color = "green")
                     })
                   }
                   if (!is.null(reactive_values$clustergroups) &
                       any(grep(
                         paste0(reactive_values$clustergroups, "2\nn ="),
                         names(LIST_DATA$gene_file)
                       ) > 0)) {
                     output$valueboxcluster2 <- renderValueBox({
                       valueBox(
                         n_distinct(LIST_DATA$gene_file[[grep(
                           paste0(reactive_values$clustergroups, "2\nn ="),
                           names(LIST_DATA$gene_file)
                         )]]$use),
                         "Group 2",
                         icon = icon("list"),
                         color = "green"
                       )
                     })
                   } else{
                     output$valueboxcluster2 <- renderValueBox({
                       valueBox(0,
                                "Group 2",
                                icon = icon("list"),
                                color = "green")
                     })
                   }
                   if (!is.null(reactive_values$clustergroups) &
                       any(grep(
                         paste0(reactive_values$clustergroups, "3\nn ="),
                         names(LIST_DATA$gene_file)
                       ) > 0)) {
                     output$valueboxcluster3 <- renderValueBox({
                       valueBox(
                         n_distinct(LIST_DATA$gene_file[[grep(
                           paste0(reactive_values$clustergroups, "3\nn ="),
                           names(LIST_DATA$gene_file)
                         )]]$use),
                         "Group 3",
                         icon = icon("list"),
                         color = "green"
                       )
                     })
                   } else{
                     output$valueboxcluster3 <- renderValueBox({
                       valueBox(0,
                                "Group 3",
                                icon = icon("list"),
                                color = "green")
                     })
                   }
                   if (!is.null(reactive_values$clustergroups) &
                       any(grep(
                         paste0(reactive_values$clustergroups, "4\nn ="),
                         names(LIST_DATA$gene_file)
                       ) > 0)) {
                     output$valueboxcluster4 <- renderValueBox({
                       valueBox(
                         n_distinct(LIST_DATA$gene_file[[grep(
                           paste0(reactive_values$clustergroups, "4\nn ="),
                           names(LIST_DATA$gene_file)
                         )]]$use),
                         "Group 4",
                         icon = icon("list"),
                         color = "green"
                       )
                     })
                   } else{
                     output$valueboxcluster4 <- renderValueBox({
                       valueBox(0,
                                "Group 4",
                                icon = icon("list"),
                                color = "green")
                     })
                   }
                   LD$gene_info <- LD$gene_info %>%
                     dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Group_|^Cluster_") &
                                                   set == input$pickerclusterfile, set, "0"))
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                {
                                  list_data_frame <- Active_list_data(LD)
                                  if (!is_empty(list_data_frame)) {
                                    ListColorSet <- brewer.pal(4,"Dark2")[1:n_distinct(distinct(list_data_frame,set,gene_list,plot_set))]
                                    LD$gene_info <- list_data_frame %>% 
                                      distinct(set,gene_list,plot_set) %>% 
                                      dplyr::mutate(mycol=ListColorSet)%>%
                                      full_join(LD$gene_info,.,by=c("set","gene_list")) %>% 
                                      dplyr::filter(!is.na(set)) %>% 
                                      dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                                      dplyr::mutate(mycol=if_else(is.na(plot_set.y),mycol.x,mycol.y)) %>% 
                                      dplyr::select(-plot_set.y,-plot_set.x,-mycol.x,-mycol.y) %>%
                                      distinct() 
                                    Apply_Cluster_Math <- ApplyMath(
                                      list_data_frame,
                                      input$myMathcluster,
                                      input$radioplotnromcluster,
                                      as.numeric(input$selectplotBinNormcluster)
                                    )
                                  }
                                })
                   
                   if (!is_empty(list_data_frame)) {
                     Plot_Cluster_Options <-
                       MakePlotOptionFrame(LD)
                     Y_Axis_Cluster_numbers <-
                       MyXSetValues(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         log_2 = input$checkboxlog2cluster
                       )
                     reactive_values$Plot_controler_cluster <-
                       GGplotLineDot(
                         Apply_Cluster_Math,
                         input$sliderplotBinRange,
                         Plot_Cluster_Options,
                         Y_Axis_Cluster_numbers,
                         reactive_values$Lines_Labels_List,
                         input$checkboxsmoothcluster,reactive_values$Plot_Options_ttest,
                         input$checkboxlog2cluster,
                         isolate(
                           YAxisLable(
                             input$myMathcluster,
                             input$radioplotnromcluster,
                             as.numeric(input$selectplotBinNorm),
                             input$checkboxsmoothcluster,
                             input$checkboxlog2cluster
                           )
                         ),
                         input$sliderplotOccupancy
                       )
                   }
                   shinyjs::show('plotcluster')
                 } else {
                   output$valueboxcluster1 <- renderValueBox({
                     valueBox(0,
                              "Cluster 1",
                              icon = icon("list"),
                              color = "green")
                   })
                   output$valueboxcluster2 <- renderValueBox({
                     valueBox(0,
                              "Cluster 2",
                              icon = icon("list"),
                              color = "green")
                   })
                   output$valueboxcluster3 <- renderValueBox({
                     valueBox(0,
                              "Cluster 3",
                              icon = icon("list"),
                              color = "green")
                   })
                   output$valueboxcluster4 <- renderValueBox({
                     valueBox(0,
                              "Cluster 4",
                              icon = icon("list"),
                              color = "green")
                   })
                   return()
                 }
               })
  
  # Create and Show cluster data table ----
  observeEvent(input$actionclusterdatatable, ignoreInit = TRUE, {
    updateTabItems(session, "clustertooltab", "Cluster 1")
    newnames <- gsub("(.{20})", "\\1... ", input$pickerclusterfile)
    if (!is.null(reactive_values$clustergroups) & any(grep(
      paste0(reactive_values$clustergroups, "1\nn ="),
      names(LIST_DATA$gene_file)
    ) > 0)) {
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$cluster1table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "1\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$full,
                           rownames = FALSE,
                           colnames = gsub("(.{22})", "\\1\n", newnames),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "1\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show("cluster1table")
    } else {
      output$cluster1table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            options = list(searching = FALSE)
          )
        )
    }
    
    if (!is.null(reactive_values$clustergroups) & any(grep(
      paste0(reactive_values$clustergroups, "2\nn ="),
      names(LIST_DATA$gene_file)
    ) > 0)) {
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$cluster2table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "2\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$full,
                           rownames = FALSE,
                           colnames = gsub("(.{22})", "\\1\n", newnames),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "2\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show("cluster2table")
    } else {
      output$cluster2table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            options = list(searching = FALSE)
          )
        )
    }
    
    if (!is.null(reactive_values$clustergroups) & any(grep(
      paste0(reactive_values$clustergroups, "3\nn ="),
      names(LIST_DATA$gene_file)
    ) > 0)) {
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$cluster3table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "3\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$full,
                           rownames = FALSE,
                           colnames = gsub("(.{22})", "\\1\n", newnames),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "3\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show("cluster3table")
    } else {
      output$cluster3table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            options = list(searching = FALSE)
          )
        )
    }
    
    if (!is.null(reactive_values$clustergroups) & any(grep(
      paste0(reactive_values$clustergroups, "4\nn ="),
      names(LIST_DATA$gene_file)
    ) > 0)) {
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$cluster4table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "4\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$full,
                           rownames = FALSE,
                           colnames = gsub("(.{22})", "\\1\n", newnames),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "4\nn ="),
                             names(LIST_DATA$gene_file)
                           )]]$info,
                           options = list(
                             pageLength = 15,
                             scrollX = TRUE,
                             scrollY = TRUE,
                             autoWidth = FALSE,
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
                   })
      shinyjs::show("cluster4table")
    } else {
      output$cluster4table <-
        DT::renderDataTable(
          datatable(
            LIST_DATA$gene_file[[1]]$empty,
            rownames = FALSE,
            colnames = strtrim(newnames, 24),
            options = list(searching = FALSE)
          )
        )
    }
    shinyjs::hide('actionclusterdatatable')
  })
  
  # creat and show cluster plot ----
  observeEvent(input$actionclusterplot, {
    print("cluster plot button")
    shinyjs::show('plotcluster')
    if (!is.null(reactive_values$clustergroups) & any(grep(
      paste0(reactive_values$clustergroups, "1\nn ="),
      names(LIST_DATA$gene_file)
    ) > 0)) {
      LD <- LIST_DATA
      LD$gene_info <- LD$gene_info %>%
        dplyr::mutate(onoff=if_else(str_detect(gene_list,"^Group_|^Cluster_") &
                                      set == input$pickerclusterfile, set, "0"))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     list_data_frame <- Active_list_data(LD)
                     if (!is_empty(list_data_frame)) {
                       ListColorSet <- brewer.pal(4,"Dark2")[1:n_distinct(distinct(list_data_frame,set,gene_list,plot_set))]
                       LD$gene_info <- list_data_frame %>% 
                         distinct(set,gene_list,plot_set) %>% 
                         dplyr::mutate(mycol=ListColorSet)%>%
                         full_join(LD$gene_info,.,by=c("set","gene_list")) %>% 
                         dplyr::filter(!is.na(set)) %>% 
                         dplyr::mutate(plot_set=if_else(is.na(plot_set.y),plot_set.x,plot_set.y)) %>% 
                         dplyr::mutate(mycol=if_else(is.na(plot_set.y),mycol.x,mycol.y)) %>% 
                         dplyr::select(-plot_set.y,-plot_set.x,-mycol.x,-mycol.y) %>%
                         distinct() 
                       Apply_Cluster_Math <- ApplyMath(
                         list_data_frame,
                         input$myMathcluster,
                         input$radioplotnromcluster,
                         as.numeric(input$selectplotBinNormcluster)
                       )
                     }
                   })
      if (!is_empty(list_data_frame)) {
        reactive_values$Plot_Cluster_Options <-
          MakePlotOptionFrame(LD)
        Y_Axis_Cluster_numbers <-
          MyXSetValues(
            Apply_Cluster_Math,
            input$sliderplotBinRange,
            log_2 = input$checkboxlog2cluster
          )
        reactive_values$Plot_controler_cluster <- GGplotLineDot(
          Apply_Cluster_Math,
          input$sliderplotBinRange,
          reactive_values$Plot_Cluster_Options,
          Y_Axis_Cluster_numbers,
          reactive_values$Lines_Labels_List,
          input$checkboxsmoothcluster,reactive_values$Plot_Options_ttest,
          input$checkboxlog2cluster,
          isolate(
            YAxisLable(
              input$myMathcluster,
              input$radioplotnromcluster,
              as.numeric(input$selectplotBinNormcluster),
              input$checkboxsmoothcluster,
              input$checkboxlog2cluster
            )
          ),
          input$sliderplotOccupancy
        )
      }
    }
  })
  
  # CDF generate gene list ----
  observeEvent(input$actioncdfdatatable, ignoreInit = TRUE, {
    if (any(grep("CDF ", names(LIST_DATA$gene_file)) > 0)) {
      newnames1 <-
        gsub("\n", " ",
             grep("CDF ",
                  names(LIST_DATA$gene_file),
                  value = TRUE))
      df <-
        dplyr::select(LIST_DATA$gene_file[[grep("CDF ", names(LIST_DATA$gene_file))]]$full, -bin, -set) %>%
        dplyr::mutate(value=round(log2(value)),5) %>% 
        tidyr::spread(., plot_set, value)
      # PI EI differenc tool
      if (length(names(df)) == 3) {
        df[paste("\'",
                 names(df[2]),
                 "\'",
                 " By ",
                 "\'",
                 names(df[3]),
                 "\'",
                 sep = "")] <- df[2] - df[3]
      }
      df <- arrange(df, df[[names(df)[2]]])
      dt <- datatable(
        df,
        colnames = gsub("(.{22})", "\\1\n", c(newnames1, names(df)[-1])),
        rownames = FALSE,
        class = 'cell-border stripe compact',
        filter = 'top',
        caption = LIST_DATA$gene_file[[grep("CDF ", names(LIST_DATA$gene_file))]]$info,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          scrollY = TRUE,
          autoWidth = FALSE,
          width = 5,
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
    } else {
      dt <- datatable(
        LIST_DATA$gene_file[[1]]$empty,
        rownames = FALSE,
        options = list(searching = FALSE)
      )
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   output$cdftable <- DT::renderDataTable(dt)
                 })
    shinyjs::hide('actioncdfdatatable')
    shinyjs::show('cdftable')
  })
  
  
  # cdf tool gene lists $use ----
  observeEvent(input$cdftable_rows_all,
               ignoreInit = TRUE,
               ignoreNULL = TRUE,
               {
                 newname <- paste("CDF n =", length(input$cdftable_rows_all))
                 oldnamenum <- grep("CDF ", names(LIST_DATA$gene_file))
                 if (newname != names(LIST_DATA$gene_file)[oldnamenum] &
                     length(input$cdftable_rows_all) != 0) {
                   print("cdf filter $use")
                   oldname <- names(LIST_DATA$gene_file)[oldnamenum]
                   names(LIST_DATA$gene_file)[oldnamenum] <<- newname
                   dt <-
                     dplyr::select(LIST_DATA$gene_file[[newname]]$full, -bin, -set) %>%
                     tidyr::spread(., plot_set, value)
                   LIST_DATA$gene_file[[newname]]$use <<-
                     tibble(gene = dt$gene[input$cdftable_rows_all])
                   LIST_DATA$gene_info <<- LIST_DATA$gene_info %>% 
                     dplyr::mutate(gene_list=if_else(gene_list == oldname,newname,gene_list))
                   reactive_values$Picker_controler <- 
                     c(names(LIST_DATA$gene_file), distinct(LIST_DATA$gene_info, set)$set)
                   df_options <- LIST_DATA$gene_info %>%
                     dplyr::filter(gene_list ==  newname) %>%
                     dplyr::mutate(set = paste(
                       sub("\n", " ", newname),
                       gsub("(.{20})", "\\1\n", plot_set),
                       sep = '\n'
                     ))
                   if (any(duplicated(df_options$mycol))) {
                     df_options$mycol <-
                       brewer.pal(8, "Set1")[1:n_distinct(df_options$set)]
                   }
                   df <-
                     inner_join(LIST_DATA$gene_file[[newname]]$full,
                                LIST_DATA$gene_file[[newname]]$use,
                                by = "gene") %>%
                     dplyr::mutate(set = paste(
                       sub("\n", " ", newname),
                       gsub("(.{20})", "\\1\n", plot_set),
                       sep = '\n'
                     ))
                   use_header <-
                     pull(distinct(df_options, myheader))
                   if (n_groups(group_by(df_options, set)) == 2 &
                       n_distinct(df$gene) > 1) {
                     tt_name <- pull(distinct(df_options, set))
                     tt <-
                       suppressWarnings(ks.test(pull(dplyr::filter(
                         df, set == tt_name[1]
                       ), value),
                       pull(dplyr::filter(
                         df, set == tt_name[2]
                       ), value)))
                     if (tt[[2]] == 0) {
                       use_header <- paste(use_header, "  p-value < 2.2e-16 ")
                     } else {
                       use_header <-
                         paste(use_header, paste("  p-value = ", format(tt[[2]], scientific = TRUE)))
                     }
                   }
                   mycdf <- GGplotC(df, df_options, use_header)
                   output$plotcdf <- renderPlot({
                     mycdf
                   })
                   output$valueboxcdf <- renderValueBox({
                     valueBox(
                       n_distinct(LIST_DATA$gene_file[[newname]]$use),
                       "Gene List 1",
                       icon = icon("list"),
                       color = "green"
                     )
                   })
                 }
               })
  
  # CDF tool action ----
  observeEvent(input$actioncdftool, ignoreInit = TRUE, {
    print("CDF tool action")
    shinyjs::hide('cdftable')
    shinyjs::hide('plotcdf')
    if (any(between(
      input$sliderbincdf1,
      input$sliderbincdf2[1],
      input$sliderbincdf2[2]
    )) |
    any(between(
      input$sliderbincdf2,
      input$sliderbincdf1[1],
      input$sliderbincdf1[2]
    ))) {
      showModal(modalDialog(
        title = "Information message",
        paste("Bins regions should not overlab, \nBins set to 1/3 2/3"),
        size = "s",
        easyClose = TRUE
      ))
      updateSliderInput(session,
                        "sliderbincdf1",
                        value = c(
                          LIST_DATA$x_plot_range[1],
                          floor(LIST_DATA$x_plot_range[2] / 4)
                        ))
      updateSliderInput(session,
                        "sliderbincdf2",
                        value = c(
                          ceiling(LIST_DATA$x_plot_range[2] / 4) + 1,
                          LIST_DATA$x_plot_range[2]
                        ))
    }
    ttt <-
      reactiveValuesToList(input)[gsub(" ", "-cdfspace2-", gsub("\n", "-cdfspace1-", names(LIST_DATA$gene_file)))]
    checkboxonoff <- list()
    for (i in names(ttt)) {
      for (tt in ttt[i]) {
        selectgenelistonoff <-
          gsub("-cdfspace2-", " ", gsub("-cdfspace1-", "\n", i))
        checkboxonoff[[selectgenelistonoff]] <-
          c(checkboxonoff[[selectgenelistonoff]], tt)
      }
    }
    if (is_empty(checkboxonoff)) {
      return()
    }
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <-
                     CumulativeDistribution(
                       LIST_DATA,
                       checkboxonoff,
                       input$sliderbincdf1[1],
                       input$sliderbincdf1[2],
                       input$sliderbincdf2[1],
                       input$sliderbincdf2[2]
                     )
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      if(LIST_DATA$STATE[2] == 0){
        LIST_DATA$STATE[2] <<- -2
      }
      shinyjs::show('actioncdfdatatable')
      shinyjs::show('plotcdf')
      newname <-
        grep("CDF ", names(LIST_DATA$gene_file), value = TRUE)
      df_options <-
        LIST_DATA$gene_info %>%
        dplyr::filter(gene_list ==  newname) %>%
        dplyr::mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{20})", "\\1\n", plot_set),
          sep = '\n'
        ))
      # fix same color problems
      if (any(duplicated(df_options$mycol))) {
        df_options$mycol <-
          brewer.pal(8, "Set1")[1:n_distinct(df_options$set)]
      }
      df <- LIST_DATA$gene_file[[newname]]$full %>%
        dplyr::mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{20})", "\\1\n", plot_set),
          sep = '\n'
        ))
      use_header <- pull(distinct(df_options, myheader))
      if (n_groups(group_by(df_options, set)) == 2 &
          n_distinct(df$gene) > 1) {
        tt_name <- pull(distinct(df_options, set))
        tt <-
          suppressWarnings(ks.test(pull(dplyr::filter(
            df, set == tt_name[1]
          ), value),
          pull(dplyr::filter(
            df, set == tt_name[2]
          ), value)))
        if (tt[[2]] == 0) {
          use_header <- paste(use_header, "  p-value < 2.2e-16 ")
        } else {
          use_header <-
            paste(use_header, paste("  p-value = ", format(tt[[2]], scientific = TRUE)))
        }
      }
      mycdf <- GGplotC(df, df_options, use_header)
      output$plotcdf <- renderPlot({
        mycdf
      })
      output$valueboxcdf <- renderValueBox({
        valueBox(
          n_distinct(df$gene),
          "Gene List",
          icon = icon("list"),
          color = "green"
        )
      })
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
      tags$head(tags$style(
        HTML(
          ".shiny-notification {
          position:fixed;
          top: calc(50%);;
          left: calc(50%);;
          }
          "
        )
      )),
      menuItem("Load Data", tabName = "loaddata", icon = icon("file")),
      menuItem("Norm data", tabName = "filenorm", icon = icon("files-o")),
      menuItem("Plot", tabName = "mainplot", icon = icon("area-chart")),
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
      menuItem("Filter Tool", tabName = "sorttool", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showsorttoolpicker",
          div(
            style = "margin-bottom: -30px;",
            sliderInput(
              "slidersortbinrange",
              label = "Select Bin Range:",
              min = 0,
              max = 80,
              value = c(0, 80)
            )),
          div(
            style = "margin-bottom: -20px;",
            awesomeCheckbox("checkboxfilterall","Filter if any",value = TRUE)
          ),
          div(
            style = "margin-bottom: -30px;",
            selectInput(
              inputId = "selectsortfile",
              label = "Select gene list to sort on",
              choices = "Load data file",
              width = "99%"
            )),
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
            label = "Select file",
            choices = "Load data file",
            multiple = F,
            options = list(title = "Select file")
          )
        )
      ),
      menuItem("CDF Tools", tabName = "cdftool", icon = icon("gears")),
      hidden(
        div(
          style = "padding-left: 15%;",
          id = "showcdftoolpicker",
          sliderInput(
            "sliderbincdf1",
            label = "Select numerator Bin Range:",
            min = 0,
            max = 80,
            value = c(0, 0)
          ),
          sliderInput(
            "sliderbincdf2",
            label = "Select denominator Bin Range:",
            min = 0,
            max = 80,
            value = c(0, 80)
          ),
          actionButton("actioncdftool", "Plot CDF")
        )
      ),
      
      hr(style = "color: #fff; background-color: #337ab7; border-color: #2e6da4;")
    )
  ),
  dashboardBody(
    useShinyjs(),
    tabItems(
      # load data tab ----
      tabItem(tabName = "loaddata",
              fluidRow(
                tags$style(
                  ".nav-tabs-custom .nav-tabs li.active a { background-color: transparent; border-color: #2e6da4; } "
                ),
                tabBox(
                  id = "loadfiles",
                  width = 5,
                  tabPanel("Table",
                           div(
                             style = "height: 200px;",
                             fluidRow(div(
                               style = "padding: 5px 15px;",
                               fileInput(
                                 "filetable",
                                 label = "Load table/URL.txt file",
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
                                                 label = "Load gene list",
                                                 accept = c('.txt'))
                                       
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
                               ))),
                  tabPanel(title = "Save",
                           div(style = "height: 200px;",
                               fluidRow(
                                 hidden(div(
                                   id = "hiddensave",
                                   style = "padding: 5px 15px;",
                                   radioGroupButtons(
                                     inputId = "radiogroupsave",
                                     choices = c("Save Gene list", 
                                                 "Save Gene list as bed",
                                                 "Save full Table file", 
                                                 "Save Compleat color - file pair"),
                                     selected = "Save Gene list",
                                     checkIcon = list(
                                       yes = tags$i(class = "fa fa-check-square", 
                                                    style = "color: steelblue"),
                                       no = tags$i(class = "fa fa-square-o", 
                                                   style = "color: steelblue"))
                                   ),
                                   downloadButton("downloadGeneList", "Save List")
                                 ))
                               )))
                ),
                
                hidden(div(
                  id = "startoff",
                  box(
                    title =  "Select Gene list",
                    width = 7,
                    height = "300px",
                    status = "primary",
                    solidHeader = T,
                    selectInput("selectgenelistoptions", "", choices = "Compleat"),
                    actionButton("actionremovegene", "Remove Gene list"),
                    fluidRow(
                      column(4, selectInput(
                        "kbrewer",
                        "color brewer theme",
                        c(choices = "select", kBrewerList)
                      )),
                      column(3, style = "padding-top:5%;",
                             actionButton("BttnNewColor", "Set color same as Compleat"))
                    )
                  ),
                  box(
                    title = "File Options",
                    solidHeader = T,
                    width = 12,
                    box(
                      title =  "Set Plot Color Options",
                      width = 4,
                      status = "primary",
                      solidHeader = T,
                      fluidRow(
                        box(
                          width = 12,
                          status = "info",
                          background = "light-blue",
                          colourInput("colourhex", "Select color HEX"),
                          tags$hr(),
                          textInput("textrgbtohex", "RGB"),
                          actionButton("actionmyrgb", "Update HEX color")
                        )
                      ),
                      selectInput("selectdot", "Select dot type", choices = kDotOptions),
                      selectInput("selectline", "Select line type", choices = kLineOptions)
                    ),
                    box(
                      title =  "Set Plot Options",
                      width = 8,
                      status = "primary",
                      solidHeader = T,
                      selectInput("selectdataoption", "", choices = "Load data file"),
                      fluidRow(column(
                        4, actionButton("actionremovefile", "Remove File(s)")
                      ),
                      column(
                        4,
                        awesomeCheckbox("checkboxremovefile",
                                      "remove all files and restart", value = FALSE)
                      )
                      ),
                      tags$hr(style = "color: #2e6da4; background-color: #2e6da4; border-color: #2e6da4;"),
                      textInput("textnickname", "Update Nickname"),
                      actionButton("actionoptions", "Set Nickname"),
                      helpText("Need to press to update")
                    )
                  )
                ))
              )),
      # Norm file tab ----
      tabItem(
        tabName = "filenorm",
        
        box(
          title = "Select files for normalization",
          width = 12,
          status = "primary",
          solidHeader = T,
          div(style = "padding-left: 15%;",
              fluidRow(
                pickerInput(
                  "pickernumerator",
                  label = "numerator",
                  width = "90%",
                  choices = "Load data file",
                  multiple = F,
                  options = list(title = "Select first file")
                )
              )),
          div(style = "padding-left: 15%;",
              fluidRow(
                column(
                  3,
                  radioGroupButtons(
                    "adddata",
                    label = "",
                    status = "primary",
                    choices = c("/", "+"),
                    selected = "/"
                  )
                ),
                column(4, style = "padding-top: 4%;",
                       actionButton("actionnorm", label = "create file"))
              )),
          div(style = "padding-left: 15%;",
              fluidRow(
                pickerInput(
                  "pickerdenominator",
                  label = "denominator",
                  width = "90%",
                  choices = "Load data file",
                  multiple = F,
                  options = list(title = "Select second file")
                )
              ),
              fluidRow(
                textInput("textnromname", "Norm file name",
                          width = "90%",))
          ),
          awesomeRadio(
            "radiogenebygene",
            label = "",
            choices = c("bin by bin", "mean of bins by mean of bins"),
            selected = "bin by bin"
          ),
          awesomeCheckbox(
            "checkboxnormzero",
            label = "replace 0 with min/2",value = FALSE),
          valueBoxOutput("valueboxnormfile")
        )
      ),
      
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
              )
              ),
              div(id = "showmainplot",  fluidRow(
                tags$style(
                  ".nav-tabs-custom .nav-tabs li.active a { background-color: transparent; border-color: #2e6da4; } "
                ),
                tabBox(
                  id = "tabboxmainplot",
                  width = 12, height = "500px",
                  tabPanel("Gene lists",
                           div(
                             box(title = "Main",
                                 width = 6,
                                 status = "primary",
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 uiOutput("DynamicGenePicker_main")
                             )),
                           hidden(
                             div(
                               id = "showpickersort",
                               box(
                                 title = "Filter (max 4 lists)",
                                 width = 6,
                                 status = "primary",
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 uiOutput("DynamicGenePicker_sort")
                               )
                             )),
                           hidden(
                             div(
                               id = "showpickercomparisons",
                               box(
                                 title = "Gene comparisons",
                                 width = 6,
                                 status = "primary",
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 uiOutput("DynamicGenePicker_comparisons")
                               )
                             )),
                           hidden(
                             div(
                               id = "showpickerratio",
                               box(
                                 title = "Ratio",
                                 width = 6,
                                 status = "primary",
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 uiOutput("DynamicGenePicker_ratio")
                               )
                             )),
                           hidden(
                             div(
                               id = "showpickercluster",
                               box(
                                 title = "Clusters/Groups",
                                 width = 6,
                                 status = "primary",
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 uiOutput("DynamicGenePicker_clusters")
                               )
                             )),
                           hidden(
                             div(
                               id = "showpickercdf",
                               box(
                                 title = "CDF",
                                 width = 6,
                                 status = "primary",
                                 solidHeader = T,
                                 collapsible = T,
                                 collapsed = F,
                                 uiOutput("DynamicGenePicker_cdf")
                               )
                             ))
                  ),
                  tabPanel("Plot Options",
                           box(
                             title = "Sliders",
                             status = "primary",
                             solidHeader = T,
                             width = 6,
                             collapsible = F,
                             collapsed = F,
                             sliderInput(
                               "sliderplotBinRange",
                               label = "Plot Bin Range:",
                               min = 0,
                               max = 80,
                               value = c(0, 80)
                             ),
                             fluidRow(
                               column(
                                 4,
                                 numericInput("numericYRangeLow", label = "Plot Y min:", value = 0)
                               ),
                               column(
                                 4,
                                 numericInput("numericYRangeHigh", label = "Plot Y max:", value = 0)
                               ),
                               column(
                                 2,
                                 style = "padding-top: 8%;",
                                 actionBttn(
                                   inputId = "actionButtonYXrange",
                                   label = "apply",
                                   style = "unite",
                                   color = "default",
                                   size = "sm"
                                 )
                               )
                             )
                           ),
                           box(
                             title = "Math and Normalization",
                             status = "primary",
                             solidHeader = T,
                             width = 6,
                             collapsible = F,
                             collapsed = F,
                             column(
                               4,
                               selectInput(
                                 "selectplotBinNorm",
                                 label = "Bin Norm:",
                                 choices = c(0:80),
                                 selected = 0
                               )
                             ),
                             column(
                               8,
                               selectInput(
                                 "selectplotnrom",
                                 label = "Set Y Normalization",
                                 choices = c("none", "relative frequency", "rel gene frequency"),
                                 selected = "none"
                               )
                             ),
                             column(4, awesomeCheckbox("checkboxsmooth", label = "smooth"), 
                                    awesomeCheckbox("checkboxlog2", label = "log2")),
                             column(8, 
                                    selectInput("myMath",
                                                label = " ",
                                                choices = c("mean", "sum", "median", "var"),
                                                selected = "mean"
                                    ))
                           )),
                  tabPanel("Lines and Labels",
                           box(
                             width = 12,
                             status = "primary",
                             solidHeader = T,
                             collapsible = FALSE,
                             collapsed = FALSE,
                             div(
                               style = "padding-left: 25px; display:inline-block;",
                               selectInput(
                                 selectize = T,
                                 "selectlineslabels",
                                 width = "200px",
                                 label = "quick set lines and labels",
                                 choices = c("Choose one" = "",
                                             kLinesandlabels)
                               )
                             ),
                             column(12,
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
                                      textInput("numerictssname", value = "TSS", label = "lable",width = "50px"),
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
                                        "numerictes",
                                        "pA bin",
                                        value = 45,
                                        min = 0,
                                        max = 100
                                      )
                                    ),
                                    div(
                                      style = "padding:2px; display:inline-block;",
                                      textInput("numerictesname", value = "pA", label = "lable",width = "50px"),
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
                                      style = "padding:2px 8px 2px 2px; display:inline-block;",
                                      numericInput(
                                        "numericlabelspaceing",
                                        "every bin",
                                        value = 5,
                                        min = 0,
                                        max = 100
                                      )
                                    ),
                                    actionButton("actionlineslabels", "UPDATE PLOT")
                             ),
                             helpText("For 543 style 0 > TSS < 5|4 < 4|3 < pA < max bin"),
                             div(
                               textInput("landlnames", "", label = "Yaxis labels"),
                               textInput("landlposition", "", label = "Yaxis lable position (numbers only)")
                             ),
                             helpText("select buttons for more options"),
                             column(
                               12,
                               div(
                                 style = "padding-left: -5px; display:inline-block;",
                                 dropdownButton(
                                   tags$h3("Set TSS Options"),
                                   
                                   selectInput(
                                     inputId = 'selecttsscolor',
                                     label = 'TSS line and lable color',
                                     choices = c("red", "green", "blue", "brown", "black", "white"),
                                     selected = "green"
                                   ),
                                   selectInput(
                                     inputId = 'selecttssline',
                                     label = 'TSS line type',
                                     choices = c("dotted", "solid"),
                                     selected = "dotted"
                                   ),
                                   icon = icon("sliders"),
                                   status = "success",
                                   tooltip = tooltipOptions(title = "TSS Options")
                                 )
                               ),
                               div(
                                 style = "padding-left: 20px; display:inline-block;",
                                 dropdownButton(
                                   tags$h3("Set 5|4 Options"),
                                   
                                   selectInput(
                                     inputId = 'selectbody1color',
                                     label = '5|4 line and lable color',
                                     choices = c("red", "green", "blue", "brown", "black", "white"),
                                     selected = "black"
                                   ),
                                   selectInput(
                                     inputId = 'selectbody1line',
                                     label = '5|4 line type',
                                     choices = c("dotted", "solid"),
                                     selected = "solid"
                                   ),
                                   icon = icon("sliders"),
                                   tooltip = tooltipOptions(title = "5|4 Options")
                                 )
                               ),
                               div(
                                 style = "padding-left: 20px; display:inline-block;",
                                 dropdownButton(
                                   tags$h3("Set 4|3 Options"),
                                   
                                   selectInput(
                                     inputId = 'selectbody2color',
                                     label = '4|3 line and lable color',
                                     choices = c("red", "green", "blue", "brown", "black", "white"),
                                     selected = "black"
                                   ),
                                   selectInput(
                                     inputId = 'selectbody2line',
                                     label = '4|3 line type',
                                     choices = c("dotted", "solid"),
                                     selected = "solid"
                                   ),
                                   icon = icon("sliders"),
                                   tooltip = tooltipOptions(title = "4|3 Options")
                                 )
                               ),
                               div(
                                 style = "padding-left: 25px; display:inline-block;",
                                 dropdownButton(
                                   tags$h3("Set TES Options"),
                                   
                                   selectInput(
                                     inputId = 'selecttescolor',
                                     label = 'TES line and lable color',
                                     choices = c("red", "green", "blue", "brown", "black", "white"),
                                     selected = "red"
                                   ),
                                   selectInput(
                                     inputId = 'selecttesline',
                                     label = 'TES line type',
                                     choices = c("dotted", "solid"),
                                     selected = "dotted"
                                   ),
                                   icon = icon("sliders"),
                                   status = "danger",
                                   tooltip = tooltipOptions(title = "TES Options")
                                 )
                               ),
                               div(
                                 style = "padding-left: 25px; display:inline-block;",
                                 dropdownButton(
                                   tags$h3("Set font Options"),
                                   
                                   numericInput(
                                     inputId = 'selectvlinesize',
                                     "Set vertcal line size",
                                     value = 2,
                                     min = .5,
                                     max = 10,
                                     step = .5
                                   ),
                                   numericInput(
                                     inputId = 'selectfontsizex',
                                     "Set X axis font size",
                                     value = 13,
                                     min = 1,
                                     max = 30,
                                     step = 1
                                   ),
                                   numericInput(
                                     inputId = 'selectfontsizey',
                                     "Set Y axis font size",
                                     value = 13,
                                     min = 1,
                                     max = 30,
                                     step = 1
                                   ),
                                   icon = icon("sliders"),
                                   status = "warning",
                                   tooltip = tooltipOptions(title = "Font Options")
                                 )
                               ),
                               div(
                                 style = "padding-left: 25px; display:inline-block;",
                                 dropdownButton(
                                   tags$h3("Set line Options"),
                                   
                                   numericInput(
                                     inputId = 'selectlinesize',
                                     "Set plot line size",
                                     value = 2.5,
                                     min = .5,
                                     max = 10,
                                     step = .5
                                   ),
                                   numericInput(
                                     inputId = 'selectlegendsize',
                                     "Set plot line size",
                                     value = 10,
                                     min = 1,
                                     max = 20,
                                     step = 1
                                   ),
                                   icon = icon("sliders"),
                                   status = "warning",
                                   tooltip = tooltipOptions(title = "Line Options")
                                 )
                               )
                             )
                           )),
                  tabPanel(title = "t.test",
                           box(
                             collapsed = F,
                             collapsible = F,
                             width = 12,
                             status = "primary",
                             solidHeader = T,
                             column(3,
                                    selectInput(inputId = "switchttest",
                                                label = "Plot t.test",
                                                choices = c("none","by files", "by lists"),
                                                selected = "none"
                                    )
                             ),
                             column(3,
                                    selectInput(inputId = "switchttesttype",
                                                label = "pick test",
                                                choices = c("t.test","ks.test", "wilcox.test"),
                                                selected = "wilcox.test")
                             ),
                             column(3,
                                    selectInput(inputId = "selectttestalt",
                                                label = "alternative",
                                                choices = c("two.sided", "less", "greater"),
                                                selected = "two.sided")
                             ),
                             column(3,
                                    selectInput(inputId = "selectttestpaired",
                                                label = "Paired",
                                                choices = c("TRUE", "FALSE"),
                                                selected = "FALSE")
                             ),
                             column(3,
                                    selectInput(inputId = "selectttestexact",
                                                label = "exact",
                                                choices = c("NULL", "TRUE", "FALSE"),
                                                selected = "FALSE")
                             ),
                             column(3,
                                    selectInput(inputId = "selectttestlog",
                                                label = "log p.value",
                                                # choices = c("none","log","-log", "log2", "log10"),
                                                choices = c("none","-log", "-log10"),
                                                selected = "-log10")
                             ),
                             column(2,
                                    selectInput("padjust",
                                                label = "p.adjust?",
                                                choices = c("NO", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                                                            "fdr", "none"),
                                                selected = "fdr")
                             ),
                             column(
                               2,
                               numericInput("numericYRangeLowpval", label = "p.value Y min:", value = 0)
                             ),
                             column(
                               2,
                               numericInput("numericYRangeHighpval", label = "p.value Y max:", value = 0)
                             ),column(
                               2,
                               actionButton("actionttest","replot")
                             ),
                             column(4,
                                    sliderInput(
                                      "sliderplotOccupancy",
                                      label = "p.value plot occupancy",
                                      min = 1,
                                      max = 3,
                                      step = 0.5,
                                      value = 1)
                             ),
                             column(6,
                                    selectInput(inputId = "selectttestitem",
                                                label = "Select to modify",
                                                choices = c("none"),
                                                selected = "none"),
                                    
                                    selectInput("selectlinettest", "Select line type", choices = c("Select", kLineOptions)),
                                    numericInput(
                                      inputId = 'selectttestlinesize',
                                      "Set plot line size",
                                      value = 2.5,
                                      min = .5,
                                      max = 10,
                                      step = .5)
                                    
                             ),
                             column(2,
                                    colourInput("selectcolorttest", "Select color HEX")
                             ),
                             column(3,
                                    numericInput("hlinettest",
                                                 "horizontal line p.val 0.05",
                                                 value = 0.05,
                                                 min = -50,
                                                 max = 50,
                                                 step = .5)
                             )
                           )))
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
                  actionButton("actiongenelists", "Compare Gene lists"),
                  helpText("Shows Intersected, Exlusive, and Total gene lists")
                ),
                box(
                  title = "Gene List Tables",
                  status = "primary",
                  solidHeader = T,
                  width = 12,
                  helpText("Needs at least 2 gene lists"),
                  actionButton("actiongenelistsdatatable", "Show gene list"),
                  tabBox(
                    id = "geneliststooltab",
                    width = 12,
                    tabPanel(
                      "Intersected Gene Lists",
                      helpText("All filtering applied to gene list usage elsewhere"),
                      DT::dataTableOutput('genelists1table')
                    ),
                    tabPanel(
                      "Total Gene Lists",
                      helpText("All filtering applied to gene list usage elsewhere"),
                      DT::dataTableOutput('genelists2table')
                    ),
                    tabPanel(
                      "Exclusive Gene Lists",
                      helpText("All filtering applied to gene list usage elsewhere"),
                      DT::dataTableOutput('genelists3table')
                    )
                  )
                ),
                fluidRow(
                  valueBoxOutput("valueboxgene1"),
                  valueBoxOutput("valueboxgene2"),
                  valueBoxOutput("valueboxgene3")
                )
              )),
      # main sort tab ----
      tabItem(
        tabName = "sorttool",
        div(
          id = "enablemainsort",
          box(
            title = "Filter by sum rank",
            status = "primary",
            solidHeader = T,
            width = 6,
            fluidRow(column(
              12,
              sliderInput(
                "slidersortpercent",
                label = "% select:",
                post = "%",
                min = 1,
                max = 100,
                value = 75
              )
            ),
            ),
            fluidRow(column(
              6,
              pickerInput(
                "selectsorttop",
                "Filter Options",
                choices = c("Top%", "Middle%", "Bottom%"),
                selected = "Middle%"
              )
            )),
            fluidRow(align="center",
                     actionButton("actionsorttool", "filter sum")
            )
          ),
          box(
            title = "filter by percentile distribution",
            status = "primary",
            solidHeader = T,
            width = 6,
            fluidRow(column(
              6,
              numericInputIcon("numericsortmin",
                               "min", 
                               value = "1",
                               max = "100", min="0",
                               step = "1",
                               icon = icon("percent")
              )),
              column(
                6,
                numericInputIcon("numericsortmax",
                                 "max", 
                                 value = "99.5",
                                 max = "100", min="0",
                                 step = "1",
                                 icon = icon("percent")
                )
              )),
            fluidRow(align="center",column(6,
                                           pickerInput(
                                             "selectsortper",
                                             "Filter Option",
                                             choices = c("min%", "between%", "max%"),
                                             selected = "min%"
                                           )
            )),
            fluidRow(align="center",
                     actionButton("actionsortper", "filter percentile")
            ),
            helpText("Hint: use 1 file to display a range of %'s")
          ),
          div(
            id = "hidesortplots",
            box(
              width = 6,
              withSpinner(plotOutput("plot1sort",height = "200px"), type = 4)
            ),
            box(
              width = 6,
              withSpinner(plotOutput("plot2sort",height = "200px"), type = 4)
            )
          ),
          div(
            id = "hidesorttable",
            box(
              title = "Filter Table",
              status = "primary",
              solidHeader = T,
              width = 12,
              actionButton("actionsortdatatable", "Show gene list"),
              helpText("All filtering applied to gene list usage elsewhere"),
              DT::dataTableOutput('sorttable')
            )
          ),
          valueBoxOutput("valueboxsort")
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
                        label = "Select numerator Bin Range:",
                        min = 1,
                        max = 80,
                        value = c(1, 1)
                      )
                    ),
                    column(
                      5,
                      sliderInput(
                        "sliderbinratio2",
                        label = "Select denominator Bin Range:",
                        min = 0,
                        max = 80,
                        value = c(0, 80)
                      )
                    )
                  ),
                  actionButton("actionratiotool", "Get fold changes"),
                  awesomeCheckbox(
                    "checkratiozero",
                    label = "replace 0 with min/2",
                    value = FALSE
                  ),
                  sliderInput(
                    "sliderRatioBinNorm",
                    label = "Bin Norm:",
                    min = 0,
                    max = 80,
                    value = 0
                  )
                ),
                box(
                  title = "box Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  collapsed = TRUE,
                  withSpinner(plotOutput("plotratio"), type = 4),
                  awesomeCheckbox("checkboxoutlierratio",
                                label = "remove outliers",
                                value = FALSE),
                  numericInput(inputId = 'textboxmaxratio',
                               "yaxis max",
                               value = 0,
                               min = 0,
                               max = 1000,
                               step = .5),
                  numericInput(inputId = 'textboxminratio',
                               "yaxis min",
                               value = 0,
                               min = 0,
                               max = 1000,
                               step = .5)
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
                      tabPanel(
                        "Up Fold Change file 1",
                        helpText("All filtering applied to gene list usage elsewhere"),
                        DT::dataTableOutput('ratio1table')
                      ),
                      tabPanel(
                        "Up Fold Change file 2",
                        helpText("All filtering applied to gene list usage elsewhere"),
                        DT::dataTableOutput('ratio2table')
                      ),
                      tabPanel(
                        "No Fold Change",
                        helpText("All filtering applied to gene list usage elsewhere"),
                        DT::dataTableOutput('ratio3table')
                      )
                    )
                  )
                ),
                fluidRow(
                  valueBoxOutput("valueboxratio1"),
                  valueBoxOutput("valueboxratio2"),
                  valueBoxOutput("valueboxratio3")
                )
              )),
      # main cluster tool tab ----
      tabItem(
        tabName = "clustertool",
        div(
          id = "enablemaincluster",
          box(
            title = "Cluster tools",
            status = "primary",
            solidHeader = T,
            width = 6,
            height = "250px",
            fluidRow(column(
              4,
              selectInput(
                inputId = "selectclusternumber",
                label = "Select number of clusters",
                choices = c(4:2),
                selected = 4,
                width = "99%"
              )
            ),
            column(
              8,
              sliderInput(
                "sliderbincluster",
                label = "Select Bin Range:",
                min = 1,
                max = 80,
                value = c(1, 80)
              )
            )),
            actionButton("actionclustertool", "Get clusters"),
            actionButton("actiongroupstool", "Get groups")
          ),
          box(
            title = "Cluster Plot Options",
            status = "primary",
            solidHeader = TRUE,
            width = 6,
            height = "250px",
            fluidRow(
              column(
                5,
                awesomeRadio(
                  "myMathcluster",
                  label =
                    " ",
                  choices = c("mean", "sum"),
                  selected = "mean"
                ),
                selectInput(
                  "selectplotBinNormcluster",
                  label = "Bin Norm:",
                  choices = c(0:80),
                  selected = 0
                ),
                actionButton("actionclusterplot", "update plot")
              ),
              awesomeRadio(
                "radioplotnromcluster",
                label = "Set Y Normalization",
                choices = c("none", "relative frequency", "rel gene frequency"),
                selected = "none"
              ),
              awesomeCheckbox("checkboxsmoothcluster", label = "smooth"),
              awesomeCheckbox("checkboxlog2cluster", label = "log2")
            )
          ),
          box(
            title = "Cluster Plot",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            withSpinner(plotOutput("plotcluster"), type = 4)
          ),
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
                tabPanel(
                  "Cluster 1",
                  helpText("All filtering applied to gene list usage elsewhere"),
                  DT::dataTableOutput('cluster1table')
                ),
                tabPanel(
                  "Cluster 2",
                  helpText("All filtering applied to gene list usage elsewhere"),
                  DT::dataTableOutput('cluster2table')
                ),
                tabPanel(
                  "Cluster 3",
                  helpText("All filtering applied to gene list usage elsewhere"),
                  DT::dataTableOutput('cluster3table')
                ),
                tabPanel(
                  "Cluster 4",
                  helpText("All filtering applied to gene list usage elsewhere"),
                  DT::dataTableOutput('cluster4table')
                )
              )
            )
          ),
          fluidRow(
            valueBoxOutput("valueboxcluster1"),
            valueBoxOutput("valueboxcluster2"),
            valueBoxOutput("valueboxcluster3"),
            valueBoxOutput("valueboxcluster4")
          )
        )
      ),
      # main CDF tool tab ----
      tabItem(
        tabName = "cdftool",
        div(
          id = "enablemaincdf",
          box(
            title = "CDF tool",
            status = "primary",
            solidHeader = T,
            collapsible = T,
            width = 12,
            box(title = "Main",
                width = 6,
                status = "primary",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicCDFPicker_main")
            ),
            box(title = "Filter",
                width = 6,
                status = "primary",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicCDFPicker_filter")
            ),
            box(title = "Gene lists",
                width = 6,
                status = "primary",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicCDFPicker_comparisons")
            ),
            box(title = "Ratio & Clusters",
                width = 6,
                status = "primary",
                solidHeader = T,
                collapsible = T,
                collapsed = F,
                uiOutput("DynamicCDFPicker_ratio"),
                uiOutput("DynamicCDFPicker_clusters")
            )
            
          ),
          box(
            title = "CDF Plot",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            withSpinner(plotOutput("plotcdf"), type = 4)
          ),
          div(
            id = "hidecdftable",
            box(
              title = "CDF Tables",
              status = "primary",
              solidHeader = T,
              collapsible = T,
              width = 12,
              actionButton("actioncdfdatatable", "Show gene list(s)"),
              helpText("All filtering applied to gene list usage elsewhere"),
              DT::dataTableOutput('cdftable')
            )
          ),
          valueBoxOutput("valueboxcdf")
        )
      )
    )
  )
)

# exicute ----
shinyApp(ui = ui, server = server)