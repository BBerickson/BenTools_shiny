# globals ----
# program for loading packages
my_packages <- function(x) {
  for (i in x) {
    #  require returns TRUE invisibly if it was able to load package
    if (!require(i , character.only = TRUE)) {
      #  If package was not able to be loaded then re-install
      install.packages(i , dependencies = TRUE)
      print(paste("installing ", i, " : please wait"))
      #  Load package after installing
      require(i , character.only = TRUE)
    }
  }
}

# run load needed pakages using my_pakages(x)
suppressPackageStartupMessages(my_packages(
  c(
    "shiny",
    "shinydashboard",
    "shinyWidgets",
    "tidyverse",
    "fastcluster",
    "shinyjs",
    "colourpicker",
    "RColorBrewer"
  )
))



LIST_DATA <- list(
  table_file = list(),
  # [[]] gene X1 X2 ...
  gene_file = list(),
  # holds $common genes from files and $gene file(s)
  gene_info = list(),
  # for holding gene file info in a list of lists, a set for $common and each $gene file(s) [c("set", "dot", "line", "color", plot?, nrom)]
  clust = list(), # Cluster holder
  x_plot_range = c(0, 0),
  STATE = c(0, "common", 0, 0) # flow control, gene list flow control, none, first time plot tab
)      

# types of dots to be used in plotting
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

# types lines to be plotted
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

# Brewer color sets to be avalible
kBrewerList <-
  c("Accent",
    "Dark2",
    "Paired",
    "Pastel1",
    "Pastel2",
    "Set1",
    "Set2",
    "Set3")

# color Brewer set that is active to use in plot
kListColorSet <- brewer.pal(8, kBrewerList[6])

# math options avalible
kMathOptions <- c("mean", "sum", "median", "var")

# functions ----
RgbToHex <- function(my_hex = NULL, my_rgb = NULL, tint = FALSE){
  if(!is.null(my_hex)){
    if (is.numeric(tint)) {
      my_rgb <- as.numeric(col2rgb(c(my_hex)))
      my_rgb <- paste(round(my_rgb + (255 - my_rgb) * tint),collapse = ",")
    } else {
    return(paste(col2rgb(c(my_hex)),collapse = ","))
    }
  }
  if(!is.null(my_rgb)){
    red_green_blue <- strsplit(my_rgb, ",")[[1]]
    if (length(red_green_blue) == 3 & sum(between(as.numeric(red_green_blue),0,255))== 3) {
      my_hex <- rgb(
        as.numeric(red_green_blue)[1],
        as.numeric(red_green_blue)[2],
        as.numeric(red_green_blue)[3],
        maxColorValue = 255
      )
    }
    return(my_hex)
  }
}

# basic test for valid colors
isColor <- function(x) {
  res <- try(col2rgb(x), silent = TRUE)
  return(!"try-error" %in% class(res))
}

# reads in file, tests, fills out info and returns list_data or gene list
LoadTableFile <- function(file_path, file_name, list_data, load_gene_list = FALSE, convert = FALSE) {
  if(length(file_name) == 1 && length(grep(".url",file_name))==1){
    file_path <- read_lines(file_path)
    file_name <- NULL
    for(i in file_path){
      file_name <- c(file_name, last(strsplit(i, "/")[[1]]))
    }
  }
  file_count <- length(list_data$table_file)
  for (x in seq_along(file_path)) {
    legend_nickname <-
      strsplit(as.character(file_name[x]), '.tab')[[1]][1]
    if (any(legend_nickname == names(list_data$table_file))) {
      showModal(modalDialog(
        title = "Information message",
        paste(file_name[x], "has already been loaded"), size = "s",
        easyClose = TRUE
      ))
      next()
    }
    setProgress(1, detail = "numbins")
    
      num_bins <-
        count_fields(file_path[x],
                     n_max = 1,
                     skip = 1,
                     tokenizer = tokenizer_tsv())
      
      if(num_bins == 1 & load_gene_list){
        tablefile <-
          suppressMessages(read_tsv(
            file_path,
            col_names = "gene",
            comment = "#",
            cols(gene = col_character())
          ))
      } else {
      if(num_bins == 6){
          col_names <- c("chr", "start", "end", "gene", "bin", "score")
        } else if(num_bins == 3){
          col_names <- c("gene", "bin", "score")
        } else {
          showModal(modalDialog(
            title = "Information message",
            " I dont know how to load this file, I use windowed bed files ", size = "s",
            easyClose = TRUE
          ))
          next()
        }
      setProgress(1, detail = "load file")
      tablefile <-
          suppressMessages(read_tsv(
            file_path[x],
            comment = "#",
            col_names = col_names
          )) %>%
          select(gene, bin, score) %>%
          mutate(set = legend_nickname) %>% na_if(Inf)
      num_bins <- n_distinct(tablefile$bin)
        if (file_count > 0 & num_bins != list_data$x_plot_range[2]) {
            showModal(modalDialog(
              title = "Information message",
              "Can't load file, different number of bins", size = "s",
              easyClose = TRUE
            ))
            next ()
          }
      }
      setProgress(2, detail = "process gene list")
      if(load_gene_list){
        gene_names <- semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
        if (n_distinct(gene_names$gene) == 0 & !convert) {

          showModal(modalDialog(
            title = "Information message",
            " No genes in common, might need to reformat gene name style, try pattern matching", size = "s",
            easyClose = TRUE
          ))
          return()

        } else if (n_distinct(gene_names$gene) == 0 & convert){
          setProgress(2, detail = "looking for gene name matches")
          gene_names <- distinct(tibble(gene = grep(paste(genefile$gene, collapse = "|"),
                                                   pull(distinct(list_data$gene_file[[1]]$use)), value = T)))
          if (n_distinct(gene_names$gene) == 0) {
            showModal(modalDialog(
              title = "Information message",
              " No genes found after pattern matching search", size = "s",
              easyClose = TRUE
            ))
            return()
          }
          showModal(modalDialog(
            title = "Information message",
            " Don't forget to save the gene list for future use", size = "s",
            easyClose = TRUE
          ))
        }
        legend_nickname <-
          paste(strsplit(as.character(file_name), '.txt')[[1]][1], "\nn = ", n_distinct(gene_names$gene), sep = "")
        setProgress(3, detail = "adding file to lists")
        list_data$gene_file[[legend_nickname]]$full <- distinct(tablefile, gene)
        list_data$gene_file[[legend_nickname]]$use <- gene_names
        nn <- length(list_data$gene_file)
        list_data$gene_info[[legend_nickname]] <-
          lapply(setNames(
            names(list_data$gene_info[[1]]),
            names(list_data$gene_info[[1]])
          ),
          function(i)
            tibble(
              set = i,
              mydot = kDotOptions[1],
              myline = kLineOptions[1],
              mycol = RgbToHex(my_hex = list_data$gene_info[[1]][[i]]$mycol, tint = nn*.15),
              onoff = 0,
              rnorm = "1"
            ))
        list_data$STATE[2] <- legend_nickname
       } else {
      
        zero_genes <-
          group_by(tablefile, gene) %>% summarise(test = sum(score, na.rm = T)) %>% filter(test != 0)
        
        tablefile <- semi_join(tablefile, zero_genes, by = "gene")
   
    if (file_count > 0) {
      gene_names <- semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
      if (n_distinct(gene_names$gene) == 0) {
        showModal(modalDialog(
          title = "Information message",
          " No genes in common ", size = "s",
          easyClose = TRUE
        ))
        break()
      }
    } else {
      gene_names <- distinct(tablefile, gene)
      list_data$x_plot_range <- c(1, num_bins)

    }
      
    my_name <- paste("common\nn =", n_distinct(gene_names$gene))
    list_data$STATE[2] <- my_name
    if (file_count > 0) {
    names(list_data$gene_file)[1] <- my_name
    names(list_data$gene_info)[1] <- my_name
    }
    color_safe <-
      (length(list_data$table_file) + 1) %% length(kListColorSet)
    if (color_safe == 0) {
      color_safe <- 1
    }
    color_select <- kListColorSet[color_safe]
    setProgress(3, detail = "build")
  
    list_data$table_file[[legend_nickname]] <- tablefile
    list_data$gene_file[[my_name]]$use <- gene_names
    list_data$gene_info[[my_name]][[legend_nickname]] <-
      # don't change the order of postions
      tibble(
        set = legend_nickname,
        mydot = kDotOptions[1],
        myline = kLineOptions[1],
        mycol = color_select,
        onoff = legend_nickname,
        rnorm = "1"
      )
    setProgress(4, detail = "adding to gene list")
    
    # generate info for new file for loaded gene list(s)
    sapply(names(list_data$gene_file), function(g) {
      if (length(list_data$gene_file[[g]]$full) > 0) {
        enesg <- semi_join(list_data$gene_file[[g]]$full, gene_names, by = "gene")
        if(n_distinct(enesg$gene) < 1){
          showModal(modalDialog(
            title = "Information message",
            " No genes in common, need to remove gene file", size = "s",
            easyClose = TRUE
          ))
        }
        my_name_g <- paste0(strsplit(g, "\nn =")[[1]][1], "\nn = ", n_distinct(enesg$gene))
        nn <- which(names(list_data$gene_info) == g)
        names(list_data$gene_file)[nn] <<- my_name_g
        names(list_data$gene_info)[nn] <<- my_name_g
        list_data$gene_file[[my_name_g]]$use <<- enesg
       
        list_data$gene_info[[my_name_g]][[legend_nickname]] <<-
          tibble(
            set = legend_nickname,
            mydot = kDotOptions[1],
            myline = kLineOptions[1],
            mycol = RgbToHex(my_hex = color_select, tint = nn*.15),
            onoff = 0,
            rnorm = "1"
          )
      }
    })
    file_count <- 1
  }
  setProgress(5, detail = "done")
  }
  list_data
  
}

# read in and match up names and change colors
LoadColorFile <- function(file_path, list_data) {
   
  num_bins <-
    count_fields(file_path,
                 n_max = 1,
                 skip = 1,
                 tokenizer = tokenizer_tsv())
  if(num_bins == 2){
    color_file <-
      suppressMessages(read_tsv(
        file_path,
        col_names = F,
        col_types = "cc")
      )
  } else if(num_bins == 1){
    num_bins <-
      count_fields(file_path,
                   n_max = 1,
                   skip = 1,
                   tokenizer = tokenizer_delim(" "))
    if(num_bins == 2){
      color_file <-
        suppressMessages(read_delim(
          delim = " ",
          file_path,
          col_names = F,
          col_types = "cc")
        )
    } else {
      print("can't load file")
      return(list_data)
    }
    # match name test color and update color
    for(i in seq_along(color_file$X1)){
      nickname <-
        strsplit(as.character(color_file$X1[i]), '.tab')[[1]][1]
      num <- grep(nickname, names(list_data$table_file), ignore.case = TRUE)

      if(length(num) > 0){
        if (suppressWarnings(!is.na(as.numeric(substr(
          color_file$X2[i], 1, 1
        )))) == TRUE) {
          red_green_blue <- strsplit(color_file$X2[i], ",")
          if (length(red_green_blue[[1]]) == 3) {
            color_file$X2[i] <- rgb(
              as.numeric(red_green_blue[[1]])[1],
              as.numeric(red_green_blue[[1]])[2],
              as.numeric(red_green_blue[[1]])[3],
              maxColorValue = 255
            )
          } else {
            color_file$X2[i] <- "black"
          }
        }
        
        if (!isColor(color_file$X2[i])) {
          color_file$X2[i] <- "black"
        }
        
       lapply(seq_along(list_data$gene_info), function(j) {
            list_data$gene_info[[j]][[nickname]]['mycol'] <<- color_file$X2[i]
          })
    }
 
  } 
  }
  list_data
}

# records check box on/off
CheckBoxOnOff <- function(check_box, list_data) {
  for(j in names(list_data)){
    for (i in names(list_data[[j]])) {
      # make function
      if (!i %in% check_box[[j]]) {
        list_data[[j]][[i]]["onoff"] <- 0
      } else {
        list_data[[j]][[i]]["onoff"] <- i
      }
    }
  }
  list_data
}

# sorts active gene list contain top % signal based on selected bins and file
SortTop <- function(list_data, list_name, file_names, start_bin, end_bin, num, topbottom) {
  if (is.null(file_names)) {
    return (data.frame(gene=NA, score=0, myper=0))
  }
  lc <- 0
  outlist <- NULL
  nick_name2 <- strsplit(list_name, "\nn =")[[1]][1]
  lapply(file_names, function(j) {
    enesg <-
      list_data$gene_file[[grep(nick_name2, names(list_data$gene_file),value = T)]]$use
    df <-
      semi_join(list_data$table_file[[j]], enesg, by = 'gene')
      apply_bins <- group_by(df, gene) %>%
        filter(bin %in% start_bin:end_bin) %>%
        summarise(mysums = sum(score, na.rm = TRUE)) %>%
        mutate(myper = percent_rank(mysums)) %>%
        ungroup()
    
    gene_count <- nrow(apply_bins)
    
    if (topbottom == "Top%") {
      num2 <- c(1, ceiling(gene_count * (num / 100)))
    } else if (topbottom == "Quick%") {
      num2 <-
        c(1, count(apply_bins, myper >= mean(myper, na.rm = TRUE))[[2]][2])
    } else {
      num2 <-
        c(ceiling((gene_count + 1) - (gene_count * (num / 100))), gene_count)
    }
    outlist2 <- arrange(apply_bins, desc(mysums)) %>%
      slice(num2[1]:num2[2])
    if (lc > 0) {
      outlist <<- inner_join(outlist, outlist2, by = 'gene')
      names(outlist)[2] <<- "mysums"
    } else {
      outlist <<- outlist2
    }
    lc <<- lc + 1
  })
  old_name <- grep("Sort", names(list_data$gene_file),value = T)
  if(length(old_name) > 0){
    list_data$gene_file[[old_name]] <- NULL
    list_data$gene_info[[old_name]] <- NULL
  }
  
  nick_name <- paste("Sort\nn =", n_distinct(outlist$gene) )#, "\n", nick_name2)
  list_data$gene_file[[nick_name]]$full <- outlist
  list_data$gene_file[[nick_name]]$use <- select(outlist, gene)
  list_data$gene_info[[nick_name]] <-
    lapply(setNames(
      names(list_data$gene_info[[1]]),
      names(list_data$gene_info[[1]])
    ),
    function(i)
      tibble(
        set = i,
        mydot = kDotOptions[1],
        myline = kLineOptions[1],
        mycol = RgbToHex(my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nick_name)]][[i]]$mycol, 
                         tint = length(list_data$gene_file)*.15),
        onoff = 0,
        rnorm = "1"
      ))
  
  list_data$STATE[c(2,4)] <- c(nick_name, 3)
  list_data  
}

# Applys math to on data in each gene list
ApplyMath <-
  function(list_data, 
           use_math, 
           gene_relative_frequency,
           checkboxrf,
           normbin,
           sel_list = NULL) {
    print("apply math fun")
    table_file = list_data$table_file
    gene_file = list_data$gene_file
    gene_info = list_data$gene_info
    list_data_frame <- NULL
    list_long_data_frame <- NULL
    if(sum(unlist(sapply(names(gene_info), function(i) sapply(gene_info[[i]], "[[",5) != 0))) == 0){
      return(NULL)
    }
      for (i in names(gene_file)) {
        # checks to see if at least one file in list is acitve
        if (sum(sapply(gene_info[[i]], "[[", 5) != 0) == 0) {
          next ()
        } else {
          if (!is.null(sel_list)) { 
            enesg <- c(sel_list, gene_file[[i]]$use)
            enesg <- data_frame(gene = enesg[duplicated(enesg)])
            if (length(enesg[[1]]) == 0) {
              break()
            }
          } else {
            enesg <- gene_file[[i]]$use
          }
          truefalse <-
            c(sapply(
              gene_info[[i]], "[[", 5
            ) != 0)
          list_data_frame[[i]] <-
            bind_rows(table_file[truefalse]) %>%
            semi_join(., enesg, by = "gene") 
        }
      
      if (is.null(names(list_data_frame))) {
        print("nothing to plot")
        return(NULL)
      }
    # applys math to pared down data file
    if (gene_relative_frequency) {
      list_long_data_frame[[i]] <- bind_rows(list_data_frame) %>%
        group_by(set, gene) %>%
        mutate(score = score / sum(score, na.rm = TRUE)) %>%
        ungroup() %>%
        group_by(set, bin) %>%
        summarise(value = get(use_math)(score, na.rm = T)) %>%
        ungroup() %>%
        mutate(., set = paste(gsub("(.{17})", "\\1\n", i), gsub("(.{17})", "\\1\n", set), sep = '\n'))
      
    } else {
      list_long_data_frame[[i]] <- bind_rows(list_data_frame) %>%
        group_by(set, bin) %>%
        summarise(value = get(use_math)(score, na.rm = T)) %>%
        ungroup() %>%
        mutate(., set = paste(gsub("(.{17})", "\\1\n", i), gsub("(.{17})", "\\1\n", set), sep = '\n'))
    }
    if(normbin > 0) {
      list_long_data_frame[[i]] <-
        group_by(list_long_data_frame[[i]], set) %>%
        mutate(value = value / nth(value, normbin)) %>%
        ungroup()
    } else if(checkboxrf){
      list_long_data_frame[[i]] <- group_by(list_long_data_frame[[i]], set) %>%
        mutate(value = value / sum(value)) %>%
        ungroup()
    }
    
        list_data_frame <- NULL
      }
    return(bind_rows(list_long_data_frame))
  }

# gather relavent plot option data
MakePlotOptionFrame <- function(list_data){
  print("plot options fun")
  gene_info <- list_data$gene_info
  list_data_frame <- NULL
  for (i in names(gene_info)) {
    # checks to see if at least one file in list is acitve
    if (sum(sapply(gene_info[[i]], "[[", 5) != 0) == 0) {
      next ()
    } else {
      truefalse <-
        c(sapply(
          gene_info[[i]], "[[", 5
        ) != 0)
      my_lines <- match(sapply(gene_info[[i]][truefalse], "[[", 3), kLineOptions)
      my_dots <- match(sapply(gene_info[[i]][truefalse], "[[", 2), kDotOptions)
     
      list_data_frame[[i]] <-
        bind_rows(gene_info[[i]][truefalse]) %>%
        mutate(myline = if_else(my_lines > 6, 0, as.double(my_lines)),
               mydot = if_else(my_dots == 1, 0, my_dots + 13),
               mysize = if_else(my_dots == 1, 0.01, 4.5),
               set = paste(gsub("(.{17})", "\\1\n", i), gsub("(.{17})", "\\1\n", set), sep = '\n'))
    }
  }
  if (!is.null(names(list_data_frame))) {
    return(bind_rows(list_data_frame))
  } else {
    print("no options")
    return(NULL)
  }
}

# Sets y lable fix
YAxisLable <- function(use_math = "mean",  relative_frequency = F, gene_relative_frequency = F, norm_bin = 0, smoothed = F){
  use_y_label <- paste(use_math, "of bin counts")
  if (gene_relative_frequency) {
    use_y_label <- paste("RF per gene :", use_y_label)
  } else if (relative_frequency) {
    use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                         "bins : RF")
  }
  if (norm_bin > 0) {
    if (gene_relative_frequency) {
      use_y_label <- paste(use_y_label, " : Norm bin ", norm_bin)
    } else {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : Normalize to bin ",
                           norm_bin)
    }
  }
  
  if(smoothed){
    use_y_label <- paste0("smoothed(", use_y_label, ")")
  }
  use_y_label
}

# Sets plot lines and lables fix

LinesLablesList <- function(body1bin = 20, 
                            body2bin = 40, 
                            tssbin = 15,
                            tesbin = 45,
                            binbp = 100,
                            totbins = 80,
                            everybin = 5){
  print("lines and lables fun")
  if(tssbin > 0 & tesbin > 0 & (body1bin == 0 | body2bin == 0)){
    mytype <- "4'"
  } else if(tssbin == 0 & tesbin > 0){
    mytype <- "3'"
  } else if(tesbin == 0 & tssbin > 0){
    mytype <- "5'"
  } else if(tssbin > 0 & tesbin > 0 & body1bin > 0 & body2bin > 0){
    mytype <- "543"
  }else {
    mytype <- "none"
  }
  everybp <- everybin * binbp
 if(everybp > 0){
  if(mytype == "543"){
    LOC1 <- rev(seq(body1bin + 1, by = -everybin, length.out = (body1bin/everybin) + 1))
    LOC1[near(LOC1, tssbin, tol = everybin -1)] <- tssbin + .5
    LOC1[near(LOC1, body1bin, tol = everybin -1)] <- body1bin + .5
    LOC2 <- seq(body2bin, by = everybin, length.out = (body2bin/everybin) + 1)
    LOC2[near(LOC2, tesbin, tol = everybin -1)] <- tesbin + .5
    LOC2[near(LOC2, body2bin, tol = everybin -1)] <- body2bin + .5
    LOCname1 <- rev(seq((body1bin-tssbin)*binbp, by = -everybp, length.out = (body1bin/everybin) + 1))
    LOCname1[near(LOC1, tssbin, tol = everybin -1)] <- "TSS"
    LOCname2 <- abs(seq(-(tesbin-body2bin)*binbp, by = everybp, length.out = (body2bin/everybin) + 1))
    LOCname2[near(LOC2, tesbin, tol = everybin -1)] <- "TES"
    use_plot_breaks <- c(LOC1, LOC2)
    use_plot_breaks_labels <- c(LOCname1, LOCname2)
    use_virtical_line <- c(tssbin, tesbin, body1bin, body2bin) + .5
  } else if (mytype == "5'"){
    use_plot_breaks <- seq(1, by = everybin, length.out = (totbins/everybin))
    use_plot_breaks[near(use_plot_breaks, tssbin, tol = everybin -1)] <- tssbin + .5
    use_plot_breaks_labels <- seq(-everybp, by = everybp, length.out = length(use_plot_breaks))
    use_plot_breaks_labels[near(use_plot_breaks, tssbin, tol = everybin -1)] <- "TSS"
    use_virtical_line <- c(tssbin, NA, NA, NA) + .5
  } else if(mytype == "3'"){
    use_plot_breaks <- seq(1, by = everybin, length.out = (totbins/everybin))
    use_plot_breaks[near(use_plot_breaks, tesbin, tol = everybin -1)] <- tesbin + .5
    use_plot_breaks_labels <- abs(seq(-everybp, by = everybp, length.out = length(use_plot_breaks)))
    #use_plot_breaks_labels <- abs(seq(-(tesbin-totbins)*binbp, by = everybp, length.out = length(use_plot_breaks)))
    use_plot_breaks_labels[near(use_plot_breaks, tesbin, tol = everybin -1)] <- "TES"
    use_virtical_line <- c(NA, tesbin, NA, NA) + .5
  } else if(mytype == "none"){
    use_plot_breaks <- seq(1, by = everybin, length.out = (totbins/everybin))
    use_plot_breaks_labels <- c(1,rev(seq(totbins, by = -everybin, length.out = (totbins/everybin))))
    use_virtical_line <- c(NA, NA, NA, NA) + .5
  } else{
    LOC1 <- rev(seq(tssbin + 1, by = -everybin, length.out = (tssbin/everybin) + 1))
    LOC1[near(LOC1, tssbin, tol = everybin -1)] <- tssbin + .5
    LOC2 <- seq(tesbin, by = everybin, length.out = (tesbin/everybin) + 1)
    LOC2[near(LOC2, tesbin, tol = everybin -1)] <- tesbin + .5
    LOCname1 <- rev(seq(binbp, by = -everybp, length.out = (tssbin/everybin) + 1))
    LOCname1[near(LOC1, tssbin, tol = everybin -1)] <- "TSS"
    LOCname2 <- abs(seq(-binbp, by = everybp, length.out = (tesbin/everybin) + 1))
    LOCname2[near(LOC2, tesbin, tol = everybin -1)] <- "TES"
    use_plot_breaks <- c(LOC1, LOC2)
    use_plot_breaks_labels <- c(LOCname1, LOCname2)
    use_virtical_line <- c(tssbin, tesbin, NA, NA) + .5
    }
 } else {
   if(mytype == "543"){
     use_plot_breaks <- c(tssbin, tesbin, body1bin, body2bin) + .5
     use_plot_breaks_labels <- c("TSS", "TES", "5|4", "4|3")
     use_virtical_line <- c(tssbin, tesbin, body1bin, body2bin) + .5
   } else if (mytype == "5'"){
     use_plot_breaks <- tssbin + .5
     use_plot_breaks_labels <- "TSS"
     use_virtical_line <- c(tssbin, NA, NA, NA) + .5
   } else if(mytype == "3'"){
     use_plot_breaks <- tesbin + .5
     use_plot_breaks_labels <- "TES"
     use_virtical_line <- c(NA, tesbin, NA, NA) + .5
   } else if(mytype == "none"){
     use_plot_breaks <- .5
     use_plot_breaks_labels <- "none"
     use_virtical_line <- c(NA, NA, NA, NA) + .5
   }else{
     use_plot_breaks <- c(tssbin, tesbin) + .5
     use_plot_breaks_labels <- c("TSS", "TES")
     use_virtical_line <- c(tssbin, tesbin, NA, NA) + .5
   }
  }

  use_virtical_line_color <- c("green", "red", "black", "black")
  
  use_virtical_line_type <- c(1,1,1,1)
  use_plot_breaks <- na_if(use_plot_breaks, 0.5)
  use_virtical_line <- na_if(use_virtical_line, 0.5)
  use_plot_breaks_labels <- use_plot_breaks_labels[!is.na(use_plot_breaks)]
  use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
  
  
  use_virtical_line_type <-
    use_virtical_line_type[!is.na(use_virtical_line)]
  use_virtical_line_color <-
    use_virtical_line_color[!is.na(use_virtical_line)]
  use_virtical_line <- use_virtical_line[!is.na(use_virtical_line)]
  
  list(myline = virtical_line_data_frame <- data.frame(
    use_virtical_line,
    use_virtical_line_type,
    use_virtical_line_color,
    stringsAsFactors = FALSE
  ),
  mybrakes = use_plot_breaks,
  mylables = use_plot_breaks_labels)

}

# lines and labels preset helper
LinesLablesPreSet <- function(mytype){
  # 5|4, 4|3, tss, tes, bp/bin, every bin
  if(mytype == "543 bins 20,20,40"){
    tt <- c(20,40,15,45,100,5)
  } else if (mytype == "543 bins 10,10,10"){ 
    tt <- c(10,20,5,25,100,5)
  } else if(mytype == "5' 1k 1k 80bins"){
    tt <- c(0,0,40,0,25,20)
  } else if(mytype == "5' .25k 10k 205bins"){
    tt <- c(0,0,5,0,50,6)
  } else if(mytype == "3'"){
    tt <- c(0,0,0,40,25,20)
  } else{
    tt <- c(0,0,15,45,100,5)
  }
  tt
}

# help get min and max from apply math data set
MyXSetValues <- function(apply_math, xBinRange) {
  
  tt <- group_by(apply_math, set) %>%
    filter(bin %in% xBinRange[1]:xBinRange[2]) %>%
    ungroup() %>%
    summarise(min(value, na.rm = T), max(value, na.rm = T)) %>% 
    unlist(.,use.names=FALSE) 
  tt <- round(c(tt, tt[1]-tt[1]*.1, tt[2]+tt[2]*.1),4) 
}

# main ggplot function
GGplotLineDot <-
function(list_long_data_frame, xBinRange, plot_options, yBinRange, line_list, use_smooth, use_y_label) {
    print("ggplot")
    
    use_col <- plot_options$mycol
    use_dot <- plot_options$mydot
    use_line <- plot_options$myline
    names(use_col) <- plot_options$set
    names(use_dot) <- plot_options$set
    names(use_line) <- plot_options$set
    legend_space <- max(1, (lengths(
      strsplit(plot_options$set, "\n")
    )))
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
    if (use_smooth) { 
      gp <- gp +
        geom_smooth(se = FALSE,
                    size = 2.5,
                    span = .2)
    } else{
      gp <- gp +
        geom_line(size = 2.5)
    }
    gp <- gp +
      geom_point(stroke = .001) +
      scale_size_manual(values = plot_options$mysize) +
      scale_color_manual(values = use_col) + 
      scale_shape_manual(values = use_dot) +
      scale_linetype_manual(values = use_line) +
      # xlab(use_x_label) + 
      ylab(use_y_label) +  # Set axis labels
      scale_x_continuous(breaks = line_list$mybrakes,
                         labels = line_list$mylables) +

      geom_vline(
        data = line_list$myline,
        aes(xintercept = use_virtical_line),
        size = 2,
        linetype = line_list$myline$use_virtical_line_type,
        color = line_list$myline$use_virtical_line_color
      ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      theme(axis.title.y = element_text(size =  15, margin = margin(2,10,2,2))) +
      theme(axis.title.x =  element_blank()) + # element_text(size =  10, vjust = .5)) +
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
        legend.text = element_text(size = 12)
      )  +
      coord_cartesian(xlim = xBinRange, ylim = unlist(yBinRange))
    suppressMessages(print(gp))
    return(suppressMessages(gp))
  }

