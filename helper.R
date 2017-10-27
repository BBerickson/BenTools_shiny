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
    "DT",
    "shinydashboard",
    "shinyWidgets",
    "shinycssloaders",
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
  clust = list(),
  # Cluster holder
  x_plot_range = c(0, 0),
  STATE = c(0, "common", 0, 0, 0) # flow control,
  # [1] 1 = at least on file has been loadded and lets reactives fill in info
  #     2 = lets reactive change tab toggle plot button
  # [2] name of most recent loaded gene list, for setting options select,
  # [3] name of most recent deleted gene list, for resetting tool datatables
  # [4] 1 = first time switching tab auto ploting
  #     2 = on/off reactive can deactivate plot options until plot button is pressed
  #     3 = picker(s) have been remade but no reploting is needed so don't show plot button
  # [5] 1 = plot norm applymath has been run dont retrigger
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
RgbToHex <- function(my_hex = NULL,
                     my_rgb = NULL,
                     tint = FALSE) {
  if (!is.null(my_hex)) {
    if (is.numeric(tint)) {
      my_rgb <- as.numeric(col2rgb(c(my_hex)))
      my_rgb <-
        paste(round(my_rgb + (255 - my_rgb) * tint), collapse = ",")
    } else {
      return(paste(col2rgb(c(my_hex)), collapse = ","))
    }
  }
  if (!is.null(my_rgb)) {
    red_green_blue <- strsplit(my_rgb, ",")[[1]]
    if (length(red_green_blue) == 3 &
        sum(between(as.numeric(red_green_blue), 0, 255)) == 3) {
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
LoadTableFile <-
  function(file_path,
           file_name,
           list_data,
           load_gene_list = FALSE,
           convert = FALSE) {
    if (length(file_name) == 1 && length(grep(".url", file_name)) == 1) {
      file_path <- read_lines(file_path)
      file_name <- NULL
      if(load_gene_list){
        file_path <- file_path[1]
        showModal(modalDialog(
          title = "Information message",
          paste("only loading first gene list"),
          size = "s",
          easyClose = TRUE
        ))
      }
      for (i in file_path) {
        file_name <- c(file_name, last(strsplit(i, "/")[[1]]))
      }
    }
    file_count <- length(list_data$table_file)
    for (x in seq_along(file_path)) {
      legend_nickname <-
        strsplit(as.character(file_name[x]), '.tab')[[1]][1]
      if (any(legend_nickname == names(list_data$table_file) | legend_nickname == names(list_data$gene_info))) {
        showModal(modalDialog(
          title = "Information message",
          paste(file_name[x], "has already been loaded"),
          size = "s",
          easyClose = TRUE
        ))
        next()
      }
      setProgress(1, detail = "numbins")
      
      num_bins <-
        count_fields(
          file_path[x],
          n_max = 1,
          tokenizer = tokenizer_tsv()
        )
      
      if (num_bins == 1 & load_gene_list) {
        tablefile <-
          suppressMessages(read_tsv(
            file_path,
            col_names = "gene",
            comment = "#",
            cols(gene = col_character())
          ))
      } else {
        if (num_bins > 6){
          tablefile <- suppressMessages(read_tsv (file_path[x], comment = "#",
                                                  col_names = c("gene", 1:(num_bins - 1)),
                                                  skip = 1) %>% 
                                          gather(., bin, score, 2:(num_bins))) %>% 
            select(gene, bin, score) %>%
            mutate(set = legend_nickname, 
                   bin = as.numeric(bin), 
                   score = as.numeric(score)) %>%
            na_if(Inf) %>%
            replace_na(list(score = 0))
          
        } else {
          if (num_bins == 6) {
          col_names <- c("chr", "start", "end", "gene", "bin", "score")
        } else if (num_bins == 3) {
          col_names <- c("gene", "bin", "score")
        } else {
          showModal(
            modalDialog(
              title = "Information message",
              " I dont know how to load this file, I use windowed bed files ",
              size = "s",
              easyClose = TRUE
            )
          )
          next()
        }
        setProgress(1, detail = "load file")
        tablefile <-
          suppressMessages(read_tsv(file_path[x],
                                    comment = "#",
                                    col_names = col_names)) %>%
          select(gene, bin, score) %>%
          mutate(set = legend_nickname) %>% na_if(Inf) %>%
          replace_na(list(score = 0))
        }
        num_bins <- n_distinct(tablefile$bin)
        all_empty_bin <- group_by(tablefile, bin) %>% summarise(mytest = sum(score) == 0) %>% pull(mytest)
        if(any(all_empty_bin)){
          showModal(
            modalDialog(
              title = "Information message",
              "!!!! All scores for at least one bin are empty !!!!",
              size = "s",
              easyClose = TRUE
            )
          )
        }
        if (file_count > 0 & num_bins != list_data$x_plot_range[2]) {
          showModal(
            modalDialog(
              title = "Information message",
              "Can't load file, different number of bins",
              size = "s",
              easyClose = TRUE
            )
          )
          next ()
        }
      }
      setProgress(2, detail = "process gene list")
      if (load_gene_list) {
        gene_names <-
          semi_join(tablefile, list_data$gene_file[[1]]$use, by = "gene") %>% distinct(gene)
        if (n_distinct(gene_names$gene) == 0 & !convert) {
          showModal(
            modalDialog(
              title = "Information message",
              " No genes in common, might need to reformat gene name style, try pattern matching",
              size = "s",
              easyClose = TRUE
            )
          )
          return()
          
        } else if (n_distinct(gene_names$gene) == 0 & convert) {
          setProgress(2, detail = "looking for gene name matches")
          gene_names <-
            distinct(tibble(gene = grep(
              paste(tablefile$gene, collapse = "|"),
              pull(distinct(list_data$gene_file[[1]]$use)),
              value = T
            )))
          if (n_distinct(gene_names$gene) == 0) {
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
        legend_nickname <-
          paste(strsplit(as.character(file_name), '.txt')[[1]][1],
                "\nn = ",
                n_distinct(gene_names$gene),
                sep = "")
        setProgress(3, detail = "adding file to lists")
        list_data$gene_file[[legend_nickname]]$full <-
          distinct(tablefile, gene)
        list_data$gene_file[[legend_nickname]]$use <- gene_names
        if(length(list_data$gene_file) > 5){
          mytint <- 0
        } else{
          mytint <- length(list_data$gene_file) * 0.1
        }
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
              mycol = RgbToHex(
                my_hex = list_data$gene_info[[1]][[i]]$mycol,
                tint = mytint
              ),
              onoff = 0,
              rnorm = "1"
            ))
        list_data$STATE[2] <- legend_nickname
      } else {
        zero_genes <-
          group_by(tablefile, gene) %>% summarise(test = sum(score, na.rm = T)) %>% filter(test != 0)
        
        tablefile <- semi_join(tablefile, zero_genes, by = "gene")
        
        if (file_count > 0) {
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
        # only have plot on if a plot has not been created yet for first 4 files
        if (list_data$STATE[4] == 0 &
            length(list_data$table_file) < 5) {
          oo <- legend_nickname
        } else {
          oo <- 0
        }
        list_data$table_file[[legend_nickname]] <- tablefile
        list_data$gene_file[[my_name]]$use <- gene_names
        list_data$gene_info[[my_name]][[legend_nickname]] <-
          # don't change the order of postions
          tibble(
            set = legend_nickname,
            mydot = kDotOptions[1],
            myline = kLineOptions[1],
            mycol = color_select,
            onoff = oo,
            rnorm = "1"
          )
        setProgress(4, detail = "adding to gene list")
        
        # generate info for new file for loaded gene list(s)
        sapply(seq_along(list_data$gene_file)[-1], function(g) {
          enesg <- inner_join(list_data$gene_file[[g]]$full, list_data$gene_file[[1]]$use, by = "gene")
          if (n_distinct(enesg$gene) < 1) {
            showModal(
              modalDialog(
                title = "Information message",
                " No genes in common, need to remove gene file",
                size = "s",
                easyClose = TRUE
              )
            )
          }
          list_data$gene_file[[g]]$use <<- select(enesg, gene)
          my_name_g <- sub("([0-9]+)", n_distinct(list_data$gene_file[[g]]$full$gene), names(list_data$gene_file)[g])
          names(list_data$gene_file)[g] <<- my_name_g
          names(list_data$gene_info)[g] <<- my_name_g
          if(length(list_data$gene_file) > 5){
            mytint <- 0
          } else{
            mytint <- length(list_data$gene_file) * 0.1
          }
          list_data$gene_info[[g]][[legend_nickname]] <<-
            tibble(
              set = legend_nickname,
              mydot = kDotOptions[1],
              myline = kLineOptions[1],
              mycol = RgbToHex(my_hex = color_select, tint = mytint),
              onoff = 0,
              rnorm = "1"
            )
        })
        file_count <- 1
      }
      setProgress(5, detail = "done")
    }
    list_data
    
  }

# read in and match up names and change colors
LoadColorFile <- function(file_path, list_data, gene_list) {
  num_bins <-
    count_fields(file_path,
                 n_max = 1,
                 skip = 1,
                 tokenizer = tokenizer_tsv())
  if (num_bins == 2) {
    color_file <-
      suppressMessages(read_tsv(file_path,
                                col_names = F,
                                col_types = "cc"))
  } else if (num_bins == 1) {
    num_bins <-
      count_fields(
        file_path,
        n_max = 1,
        skip = 1,
        tokenizer = tokenizer_delim(" ")
      )
    if (num_bins == 2) {
      color_file <-
        suppressMessages(read_delim(
          delim = " ",
          file_path,
          col_names = F,
          col_types = "cc"
        ))
    } else {
      print("can't load file")
      return(list_data)
    }
    # match name test color and update color
    for (i in seq_along(color_file$X1)) {
      nickname <-
        strsplit(as.character(color_file$X1[i]), '.tab')[[1]][1]
      num <-
        grep(
          nickname,
          names(list_data$table_file),
          ignore.case = TRUE,
          value = T
        )
      
      if (length(num) > 0) {
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
        
        list_data$gene_info[[gene_list]][[num]]['mycol'] <-
          color_file$X2[i]
        
      }
      
    }
  }
  list_data
}

# records check box on/off
CheckBoxOnOff <- function(check_box, list_data) {
  for (j in names(list_data)) {
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

# make normalized file ... devide one by the other
MakeNormFile <- function(list_data, nom, dnom, gbyg, nodivzero) {
  if (nom != "" && dnom != "") {
    mynom <- list_data$table_file[[nom]]
    mydom <- list_data$table_file[[dnom]]
    myname <- "gene_by_gene"
    setProgress(1, detail = "Gathering data")
    if(!gbyg){
      myname <- "mean_norm"
      mynom <- group_by(mynom, bin, set) %>% mutate(score=mean(score, na.rm = TRUE)) %>% ungroup()
      mydom <- group_by(mydom, bin, set) %>% mutate(score=mean(score, na.rm = TRUE)) %>% ungroup()
    }
    if(nodivzero){
      myname <- paste0(myname, "-0_min/2")
    # find min value /2 to replace 0s
      new_gene_list <- inner_join(mynom, mydom, by = c("gene", "bin")) %>%
        na_if(0)
    new_min_for_na <-
      min(c(new_gene_list$score.x, new_gene_list$score.y), na.rm = TRUE) / 2
    # replace 0's with min/2
    new_gene_list <- replace_na(new_gene_list, list(score.y = new_min_for_na, score.x = new_min_for_na))
    } else {
      new_gene_list <- inner_join(mynom, mydom, by = c("gene", "bin")) 
      new_gene_list <-
        group_by(new_gene_list, gene) %>%
        summarise(test = sum(score.x, score.y)) %>% 
        filter(!is.na(test)) %>% 
        semi_join(new_gene_list, ., by = "gene")
    }
    legend_nickname <- paste0(nom, "/\n", dnom,":\n", myname)
    color_safe <-
      (length(list_data$table_file) + 1) %% length(kListColorSet)
    if (color_safe == 0) {
      color_safe <- 1
    }
    gene_names <- semi_join(list_data$gene_file[[1]]$use, new_gene_list, by = "gene")
    my_name <- paste("common\nn =", n_distinct(gene_names$gene))
    list_data$STATE[2] <- my_name
    names(list_data$gene_file)[1] <- my_name
    names(list_data$gene_info)[1] <- my_name
    color_select <- kListColorSet[color_safe]
    setProgress(2, detail = "building new data")
  list_data$table_file[[legend_nickname]] <-
    transmute(
      new_gene_list,
      gene = gene,
      bin = bin,
      set = legend_nickname,
      score = score.x / score.y
    ) %>% na_if(Inf) %>% replace_na(list(score = 0))
  list_data$gene_file[[my_name]]$use <- semi_join(list_data$gene_file[[my_name]]$use, 
                                                  list_data$table_file[[legend_nickname]], by = "gene")
  list_data$gene_info[[my_name]][[legend_nickname]] <-
    # don't change the order of postions
    tibble(
      set = legend_nickname,
      mydot = kDotOptions[1],
      myline = kLineOptions[1],
      mycol = color_select,
      onoff = 0,
      rnorm = "1"
    )

  # generate info for new file for loaded gene list(s)
  setProgress(3, detail = "building gene lists info")
  sapply(seq_along(list_data$gene_file)[-1], function(g) {
    enesg <- inner_join(list_data$gene_file[[g]]$full, list_data$gene_file[[1]]$use, by = "gene")
    if (n_distinct(enesg$gene) < 1) {
      showModal(
        modalDialog(
          title = "Information message",
          " No genes in common, need to remove gene file",
          size = "s",
          easyClose = TRUE
        )
      )
    }
    list_data$gene_file[[g]]$use <<- select(enesg, gene)
    my_name_g <- sub("([0-9]+)", n_distinct(list_data$gene_file[[g]]$full$gene), names(list_data$gene_file)[g])
    names(list_data$gene_file)[g] <<- my_name_g
    names(list_data$gene_info)[g] <<- my_name_g
    if(length(list_data$gene_file) > 5){
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
    list_data$gene_info[[g]][[legend_nickname]] <<-
      tibble(
        set = legend_nickname,
        mydot = kDotOptions[1],
        myline = kLineOptions[1],
        mycol = RgbToHex(my_hex = color_select, tint = mytint),
        onoff = 0,
        rnorm = "1"
      )
  })
  }
  setProgress(5, detail = "Done")
  list_data
}

# removes gene list
RemoveGeneList <-
  function(list_data, list_name) {
      list_data$gene_file[[list_name]] <- NULL
      list_data$gene_info[[list_name]] <- NULL
      if (list_data$STATE[4] == 0) {
        list_data$STATE[2] <- names(list_data$gene_file)[1]
      } else{
        list_data$STATE[c(2, 4)] <- c(names(list_data$gene_file)[1], 3)
      }
      list_data$STATE[3] <- strtrim(list_name, 6)
      
    list_data
  }

# removes data file
RemoveFile <- function(list_data, file_name, remove_all){
  if (length(list_data$table_file) > 1 & !remove_all) {
    # remove tool gene lists ? TODO
    sapply(names(list_data$gene_file), function(g) {
      list_data$gene_file[[g]][[file_name]] <<- NULL
      list_data$gene_info[[g]][[file_name]] <<- NULL
    })
    list_data$table_file[[file_name]] <- NULL
    list_data$gene_file[[1]]$use <- distinct(list_data$table_file[[1]], gene)
    sapply(list_data$table_file[-1], function(i) {
      list_data$gene_file[[1]]$use <<- semi_join(list_data$gene_file[[1]]$use, i, by="gene")
      })
  my_name <- paste("common\nn =", n_distinct(list_data$gene_file[[1]]$use$gene))
  names(list_data$gene_file)[1] <- my_name
  names(list_data$gene_info)[1] <- my_name
  
  sapply(seq_along(list_data$gene_file)[-1], function(g) {
      list_data$gene_file[[g]]$full <<- inner_join(list_data$gene_file[[g]]$full, list_data$gene_file[[1]]$use, by = "gene")
      list_data$gene_file[[g]]$use <<- select(list_data$gene_file[[g]]$full, gene)
      my_name_g <- sub("([0-9]+)", n_distinct(list_data$gene_file[[g]]$full$gene), names(list_data$gene_file)[g])
      names(list_data$gene_file)[g] <<- my_name_g
      names(list_data$gene_info)[g] <<- my_name_g
  })
  
  if (list_data$STATE[4] == 0) {
    list_data$STATE[2] <- names(list_data$gene_file)[1]
  } else{
    list_data$STATE[c(2, 4)] <- c(names(list_data$gene_file)[1], 2)
  }
  } else {
    list_data <- list(
      table_file = list(),
      gene_file = list(),
      gene_info = list(),
      clust = list(),
      x_plot_range = c(0, 0),
      STATE = c(0, "common", 0, 0) 
      )
  }
  list_data
}

# inclusive, exclusive and intersected gene lists
IntersectGeneLists <- function(list_data, list_name){
  if(is.null(list_name)){
    return(NULL)
  } 
  setProgress(1, detail = paste("building list"))
  outlist <- NULL
  lapply(list_name, function(j){
    outlist[[j]] <<- list_data$gene_file[[j]]$use
    })
  outlist <- bind_rows(outlist)
  for(rr in grep("Gene_List_", names(LIST_DATA$gene_file), value = T)){
    if (length(rr) > 0) {
      list_data$gene_file[[rr]] <- NULL
      list_data$gene_info[[rr]] <- NULL
    }
  }
  
  nick_name <- NULL
  setProgress(2, detail = paste("building inclusive list"))
  inclusive <- distinct(outlist)
  if(n_distinct(inclusive$gene) > 0){
    nick_name1 <-
      paste("Gene_List_inclusive\nn =", n_distinct(inclusive$gene))
    nick_name <- c(nick_name, nick_name1)
    list_data$gene_file[[nick_name1]]$full <- inclusive
    list_data$gene_file[[nick_name1]]$use <- select(inclusive, gene)
    list_data$gene_file[[nick_name1]]$info <-
      paste(
        "Gene_List_inclusive",
        "from",
        paste(list_name, collapse = " and "),
        Sys.Date()
      )
  }
  setProgress(3, detail = paste("building intersect list"))
  intersect <- filter(outlist, duplicated(gene))
  if(n_distinct(intersect$gene) > 0){
    nick_name1 <-
      paste("Gene_List_intersect\nn =", n_distinct(intersect$gene))
    nick_name <- c(nick_name, nick_name1)
    list_data$gene_file[[nick_name1]]$full <- intersect
    list_data$gene_file[[nick_name1]]$use <- select(intersect, gene)
    list_data$gene_file[[nick_name1]]$info <-
      paste(
        "Gene_List_intersect",
        "from",
        paste(list_name, collapse = " and "),
        Sys.Date()
      )

  setProgress(4, detail = paste("building exclusive list"))
  exclusive <- anti_join(inclusive, intersect, by="gene")
  if(n_distinct(exclusive$gene) > 0){
    nick_name1 <-
      paste("Gene_List_exclusive\nn =", n_distinct(exclusive$gene))
    nick_name <- c(nick_name, nick_name1)
    list_data$gene_file[[nick_name1]]$full <- exclusive
    list_data$gene_file[[nick_name1]]$use <- select(exclusive, gene)
    list_data$gene_file[[nick_name1]]$info <-
      paste(
        "Gene_List_exclusive",
        "from",
        paste(list_name, collapse = " and "),
        Sys.Date()
      )
  }
  }
  setProgress(5, detail = "finishing up")
  if(length(list_data$gene_file) > 5){
    mytint <- 0
  } else{
    mytint <- length(list_data$gene_file) * 0.1
  }
  for (nn in nick_name) {
    list_data$gene_info[[nn]] <-
      lapply(setNames(
        names(list_data$gene_info[[1]]),
        names(list_data$gene_info[[1]])
      ),
      function(i)
        tibble(
          set = i,
          mydot = kDotOptions[1],
          myline = kLineOptions[1],
          mycol = RgbToHex(
            my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nn)]][[i]]$mycol,
            tint = mytint),
          onoff = 0,
          rnorm = "1"
        ))
  }
  if (list_data$STATE[4] != 0) {
    list_data$STATE[4] <- 3
  }
  list_data
}

# sorts active gene list contain top % signal based on selected bins and file
SortTop <-
  function(list_data,
           list_name,
           file_names,
           start_bin,
           end_bin,
           num,
           topbottom) {
    if (is.null(file_names)) {
      return(NULL)
    }
    lc <- 0
    outlist <- NULL
    lapply(file_names, function(j) {
      setProgress(lc + 1, detail = paste("sorting", j))
      enesg <-
        list_data$gene_file[[list_name]]$use
      apply_bins <-
        semi_join(list_data$table_file[[j]], enesg, by = 'gene')
      apply_bins <- group_by(apply_bins, gene) %>%
        filter(bin %in% start_bin:end_bin) %>%
        summarise(mysums = sum(score, na.rm = TRUE)) %>%
        mutate(myper = as.numeric(strtrim(cume_dist(mysums), 5)))
      
      gene_count <- nrow(apply_bins)
      
      if (topbottom == "Top%") {
        num2 <- c(1, ceiling(gene_count * (num / 100)))
        topbottom <- paste(topbottom, paste0(num, "%"))
      } else if (topbottom == "Quick%") {
        num2 <-
          c(1, count(apply_bins, myper >= mean(myper, na.rm = TRUE))[[2]][2])
      } else {
        num2 <-
          c(ceiling((gene_count + 1) - (gene_count * (num / 100))), gene_count)
        topbottom <- paste(topbottom, paste0(num, "%"))
      }
      nickname <- list_data$gene_info[[1]][[j]]$set
      outlist2 <- arrange(apply_bins, desc(mysums)) %>%
        mutate(!!nickname := myper) %>%
        select(gene,!!nickname) %>%
        slice(num2[1]:num2[2])
      if (lc > 0) {
        outlist <<- inner_join(outlist, outlist2, by = 'gene')
      } else {
        outlist <<- outlist2
      }
      lc <<- lc + 1
    })
    old_name <- grep("Sort", names(list_data$gene_file))
    if (length(old_name) > 0) {
      list_data$gene_file[[old_name]] <- NULL
      list_data$gene_info[[old_name]] <- NULL
    }
    setProgress(lc + 2, detail = "building list")
    nick_name <- strtrim(gsub("(.{30})", "\\1... ", paste0("Sort\nn = ", n_distinct(outlist$gene), "-", list_name)),33)
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <-
      paste("Sort",
            topbottom,
            "bins",
            start_bin,
            "to",
            end_bin,
            "from",
            list_name,
            paste(file_names, collapse = " "),
            Sys.Date())
    if(length(list_data$gene_file) > 5){
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
    list_data$gene_info[[nick_name]] <-
      lapply(setNames(names(list_data$gene_info[[1]]),
                      names(list_data$gene_info[[1]])),
             function(i)
               tibble(
                 set = i,
                 mydot = kDotOptions[1],
                 myline = kLineOptions[1],
                 mycol = RgbToHex(
                   my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nick_name)]][[i]]$mycol,
                   tint = mytint
                 ),
                 onoff = 0,
                 rnorm = "1"
               ))
    list_data$STATE[2] <- nick_name
    if (list_data$STATE[4] != 0) {
      list_data$STATE[4] <- 3
    }
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
           nodivzero) {
    if (ratio1file == "" | ratio2file == "") {
      return()
    }
    setProgress(1, detail = paste("dividing one by the other"))
    lc <- 0
    outlist <- NULL
    lapply(c(ratio1file, ratio2file), function(j) {
      enesg <-
        list_data$gene_file[[list_name]]$use
      df <- semi_join(list_data$table_file[[j]], enesg, by = 'gene') %>%
        group_by( gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T)) %>%
        na_if(0) %>% ungroup()
      if(nodivzero){
        # find min value /2 to replace 0s
        new_min_for_na <-
          min(c(df$sum1, df$sum2), na.rm = TRUE) / 2
        df <-
          group_by(df, gene) %>%
          replace_na(list(sum1 = new_min_for_na, sum2 = new_min_for_na)) %>%
          ungroup()
      } else {
        if (start2_bin == 0 | end2_bin == 0) {
          df <-
            group_by(df, gene) %>%
            summarise(test = sum(sum1)) %>% 
            filter(!is.na(test)) %>% 
            semi_join(df, ., by = "gene") %>%
            ungroup()
        } else {
          df <-
            group_by(df, gene) %>%
            summarise(test = sum(sum1, sum2)) %>% 
            filter(!is.na(test)) %>% 
            semi_join(df, ., by = "gene") %>%
            ungroup()
        }
      }
      lc <<- lc + 1
      if (start2_bin == 0 | end2_bin == 0) {
        outlist[[lc]] <<- df
      } else {
        outlist[[lc]] <<-
          transmute(df, gene = gene, sum1 = sum1 / sum2) %>%
          na_if(Inf) %>%
          replace_na(list(score = 0))
      }
      if (lc > 1) {
        outlist[[1]] <<-
          inner_join(outlist[[1]], outlist[[lc]], by = 'gene') %>%
          transmute(gene = gene, Ratio = sum1.x / sum1.y) %>%
          na_if(Inf) %>%
          replace_na(list(Ratio = 0)) %>%
          arrange(desc(Ratio)) %>% select(gene, Ratio)
        
      }
    })
    
    for(rr in grep("Ratio_", names(LIST_DATA$gene_file), value = T)){
      if (length(rr) > 0) {
        list_data$gene_file[[rr]] <- NULL
        list_data$gene_info[[rr]] <- NULL
      }
    }
    
    nick_name <- NULL
    setProgress(2, detail = paste("building list", ratio1file))
    upratio <- filter(outlist[[1]], Ratio < 1 / num)
    if(n_distinct(upratio$gene) > 0){
    nick_name1 <-
      paste("Ratio_Up_file1\nn =", n_distinct(upratio$gene))
    nick_name <- c(nick_name, nick_name1)
    list_data$gene_file[[nick_name1]]$full <- upratio
    list_data$gene_file[[nick_name1]]$use <- select(upratio, gene)
    list_data$gene_file[[nick_name1]]$info <-
      paste(
        "Ratio_Up_file1",
        ratio1file,
        "/",
        ratio2file,
        "bins",
        start1_bin,
        "to",
        end1_bin,
        "/",
        start2_bin,
        end2_bin,
        "fold change cut off",
        num,
        "0  to min/2?",
        nodivzero,
        "from",
        list_name,
        "gene list",
        Sys.Date()
      )
    }
    setProgress(3, detail = paste("building list", ratio2file))
    upratio <- filter(outlist[[1]], Ratio > num)
    if(n_distinct(upratio$gene) > 0){
    nick_name2 <-
      paste("Ratio_Up_file2\nn =", n_distinct(upratio$gene))
    nick_name <- c(nick_name, nick_name2)
    list_data$gene_file[[nick_name2]]$full <- upratio
    list_data$gene_file[[nick_name2]]$use <- select(upratio, gene)
    list_data$gene_file[[nick_name2]]$info <-
      paste(
        "Ratio_Up_file2",
        ratio2file,
        "/",
        ratio1file,
        "bins",
        start1_bin,
        "to",
        end1_bin,
        "/",
        start2_bin,
        end2_bin,
        "fold change cut off",
        num,
        "from",
        list_name,
        "gene list",
        Sys.Date()
      )
    }
    setProgress(4, detail = paste("building list: no change"))
    upratio <- filter(outlist[[1]], Ratio <= num & Ratio >= 1 / num)
    if(n_distinct(upratio$gene) > 0){
    nick_name3 <-
      paste("Ratio_No_Diff\nn =", n_distinct(upratio$gene))
    nick_name <- c(nick_name, nick_name3)
    list_data$gene_file[[nick_name3]]$full <- upratio
    list_data$gene_file[[nick_name3]]$use <- select(upratio, gene)
    list_data$gene_file[[nick_name3]]$info <-
      paste(
        "Ratio_No_Diff",
        ratio1file,
        "/",
        ratio2file,
        "bins",
        start1_bin,
        "to",
        end1_bin,
        "/",
        start2_bin,
        end2_bin,
        "fold change cut off",
        num,
        "from",
        list_name,
        "gene list",
        Sys.Date()
      )
    }
    if(length(list_data$gene_file) > 5){
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
    for (nn in nick_name) {
      list_data$gene_info[[nn]] <-
        lapply(setNames(
          names(list_data$gene_info[[1]]),
          names(list_data$gene_info[[1]])
        ),
        function(i)
          tibble(
            set = i,
            mydot = kDotOptions[1],
            myline = kLineOptions[1],
            mycol = RgbToHex(
              my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nn)]][[i]]$mycol,
              tint = mytint),
            onoff = 0,
            rnorm = "1"
          ))
    }
    setProgress(5, detail = "finishing up")
    if (list_data$STATE[4] != 0) {
      list_data$STATE[4] <- 3
    }
    list_data
  }

# Change the number of clusters
ClusterNumList <- function(list_data,list_name,
                           clusterfile, start_bin,
                           end_bin, num, myname) {

  if(is_empty(list_data$clust)){
    return(NULL)
  }
  setProgress(3, detail = "spliting into clusters")
  for(rr in grep("Cluster_", names(list_data$gene_file), value = T)){
    if (length(rr) > 0) {
      list_data$gene_file[[rr]] <- NULL
      list_data$gene_info[[rr]] <- NULL
    }
  }
  for(rr in grep("Group_", names(list_data$gene_file), value = T)){
    if (length(rr) > 0) {
      list_data$gene_file[[rr]] <- NULL
      list_data$gene_info[[rr]] <- NULL
    }
  }
  if(myname == "Cluster_"){
    gene_list <- mutate(list_data$clust$use, cm = cutree(list_data$clust$cm, num))
  } else {
    gene_list <- mutate(list_data$clust$full, cm = ntile(cm, as.numeric(num)))
  }
  for(nn in 1:num){
    outlist <- filter(gene_list, cm == nn)
    nick_name <- paste(paste0(myname, nn, "\nn ="), n_distinct(outlist$gene))
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <-
      paste(nick_name,
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
            Sys.Date())
    setProgress(4, detail = paste("finishing cluster", nn))
    if(length(list_data$gene_file) > 5){
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
    list_data$gene_info[[nick_name]] <-
      lapply(setNames(names(list_data$gene_info[[1]]),
                      names(list_data$gene_info[[1]])),
             function(i)
               tibble(
                 set = i,
                 mydot = kDotOptions[1],
                 myline = kLineOptions[1],
                 mycol = RgbToHex(
                   my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nick_name)]][[i]]$mycol,
                   tint = mytint),
                 onoff = 0,
                 rnorm = "1"
               ))
  }
  list_data$STATE[2] <- grep(paste0(myname, "1"), names(list_data$gene_file), value = T)
  if (list_data$STATE[4] != 0) {
    list_data$STATE[4] <- 3
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
    return(NULL)
  }
  setProgress(1, detail = paste("gathering data"))
  df <- semi_join(list_data$table_file[[clusterfile]], 
                  list_data$gene_file[[list_name]]$use, by = 'gene')
  setProgress(2, detail = "hierarchical clustering using ward method")
  list_data$clust <- list()
  list_data$clust$cm <- hclust.vector(as.data.frame(spread(df, bin, score))[, c((start_bin:end_bin) + 2)], method = "ward")
  list_data$clust$use <- distinct(df, gene)
  if (list_data$STATE[4] != 0) {
    list_data$STATE[4] <- 3
  }
  list_data
}

# finds 2 - 4 groups from the one active file, plotting the patterns and displaying the gene lists
FindGroups <- function(list_data,
                         list_name,
                         clusterfile,
                         start_bin,
                         end_bin,
                         num) {
  if (clusterfile == "") {
    return(NULL)
  }
  setProgress(1, detail = paste("gathering data"))
  df <- semi_join(list_data$table_file[[clusterfile]], 
                  list_data$gene_file[[list_name]]$use, by = 'gene')
  setProgress(2, detail = "finding groups")
  list_data$clust <- list()
  list_data$clust$full <- group_by(df, gene) %>% 
    filter(bin %in% start_bin:end_bin) %>%
    summarise(cm = sum(score, na.rm = TRUE))
  list_data$clust$use <- distinct(df, gene)
  if (list_data$STATE[4] != 0) {
    list_data$STATE[4] <- 3
  }
  list_data
}

# Cumulative Distribution data prep 
CumulativeDistribution <-
  function(list_data,
           list_name,
           cdffile,
           start1_bin,
           end1_bin,
           start2_bin,
           end2_bin,
           bottom_per,
           top_per,
           nodivzero) {
    if (is.null(cdffile)) {
      return()
    }
    setProgress(1, detail = paste("dividing one by the other"))
    gene_count <- n_distinct(list_data$gene_file[[list_name]]$use$gene)
    num <- c(ceiling(gene_count * bottom_per/100), ceiling(gene_count * top_per/100))
    outlist <- NULL
    lapply(cdffile, function(j) {
      df <- semi_join(list_data$table_file[[j]], list_data$gene_file[[list_name]]$use, by = 'gene') %>%
        group_by(gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T)) %>%
        na_if(0) %>% ungroup()
      if(nodivzero){
        # find min value /2 to replace 0s
        new_min_for_na <-
          min(c(df$sum1, df$sum2), na.rm = TRUE) / 2
        df <-
          group_by(df, gene) %>%
          replace_na(list(sum1 = new_min_for_na, sum2 = new_min_for_na)) %>%
          ungroup()
      } else {
        df <- group_by(df, gene) %>%
            summarise(test = sum(sum1, sum2)) %>% 
            filter(!is.na(test)) %>% 
            semi_join(df, ., by = "gene") %>%
            ungroup()
      }
      outlist[[j]] <<- transmute(df, gene = gene, value = sum1 / sum2) %>%
          na_if(Inf) %>%
          replace_na(list(value = 0)) %>%
            arrange(desc(value)) %>%
          mutate(bin = row_number(), set = j)
    })
    
    outlist <- bind_rows(outlist)
    for(rr in grep("CDF\nn", names(LIST_DATA$gene_file), value = T)){
      if (length(rr) > 0) {
        list_data$gene_file[[rr]] <- NULL
        list_data$gene_info[[rr]] <- NULL
      }
    }
    
    setProgress(2, detail = paste("building list"))
    gene_list <- group_by(outlist, gene) %>% filter(all(between(bin, num[1], num[2]))) %>% 
      distinct(gene) %>%
      ungroup()
    if(n_distinct(gene_list$gene) > 0){
      nick_name1 <-
        paste("CDF\nn =", n_distinct(gene_list$gene))
      list_data$gene_file[[nick_name1]]$full <- outlist 
      list_data$gene_file[[nick_name1]]$use <- gene_list
      list_data$gene_file[[nick_name1]]$info <-
        paste(
          "CDF",
          "bins",
          start1_bin,
          "to",
          end1_bin,
          "/",
          start2_bin,
          "to",
          end2_bin,
          "% filter",
          bottom_per,
          "to",
          top_per,
          "0  to min/2?",
          nodivzero,
          "from",
          list_name,
          "gene list",
          paste(cdffile, collapse = " "),
          Sys.Date()
        )
    }
    if (sum(start1_bin, end1_bin) > sum(start2_bin, end2_bin)) {
      use_header <- "Log2 EI Cumulative plot"
    } else {
    use_header <- "Log2 PI Cumulative plot"
    }
    if(length(list_data$gene_file) > 5){
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
      list_data$gene_info[[nick_name1]] <-
        lapply(setNames(
          names(list_data$gene_info[[1]]),
          names(list_data$gene_info[[1]])
        ),
        function(i)
          tibble(
            set = i,
            mydot = kDotOptions[1],
            myline = kLineOptions[1],
            mycol = RgbToHex(
              my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nick_name1)]][[i]]$mycol,
              tint = mytint),
            onoff = 0,
            rnorm = "1",
            myheader = use_header
          ))
    setProgress(5, detail = "finishing up")
    list_data$STATE[2] <- nick_name1
    if (list_data$STATE[4] != 0) {
      list_data$STATE[4] <- 3
    }
    list_data
  }

# Applys math to on data in each gene list
ApplyMath <-
  function(list_data,
           use_math,
           relative_frequency,
           normbin,
           sel_list = NULL) {
    print("apply math fun")
    table_file = list_data$table_file
    gene_file = list_data$gene_file
    gene_info = list_data$gene_info
    list_data_frame <- NULL
    list_long_data_frame <- NULL
    if (sum(unlist(sapply(names(gene_info), function(i)
      sapply(gene_info[[i]], "[[", 5) != 0))) == 0) {
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
          c(sapply(gene_info[[i]], "[[", 5) != 0)
        list_data_frame[[i]] <-
          bind_rows(table_file[truefalse]) %>%
          semi_join(., enesg, by = "gene")
      }
      
      if (is.null(names(list_data_frame))) {
        print("nothing to plot")
        return(NULL)
      }
      # applys math to pared down data file
      if (relative_frequency == "rel gene frequency") {
        list_long_data_frame[[i]] <- bind_rows(list_data_frame) %>%
          group_by(set, gene) %>%
          mutate(score = score / sum(score, na.rm = TRUE)) %>%
          ungroup() %>%
          group_by(set, bin) %>%
          summarise(value = get(use_math)(score, na.rm = T)) %>%
          ungroup() %>%
          mutate(., set = paste(
            gsub("(.{17})", "\\1\n", i),
            gsub("(.{17})", "\\1\n", set),
            sep = '\n'
          ))
        
      } else {
        list_long_data_frame[[i]] <- bind_rows(list_data_frame) %>%
          group_by(set, bin) %>%
          summarise(value = get(use_math)(score, na.rm = T)) %>%
          ungroup() %>%
          mutate(., set = paste(
            gsub("(.{17})", "\\1\n", i),
            gsub("(.{17})", "\\1\n", set),
            sep = '\n'
          ))
      }
      if (normbin > 0) {
        list_long_data_frame[[i]] <-
          group_by(list_long_data_frame[[i]], set) %>%
          mutate(value = value / nth(value, normbin)) %>%
          ungroup()
      } else if (relative_frequency == "relative frequency") {
        list_long_data_frame[[i]] <-
          group_by(list_long_data_frame[[i]], set) %>%
          mutate(value = value / sum(value)) %>%
          ungroup()
      }
      
      list_data_frame <- NULL
    }
    return(bind_rows(list_long_data_frame))
  }

# gather relavent plot option data
MakePlotOptionFrame <- function(list_data) {
  print("plot options fun")
  gene_info <- list_data$gene_info
  list_data_frame <- NULL
  for (i in names(gene_info)) {
    # checks to see if at least one file in list is acitve
    if (sum(sapply(gene_info[[i]], "[[", 5) != 0) == 0) {
      next ()
    } else {
      truefalse <-
        c(sapply(gene_info[[i]], "[[", 5) != 0)
      my_lines <-
        match(sapply(gene_info[[i]][truefalse], "[[", 3), kLineOptions)
      my_dots <-
        match(sapply(gene_info[[i]][truefalse], "[[", 2), kDotOptions)
      
      list_data_frame[[i]] <-
        bind_rows(gene_info[[i]][truefalse]) %>%
        mutate(
          myline = if_else(my_lines > 6, 0, as.double(my_lines)),
          mydot = if_else(my_dots == 1, 0, my_dots + 13),
          mysize = if_else(my_dots == 1, 0.01, 4.5),
          set = paste(
            gsub("(.{17})", "\\1\n", i),
            gsub("(.{17})", "\\1\n", set),
            sep = '\n'
          )
        )
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
YAxisLable <-
  function(use_math = "mean",
           relative_frequency = "none",
           norm_bin = 0,
           smoothed = F) {
    use_y_label <- paste(use_math, "of bin counts")
    if (relative_frequency == "rel gene frequency") {
      use_y_label <- paste("RF per gene :", use_y_label)
    } else if (relative_frequency == "relative frequency") {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : RF")
    }
    if (norm_bin > 0) {
      if (relative_frequency== "rel gene frequency") {
        use_y_label <- paste(use_y_label, " : Norm bin ", norm_bin)
      } else {
        use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                             "bins : Normalize to bin ",
                             norm_bin)
      }
    }
    
    if (smoothed) {
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
                            everybin = 5) {
  print("lines and lables fun")
  if (tssbin > 0 & tesbin > 0 & (body1bin == 0 | body2bin == 0)) {
    mytype <- "4'"
  } else if (tssbin == 0 & tesbin > 0) {
    mytype <- "3'"
  } else if (tesbin == 0 & tssbin > 0) {
    mytype <- "5'"
  } else if (tssbin > 0 & tesbin > 0 & body1bin > 0 & body2bin > 0) {
    mytype <- "543"
  } else {
    mytype <- "none"
  }
  everybp <- everybin * binbp
  if (everybp > 0) {
    if (mytype == "543") {
      LOC1 <-
        rev(seq(
          body1bin + 1,
          by = -everybin,
          length.out = (body1bin / everybin) + 1
        ))
      LOC1[near(LOC1, tssbin, tol = everybin - 1)] <- tssbin + .5
      LOC1[near(LOC1, body1bin, tol = everybin - 1)] <-
        body1bin + .5
      LOC2 <-
        seq(body2bin,
            by = everybin,
            length.out = (body2bin / everybin) + 1)
      LOC2[near(LOC2, tesbin, tol = everybin - 1)] <- tesbin + .5
      LOC2[near(LOC2, body2bin, tol = everybin - 1)] <-
        body2bin + .5
      LOCname1 <-
        rev(seq((body1bin - tssbin) * binbp,
                by = -everybp,
                length.out = (body1bin / everybin) + 1
        ))
      LOCname1[near(LOC1, tssbin, tol = everybin - 1)] <- "TSS"
      LOCname2 <-
        abs(seq(
          -(tesbin - body2bin) * binbp,
          by = everybp,
          length.out = (body2bin / everybin) + 1
        ))
      LOCname2[near(LOC2, tesbin, tol = everybin - 1)] <- "TES"
      use_plot_breaks <- c(LOC1, LOC2)
      use_plot_breaks_labels <- c(LOCname1, LOCname2)
      use_virtical_line <-
        c(tssbin, tesbin, body1bin, body2bin) + .5
    } else if (mytype == "5'") {
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
      use_plot_breaks[near(use_plot_breaks, tssbin, tol = everybin - 1)] <-
        tssbin + .5
      use_plot_breaks_labels <-
        seq(-tssbin * binbp,
            by = everybp,
            length.out = length(use_plot_breaks))
      use_plot_breaks_labels[near(use_plot_breaks, tssbin, tol = everybin -
                                    1)] <- "TSS"
      use_virtical_line <- c(tssbin, NA, NA, NA) + .5
    } else if (mytype == "3'") {
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
      use_plot_breaks[near(use_plot_breaks, tesbin, tol = everybin - 1)] <-
        tesbin + .5
      use_plot_breaks_labels <-
        abs(seq(
          -(tesbin - totbins) * binbp,
          by = everybp,
          length.out = length(use_plot_breaks)
        ))
      use_plot_breaks_labels[near(use_plot_breaks, tesbin, tol = everybin -
                                    1)] <- "TES"
      use_virtical_line <- c(NA, tesbin, NA, NA) + .5
    } else if (mytype == "none") {
      use_plot_breaks <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin)+1)
      use_plot_breaks_labels <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin)+1)
      use_virtical_line <- c(NA, NA, NA, NA) + .5
    } 
  } else {
    if (mytype == "543") {
      use_plot_breaks <- c(tssbin, tesbin, body1bin, body2bin) + .5
      use_plot_breaks_labels <- c("TSS", "TES", "5|4", "4|3")
      use_virtical_line <-
        c(tssbin, tesbin, body1bin, body2bin) + .5
    } else if (mytype == "5'") {
      use_plot_breaks <- tssbin + .5
      use_plot_breaks_labels <- "TSS"
      use_virtical_line <- c(tssbin, NA, NA, NA) + .5
    } else if (mytype == "3'") {
      use_plot_breaks <- tesbin + .5
      use_plot_breaks_labels <- "TES"
      use_virtical_line <- c(NA, tesbin, NA, NA) + .5
    } else if (mytype == "none") {
      use_plot_breaks <- .5
      use_plot_breaks_labels <- "none"
      use_virtical_line <- c(NA, NA, NA, NA) + .5
    } 
  }
  
  use_virtical_line_color <- c("green", "red", "black", "black")
  
  use_virtical_line_type <- c(1, 1, 1, 1)
  use_plot_breaks <- na_if(use_plot_breaks, 0.5)
  use_virtical_line <- na_if(use_virtical_line, 0.5)
  use_plot_breaks_labels <-
    use_plot_breaks_labels[!is.na(use_plot_breaks)]
  use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
  
  
  use_virtical_line_type <-
    use_virtical_line_type[!is.na(use_virtical_line)]
  use_virtical_line_color <-
    use_virtical_line_color[!is.na(use_virtical_line)]
  use_virtical_line <- use_virtical_line[!is.na(use_virtical_line)]
  list(
    myline = virtical_line_data_frame <- data.frame(
      use_virtical_line,
      use_virtical_line_type,
      use_virtical_line_color,
      stringsAsFactors = FALSE
    ),
    mybrakes = use_plot_breaks,
    mylables = use_plot_breaks_labels
  )
}

# lines and labels preset helper
LinesLablesPreSet <- function(mytype) {
  # 5|4, 4|3, tss, tes, bp/bin, every bin
  if (mytype == "543 bins 20,20,40") {
    tt <- c(20, 40, 15, 45, 100, 5)
  } else if (mytype == "543 bins 10,10,10") {
    tt <- c(10, 20, 5, 25, 100, 5)
  } else if (mytype == "5' 1k 1k 80bins") {
    tt <- c(0, 0, 40, 0, 25, 20)
  } else if (mytype == "5' .25k 10k 205bins") {
    tt <- c(0, 0, 5, 0, 50, 6)
  } else if (mytype == "3'") {
    tt <- c(0, 0, 0, 40, 25, 20)
  } else{
    tt <- c(0, 0, 15, 45, 100, 5)
  }
  tt
}

# help get min and max from apply math data set
MyXSetValues <- function(apply_math, xBinRange) {
  tt <- group_by(apply_math, set) %>%
    filter(bin %in% xBinRange[1]:xBinRange[2]) %>%
    ungroup() %>%
    summarise(min(value, na.rm = T), max(value, na.rm = T)) %>%
    unlist(., use.names = FALSE)
  tt <- round(c(tt, tt[1] - tt[1] * .1, tt[2] + tt[2] * .1), 4)
}

# main ggplot function
GGplotLineDot <-
  function(list_long_data_frame,
           xBinRange,
           plot_options,
           yBinRange,
           line_list,
           use_smooth,
           use_y_label) {
    print("ggplot")
    
    use_col <- plot_options$mycol
    use_dot <- plot_options$mydot
    use_line <- plot_options$myline
    names(use_col) <- plot_options$set
    names(use_dot) <- plot_options$set
    names(use_line) <- plot_options$set
    legend_space <- max(1, (lengths(strsplit(
      plot_options$set, "\n"
    ))))
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
      theme(axis.title.y = element_text(size =  15, margin = margin(2, 10, 2, 2))) +
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

# main ggplot function
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