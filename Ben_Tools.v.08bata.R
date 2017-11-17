# Created by Benjamin Erickson BBErickson@gmail.com

# program for loading packages ----
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
  # library("devtools")
  # install_version("shinyWidgets", version = "0.3.4")
}

# run load needed pakages using my_pakages(x) ----
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

# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 50MB. ----
options(shiny.maxRequestSize = 50 * 1024 ^ 2)

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

# LIST_DATA ----
LIST_DATA <<- list(
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

# Brewer color sets to be avalible ----
kBrewerList <-
  c("Accent",
    "Dark2",
    "Paired",
    "Pastel1",
    "Pastel2",
    "Set1",
    "Set2",
    "Set3")

# color Brewer set that is active to use in plot ----
kListColorSet <- brewer.pal(8, kBrewerList[6])

# math options avalible ----
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
    if (length(file_name) == 1 &&
        length(grep(".url", file_name)) == 1) {
      file_path <- read_lines(file_path)
      file_name <- NULL
      if (load_gene_list) {
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
        strsplit(strsplit(as.character(file_name[x]), '.tab')[[1]][1], '\\._')[[1]][1]
      if(grepl("^543\\.", legend_nickname)){
        legend_nickname <- strsplit(legend_nickname, "^543\\.")[[1]][2]
      } else if (grepl("^5\\.", legend_nickname)){
        legend_nickname <- strsplit(legend_nickname, "^5\\.")[[1]][2]
      } else if (grepl("^4\\.", legend_nickname)){
        legend_nickname <- strsplit(legend_nickname, "^4\\.")[[1]][2]
      } else if (grepl("^3\\.", legend_nickname)){
        legend_nickname <- strsplit(legend_nickname, "^3\\.")[[1]][2]
      }
      
      if (any(legend_nickname == names(list_data$table_file)) |
          any(legend_nickname == names(list_data$gene_info))) {
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
        count_fields(file_path[x],
                     n_max = 1,
                     tokenizer = tokenizer_tsv())
      
      if (num_bins == 1 & load_gene_list) {
        tablefile <-
          suppressMessages(read_tsv(
            file_path,
            col_names = "gene",
            comment = "#",
            cols(gene = col_character())
          ))
      } else {
        if (num_bins > 6) {
          tablefile <- suppressMessages(
            read_tsv (
              file_path[x],
              comment = "#",
              col_names = c("gene", 1:(num_bins - 1)),
              skip = 1
            ) %>%
              gather(., bin, score, 2:(num_bins))
          ) %>%
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
        all_empty_bin <-
          group_by(tablefile, bin) %>% summarise(mytest = sum(score) == 0) %>% pull(mytest)
        if (any(all_empty_bin)) {
          showModal(
            modalDialog(
              title = "Information message",
              "!!!! All scores for at least one bin are empty !!!!",
              size = "s",
              easyClose = TRUE
            )
          )
        }
        if (file_count > 0 &
            num_bins != list_data$x_plot_range[2]) {
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
        list_data$gene_file[[legend_nickname]]$info <-
          paste("Loaded gene list from file",
                legend_nickname,
                Sys.Date())
        if (length(list_data$gene_file) > 5) {
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
              mycol = RgbToHex(my_hex = list_data$gene_info[[1]][[i]]$mycol,
                               tint = mytint),
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
          enesg <-
            inner_join(list_data$gene_file[[g]]$full,
                       list_data$gene_file[[1]]$use,
                       by = "gene")
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
          my_name_g <-
            sub(
              "([0-9]+)",
              n_distinct(list_data$gene_file[[g]]$full$gene),
              names(list_data$gene_file)[g]
            )
          names(list_data$gene_file)[g] <<- my_name_g
          names(list_data$gene_info)[g] <<- my_name_g
          if (length(list_data$gene_file) > 5) {
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
    if (!gbyg) {
      myname <- "mean_norm"
      mynom <-
        group_by(mynom, bin, set) %>% mutate(score = mean(score, na.rm = TRUE)) %>% ungroup()
      mydom <-
        group_by(mydom, bin, set) %>% mutate(score = mean(score, na.rm = TRUE)) %>% ungroup()
    }
    if (nodivzero) {
      myname <- paste0(myname, "-0_min/2")
      # find min value /2 to replace 0s
      new_gene_list <-
        inner_join(mynom, mydom, by = c("gene", "bin")) %>%
        na_if(0)
      new_min_for_na <-
        min(c(new_gene_list$score.x, new_gene_list$score.y),
            na.rm = TRUE) / 2
      # replace 0's with min/2
      new_gene_list <-
        replace_na(new_gene_list,
                   list(score.y = new_min_for_na, score.x = new_min_for_na))
    } else {
      new_gene_list <- inner_join(mynom, mydom, by = c("gene", "bin"))
      new_gene_list <-
        group_by(new_gene_list, gene) %>%
        summarise(test = sum(score.x, score.y)) %>%
        filter(!is.na(test)) %>%
        semi_join(new_gene_list, ., by = "gene")
    }
    legend_nickname <- paste0(nom, "/\n", dnom, ":\n", myname)
    color_safe <-
      (length(list_data$table_file) + 1) %% length(kListColorSet)
    if (color_safe == 0) {
      color_safe <- 1
    }
    gene_names <-
      semi_join(list_data$gene_file[[1]]$use, new_gene_list, by = "gene")
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
    list_data$gene_file[[my_name]]$use <-
      semi_join(list_data$gene_file[[my_name]]$use,
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
      enesg <-
        inner_join(list_data$gene_file[[g]]$full,
                   list_data$gene_file[[1]]$use,
                   by = "gene")
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
      my_name_g <-
        sub(
          "([0-9]+)",
          n_distinct(list_data$gene_file[[g]]$full$gene),
          names(list_data$gene_file)[g]
        )
      names(list_data$gene_file)[g] <<- my_name_g
      names(list_data$gene_info)[g] <<- my_name_g
      if (length(list_data$gene_file) > 5) {
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
RemoveFile <- function(list_data, file_name, remove_all) {
  if (length(list_data$table_file) > 1 & !remove_all) {
    # remove tool gene lists ? TODO
    sapply(names(list_data$gene_file), function(g) {
      list_data$gene_file[[g]][[file_name]] <<- NULL
      list_data$gene_info[[g]][[file_name]] <<- NULL
    })
    list_data$table_file[[file_name]] <- NULL
    list_data$gene_file[[1]]$use <-
      distinct(list_data$table_file[[1]], gene)
    sapply(list_data$table_file[-1], function(i) {
      list_data$gene_file[[1]]$use <<-
        semi_join(list_data$gene_file[[1]]$use, i, by = "gene")
    })
    my_name <-
      paste("common\nn =", n_distinct(list_data$gene_file[[1]]$use$gene))
    names(list_data$gene_file)[1] <- my_name
    names(list_data$gene_info)[1] <- my_name
    
    sapply(seq_along(list_data$gene_file)[-1], function(g) {
      list_data$gene_file[[g]]$full <<-
        inner_join(list_data$gene_file[[g]]$full,
                   list_data$gene_file[[1]]$use,
                   by = "gene")
      list_data$gene_file[[g]]$use <<-
        select(list_data$gene_file[[g]]$full, gene)
      my_name_g <-
        sub(
          "([0-9]+)",
          n_distinct(list_data$gene_file[[g]]$full$gene),
          names(list_data$gene_file)[g]
        )
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
IntersectGeneLists <- function(list_data, list_name) {
  if (is.null(list_name)) {
    return(NULL)
  }
  setProgress(1, detail = paste("building list"))
  outlist <- NULL
  lapply(list_name, function(j) {
    outlist[[j]] <<- list_data$gene_file[[j]]$use
  })
  outlist <- bind_rows(outlist)
  for (rr in grep("Gene_List_", names(LIST_DATA$gene_file), value = T)) {
    if (length(rr) > 0) {
      list_data$gene_file[[rr]] <- NULL
      list_data$gene_info[[rr]] <- NULL
    }
  }
  
  nick_name <- NULL
  setProgress(2, detail = paste("building inclusive list"))
  inclusive <- distinct(outlist)
  if (n_distinct(inclusive$gene) > 0) {
    nick_name1 <-
      paste("Gene_List_inclusive\nn =", n_distinct(inclusive$gene))
    nick_name <- c(nick_name, nick_name1)
    list_data$gene_file[[nick_name1]]$full <- inclusive
    list_data$gene_file[[nick_name1]]$use <- select(inclusive, gene)
    list_data$gene_file[[nick_name1]]$info <-
      paste("Gene_List_inclusive",
            "from",
            paste(list_name, collapse = " and "),
            Sys.Date())
  }
  setProgress(3, detail = paste("building intersect list"))
  intersect <- filter(outlist, duplicated(gene))
  if (n_distinct(intersect$gene) > 0) {
    nick_name1 <-
      paste("Gene_List_intersect\nn =", n_distinct(intersect$gene))
    nick_name <- c(nick_name, nick_name1)
    list_data$gene_file[[nick_name1]]$full <- intersect
    list_data$gene_file[[nick_name1]]$use <- select(intersect, gene)
    list_data$gene_file[[nick_name1]]$info <-
      paste("Gene_List_intersect",
            "from",
            paste(list_name, collapse = " and "),
            Sys.Date())
    
    setProgress(4, detail = paste("building exclusive list"))
    exclusive <- anti_join(inclusive, intersect, by = "gene")
    if (n_distinct(exclusive$gene) > 0) {
      nick_name1 <-
        paste("Gene_List_exclusive\nn =", n_distinct(exclusive$gene))
      nick_name <- c(nick_name, nick_name1)
      list_data$gene_file[[nick_name1]]$full <- exclusive
      list_data$gene_file[[nick_name1]]$use <-
        select(exclusive, gene)
      list_data$gene_file[[nick_name1]]$info <-
        paste("Gene_List_exclusive",
              "from",
              paste(list_name, collapse = " and "),
              Sys.Date())
    }
  }
  setProgress(5, detail = "finishing up")
  if (length(list_data$gene_file) > 5) {
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
          mycol = RgbToHex(my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nn)]][[i]]$mycol,
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
        mutate(myper = as.numeric(strtrim(cume_dist(mysums), 5))) %>%
        arrange(desc(mysums))
      
      gene_count <- nrow(apply_bins)
      
      if (topbottom == "Top%") {
        num2 <- c(1, ceiling(gene_count * (num / 100)))
        topbottom <- paste(topbottom, paste0(num, "%"))
      } else if (topbottom == "Middle%") {
        num2 <-
          c(count(apply_bins, myper >= max(.5, num / 100))[[2]][2],
            count(apply_bins, myper <= min(.5, (100 - num) / 100))[[2]][1])
      } else {
        num2 <-
          c(ceiling((gene_count + 1) - (gene_count * (num / 100))), gene_count)
        topbottom <- paste(topbottom, paste0(num, "%"))
      }
      nickname <- list_data$gene_info[[1]][[j]]$set
      outlist2 <- mutate(apply_bins,!!nickname := myper) %>%
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
    if (length(outlist$gene) == 0) {
      return(NULL)
    }
    setProgress(lc + 2, detail = "building list")
    nick_name <-
      strtrim(gsub(
        "(.{30})",
        "\\1... ",
        paste0("Sort\nn = ", n_distinct(outlist$gene), "-", list_name)
      ), 33)
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- select(outlist, gene)
    list_data$gene_file[[nick_name]]$info <-
      paste(
        "Sort",
        topbottom,
        "bins",
        start_bin,
        "to",
        end_bin,
        "from",
        list_name,
        paste(file_names, collapse = " "),
        Sys.Date()
      )
    if (length(list_data$gene_file) > 5) {
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
                 mycol = RgbToHex(my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nick_name)]][[i]]$mycol,
                                  tint = mytint),
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
      df <-
        semi_join(list_data$table_file[[j]], enesg, by = 'gene') %>%
        group_by(gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T)) %>%
        na_if(0) %>% ungroup()
      if (nodivzero) {
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
    
    for (rr in grep("Ratio_", names(LIST_DATA$gene_file), value = T)) {
      if (length(rr) > 0) {
        list_data$gene_file[[rr]] <- NULL
        list_data$gene_info[[rr]] <- NULL
      }
    }
    
    nick_name <- NULL
    setProgress(2, detail = paste("building list", ratio1file))
    upratio <- filter(outlist[[1]], Ratio < 1 / num)
    if (n_distinct(upratio$gene) > 0) {
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
    if (n_distinct(upratio$gene) > 0) {
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
    if (n_distinct(upratio$gene) > 0) {
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
    if (length(list_data$gene_file) > 5) {
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
            mycol = RgbToHex(my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nn)]][[i]]$mycol,
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
ClusterNumList <- function(list_data,
                           list_name,
                           clusterfile,
                           start_bin,
                           end_bin,
                           num,
                           myname) {
  if (is_empty(list_data$clust)) {
    return(NULL)
  }
  setProgress(3, detail = "spliting into clusters")
  for (rr in grep("Cluster_", names(list_data$gene_file), value = T)) {
    if (length(rr) > 0) {
      list_data$gene_file[[rr]] <- NULL
      list_data$gene_info[[rr]] <- NULL
    }
  }
  for (rr in grep("Group_", names(list_data$gene_file), value = T)) {
    if (length(rr) > 0) {
      list_data$gene_file[[rr]] <- NULL
      list_data$gene_info[[rr]] <- NULL
    }
  }
  if (myname == "Cluster_") {
    gene_list <-
      mutate(list_data$clust$use, cm = cutree(list_data$clust$cm, num))
  } else {
    gene_list <-
      mutate(list_data$clust$full, cm = ntile(cm, as.numeric(num)))
  }
  for (nn in 1:num) {
    outlist <- filter(gene_list, cm == nn)
    nick_name <-
      paste(paste0(myname, nn, "\nn ="), n_distinct(outlist$gene))
    list_data$gene_file[[nick_name]]$full <- outlist
    list_data$gene_file[[nick_name]]$use <- select(outlist, gene)
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
    setProgress(4, detail = paste("finishing cluster", nn))
    if (length(list_data$gene_file) > 5) {
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
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
                           tint = mytint),
          onoff = 0,
          rnorm = "1"
        ))
  }
  list_data$STATE[2] <-
    grep(paste0(myname, "1"), names(list_data$gene_file), value = T)
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
  list_data$clust$cm <-
    hclust.vector(as.data.frame(spread(df, bin, score))[, c((start_bin:end_bin) + 2)], method = "ward")
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
    gene_count <-
      n_distinct(list_data$gene_file[[list_name]]$use$gene)
    num <-
      c(ceiling(gene_count * bottom_per / 100),
        ceiling(gene_count * top_per / 100))
    outlist <- NULL
    lapply(cdffile, function(j) {
      df <-
        semi_join(list_data$table_file[[j]], list_data$gene_file[[list_name]]$use, by = 'gene') %>%
        group_by(gene) %>%
        summarise(sum1 = sum(score[start1_bin:end1_bin],	na.rm = T),
                  sum2 = sum(score[start2_bin:end2_bin],	na.rm = T)) %>%
        na_if(0) %>% ungroup()
      if (nodivzero) {
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
      outlist[[j]] <<-
        transmute(df, gene = gene, value = sum1 / sum2) %>%
        na_if(Inf) %>%
        replace_na(list(value = 0)) %>%
        arrange(desc(value)) %>%
        mutate(bin = row_number(), set = j)
    })
    
    outlist <- bind_rows(outlist)
    for (rr in grep("CDF\nn", names(LIST_DATA$gene_file), value = T)) {
      if (length(rr) > 0) {
        list_data$gene_file[[rr]] <- NULL
        list_data$gene_info[[rr]] <- NULL
      }
    }
    
    setProgress(2, detail = paste("building list"))
    gene_list <-
      group_by(outlist, gene) %>% filter(all(between(bin, num[1], num[2]))) %>%
      distinct(gene) %>%
      ungroup()
    if (n_distinct(gene_list$gene) > 0) {
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
    if (length(list_data$gene_file) > 5) {
      mytint <- 0
    } else{
      mytint <- length(list_data$gene_file) * 0.1
    }
    list_data$gene_info[[nick_name1]] <-
      lapply(setNames(names(list_data$gene_info[[1]]),
                      names(list_data$gene_info[[1]])),
             function(i)
               tibble(
                 set = i,
                 mydot = kDotOptions[1],
                 myline = kLineOptions[1],
                 mycol = RgbToHex(my_hex = list_data$gene_info[[sum(names(list_data$gene_info) != nick_name1)]][[i]]$mycol,
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
    setProgress(1, detail = paste("Gathering info"))
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
      setProgress(1+ length(list_data_frame), detail = paste("applying math to ", i))
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
      if (relative_frequency == "rel gene frequency") {
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
  if (tssbin == 0 & tesbin > 0) {
    mytype <- "3'"
  } else if (tesbin == 0 & tssbin > 0) {
    mytype <- "5'"
  } else if (tssbin > 0 &
             tesbin > 0 & body1bin > 0 & body2bin > 0) {
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
            length.out = (totbins / everybin) + 1)
      use_plot_breaks_labels <-
        seq(1,
            by = everybin,
            length.out = (totbins / everybin) + 1)
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
MyXSetValues <- function(apply_math, xBinRange, yBinRange) {
  tt <- group_by(apply_math, set) %>%
    filter(bin %in% xBinRange[1]:xBinRange[2]) %>%
    ungroup() %>%
    summarise(min(value, na.rm = T), max(value, na.rm = T)) %>%
    unlist(., use.names = FALSE)
  tt <-
    c(tt[1] + (tt[1] * (yBinRange[1] / 100)), tt[2] + (tt[2] * ((yBinRange[2] -
                                                                   100) / 100)))
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

# server ----
server <- function(input, output, session) {
  # remove on non-local deployment
  session$onSessionEnded(stopApp)
  
  # reacive values ----
  reactive_values <- reactiveValues(
    pickerfile_controler = "",
    Y_Axis_Lable = NULL,
    Y_Axis_numbers = NULL,
    Lines_Lables_List = NULL,
    Apply_Math = NULL,
    Plot_Options = NULL,
    Plot_controler = NULL,
    Picker_controler = NULL,
    Y_Axis_plot = 0
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
    toggle("showcdftoolpicker",
           condition = (input$tabs == "cdftool" &
                          LIST_DATA$STATE[1] != 0))
    toggle("showgenelistspicker",
           condition = (input$tabs == "genelists" &
                          length(LIST_DATA$gene_file) > 1))
    if (input$tabs == "filenorm" & LIST_DATA$STATE[1] != 0) {
      updatePickerInput(
        session,
        "pickernumerator",
        choices = names(LIST_DATA$table_file),
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[1]], "[[", 4)
        ), sep = ":"))
      )
      updatePickerInput(
        session,
        "pickerdenominator",
        choices = names(LIST_DATA$table_file),
        choicesOpt = list(style = paste("color", c(
          sapply(LIST_DATA$gene_info[[1]], "[[", 4)
        ), sep = ":"))
      )
      if (LIST_DATA$STATE[4] != 0) {
        LIST_DATA$STATE[4] <<- 3
      }
    }
    if (input$tabs == "genelists" &
        length(LIST_DATA$gene_file) != 0) {
      hide('actiongenelistsdatatable')
      if (length(LIST_DATA$gene_file) > 1) {
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
      if (LIST_DATA$STATE[4] != 0) {
        LIST_DATA$STATE[4] <<- 3
      }
    }
    if (input$tabs == "sorttool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectsortfile
      og <- input$pickersortfile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- names(LIST_DATA$gene_file)[1]
      } else if (!all(og %in% names(LIST_DATA$table_file)) |
                 LIST_DATA$STATE[3] == "Sort\nn") {
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
      if (LIST_DATA$STATE[4] != 0) {
        LIST_DATA$STATE[4] <<- 3
      }
    }
    
    if (input$tabs == "ratiotool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectratiofile
      og <- c(input$pickerratio1file, input$pickerratio2file)
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- names(LIST_DATA$gene_file)[1]
      } else if (!all(og %in% names(LIST_DATA$table_file)) |
                 LIST_DATA$STATE[3] == "Ratio_") {
        og <- c("", "")
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
          value = c(0, 0)
        )
        hide('actionratiodatatable')
      }
      if (LIST_DATA$STATE[4] != 0) {
        LIST_DATA$STATE[4] <<- 3
      }
    }
    
    if (input$tabs == "clustertool" & LIST_DATA$STATE[1] != 0) {
      ol <- input$selectclusterfile
      og <- input$pickerclusterfile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- names(LIST_DATA$gene_file)[1]
      } else if (!all(og %in% names(LIST_DATA$table_file)) |
                 LIST_DATA$STATE[3] == "Cluste") {
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
      if (LIST_DATA$STATE[4] != 0) {
        LIST_DATA$STATE[4] <<- 3
      }
    }
    
    if (input$tabs == "cdftool" & LIST_DATA$STATE[1] != 0) {
      ol1 <- input$selectcdffile1
      og1 <- input$pickercdffile1
      if (!ol1 %in% names(LIST_DATA$gene_file)) {
        ol1 <- names(LIST_DATA$gene_file)[1]
      } else if (!all(og1 %in% names(LIST_DATA$table_file)) |
                 LIST_DATA$STATE[3] == "CDF\nn") {
        og1 <- ""
      }
      updateSelectInput(
        session,
        "selectcdffile1",
        choices = names(LIST_DATA$gene_file),
        selected = ol1
      )
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
          value = c(
            LIST_DATA$x_plot_range[1],
            floor(LIST_DATA$x_plot_range[2] / 3)
          )
        )
        updateSliderInput(
          session,
          "sliderbincdf2",
          min = LIST_DATA$x_plot_range[1],
          max = LIST_DATA$x_plot_range[2],
          value = c(
            ceiling(LIST_DATA$x_plot_range[2] / 3),
            LIST_DATA$x_plot_range[2]
          )
        )
        
        hide('actioncdfdatatable')
      }
      if (LIST_DATA$STATE[4] != 0) {
        LIST_DATA$STATE[4] <<- 3
      }
    }
    
    toggle(
      "selectlineslablesshow",
      condition = (input$tabs == "mainplot" &
                     LIST_DATA$STATE[1] != 0)
    )
    # first time switch tab auto plot
    if (input$tabs == "mainplot" & LIST_DATA$STATE[1] != 0) {
      reactive_values$Picker_controler <-
        sapply(LIST_DATA, function(i)
          (names(i)))
      if (LIST_DATA$STATE[4] == 0) {
        withProgress(message = 'Calculation in progress',
                     detail = 'This may take a while...',
                     value = 0,
                     {
        reactive_values$Apply_Math <-
          ApplyMath(
            LIST_DATA,
            input$myMath,
            input$radioplotnrom,
            as.numeric(input$sliderplotBinNorm)
          )
                     })
        if (!is.null(reactive_values$Apply_Math)) {
          reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
          enable("showmainplot")
          print(LIST_DATA$STATE)
          LIST_DATA$STATE[1] <<- 1
          LIST_DATA$STATE[4] <<- 1
          print(LIST_DATA$STATE)
        } else{
          disable("showmainplot")
        }
      } else {
        toggle(
          "actionmyplotshow",
          condition = (input$tabs == "mainplot" &
                         LIST_DATA$STATE[1] == 2)
        )
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
          reactive_values$Y_Axis_numbers,
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
                 if (LIST_DATA$STATE[1] == 0) {
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
  observeEvent(input$actionnormfactor, ignoreInit = TRUE, {
    if (!is.na(input$normfactor) &
        !input$normfactor %in% c(0, 1) &
        input$normfactor != LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"]) {
      print("norm")
      LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"] <<-
        as.character(input$normfactor)
      LIST_DATA$table_file[[input$radiodataoption]] <<-
        mutate(LIST_DATA$table_file[[input$radiodataoption]],
               score = score / as.numeric(input$normfactor))
      LIST_DATA$STATE[4] <<- 0
    } else if (!is.na(input$normfactor)) {
      updateNumericInput(session,
                         "normfactor",
                         value = as.numeric(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["rnorm"]))
    }
  })
  
  # record new nickname and norm factor ---- 
  observeEvent(input$actionoptions, ignoreInit = TRUE, {
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
    
    LIST_DATA$STATE[4] <<- 0
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
              reactive_values$Y_Axis_numbers,
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
              reactive_values$Y_Axis_numbers,
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
              reactive_values$Y_Axis_numbers,
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
  
  # records check box on/off ----
  observeEvent(reactive_values$picker,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 print(LIST_DATA$STATE)
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
                       selectgenelistonoff <-
                         gsub("-bensspace2-", " ", gsub("-bensspace1-", "\n", i))
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
                     reactive_values$Apply_Math <- NULL
                     disable("showmainplot")
                   } else {
                     LIST_DATA$STATE[4] <<- 2
                   }
                 }
               })
  
  # plots when action button is pressed ----
  observeEvent(input$actionmyplot, ignoreInit = TRUE, {
    print("plot button")
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   reactive_values$Apply_Math <-
      ApplyMath(
        LIST_DATA,
        input$myMath,
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
      )
                 })
    if (!is.null(reactive_values$Apply_Math)) {
      reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
      enable("showmainplot")
    } else{
      disable("showmainplot")
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
    reactive_values$Y_Axis_numbers <-
      MyXSetValues(reactive_values$Apply_Math,
                   input$sliderplotBinRange,
                   input$sliderplotYRange)
    updateSliderInput(session,
                      "sliderplotYRange",
                      value = c(0, 100))
    # Forces update if y Asis values stay the same
    if (all(c(0, 100) == input$sliderplotYRange)) {
      reactive_values$Y_Axis_plot <- reactive_values$Y_Axis_plot + 1
    }
  })
  
  # renders plot ----
  output$plot <- renderPlot({
    reactive_values$Plot_controler
  })
  
  # updates norm applymath ----
  observeEvent(c(input$myMath,
                 input$sliderplotBinNorm,
                 input$radioplotnrom),
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
                   withProgress(message = 'Calculation in progress',
                                detail = 'This may take a while...',
                                value = 0,
                                { reactive_values$Apply_Math <-
                     ApplyMath(
                       LIST_DATA,
                       input$myMath,
                       input$radioplotnrom,
                       as.numeric(input$sliderplotBinNorm)
                     )
                                })
                 }
               })
  
  # y slider is trigger ----
  observeEvent(input$sliderplotYRange, ignoreInit = T, {
    print("y slider")
    reactive_values$Y_Axis_numbers <-
      MyXSetValues(reactive_values$Apply_Math,
                   input$sliderplotBinRange,
                   input$sliderplotYRange)
    reactive_values$Y_Axis_plot <- reactive_values$Y_Axis_plot + 1
  })
  
  # plots when bin slider or other triggers is triggered ----
  observeEvent(
    c(
      reactive_values$Lines_Lables_List,
      input$sliderplotBinRange,
      reactive_values$Y_Axis_plot
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
            reactive_values$Y_Axis_numbers,
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
        reactive_values$Y_Axis_numbers,
        reactive_values$Lines_Lables_List,
        input$checkboxsmooth,
        reactive_values$Y_Axis_Lable
      )
  })
  
  # quick color set change ----
  observeEvent(input$kbrewer, ignoreInit = TRUE, {
    print("kbrewer")
    if (!is.null(LIST_DATA$gene_info[[1]]) &
        input$kbrewer != "select") {
      kListColorSet <<- brewer.pal(8, input$kbrewer)
      common_name <- names(LIST_DATA$gene_info)[1]
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
      updateColourInput(session, "colourhex", value =
                          paste(LIST_DATA$gene_info[[input$selectgenelistoptions]][[input$radiodataoption]]["mycol"]))
      if (!is.null(reactive_values$Apply_Math)) {
        reactive_values$Plot_Options <- MakePlotOptionFrame(LIST_DATA)
        reactive_values$Plot_controler <-
          GGplotLineDot(
            reactive_values$Apply_Math,
            input$sliderplotBinRange,
            reactive_values$Plot_Options,
            reactive_values$Y_Axis_numbers,
            reactive_values$Lines_Lables_List,
            input$checkboxsmooth,
            reactive_values$Y_Axis_Lable
          )
      }
      updateSelectInput(session, "kbrewer", selected = "select")
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
      output$DynamicGenePicker <- renderUI({
        pickerlist
      })
    }
  )
  
  # Remove data file ----
  observeEvent(input$actionremovefile, ignoreInit = TRUE, {
    print("remove file")
    LIST_DATA <<-
      RemoveFile(LIST_DATA,
                 input$radiodataoption,
                 input$checkboxremovefile)
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
      updateSelectInput(session,
                        "selectgenelistoptions",
                        choices = 'Load Data File')
      hide("filegene1")
      hide("checkboxconvert")
      hide("downloadGeneList")
      hide("checkboxsavesplit")
      hide("filecolor")
      hide("showmainplot")
      hide("startoff")
      shinyjs::addClass(selector = "body", class = "sidebar-collapse")
    }
    
    reactive_values$Picker_controler <- names(LIST_DATA$table_file)
  })
  
  # Remove gene list ----
  observeEvent(input$actionremovegene, ignoreInit = TRUE, {
    print("remove gene list")
    LIST_DATA <<-
      RemoveGeneList(LIST_DATA, input$selectgenelistoptions)
    updateSelectInput(
      session,
      "selectgenelistoptions",
      choices = names(LIST_DATA$gene_info),
      selected = LIST_DATA$STATE[2]
    )
    
  })
  
  # create norm file ----
  observeEvent(input$actionnorm, ignoreInit = TRUE, {
    LIST_DATA <<- MakeNormFile(
      LIST_DATA,
      input$pickernumerator,
      input$pickerdenominator,
      input$checkboxnormmean,
      input$checkboxnormzero
    )
    updatePickerInput(session,
                      "pickernumerator", selected = "")
    updatePickerInput(session,
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
  observeEvent(input$actiongenelists,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 if (is.null(input$actiongenelists)) {
                   hide('genelists1table')
                   hide('genelists2table')
                   hide('genelists3table')
                 }
               })
  
  
  # Gene action ----
  observeEvent(input$actiongenelists, {
    print("gene lists action")
    hide('actiongenelistsdatatable')
    hide('genelists1table')
    hide('genelists2table')
    hide('genelists3table')
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
                   LD <- IntersectGeneLists(LIST_DATA,
                                            input$pickergenelists)
                 })
    if (!is_empty(LD$table_file)) {
      LIST_DATA <<- LD
      glo <- input$selectgenelistoptions
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
      show('actiongenelistsdatatable')
    } else {
      return()
    }
  })
  
  # Gene lists show gene list ----
  observeEvent(input$actiongenelistsdatatable, ignoreInit = TRUE, {
    print("generiate gene lists table")
    hide('actiongenelistsdatatable')
    if (any(grep("Gene_List_intersect\nn =", names(LIST_DATA$gene_info)) >
            0)) {
      newnames1 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_intersect\nn =",
               names(LIST_DATA$gene_info),
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
                                                     names(LIST_DATA$gene_info))]]$full,
                           rownames = FALSE,
                           colnames = newnames1,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Gene_List_intersect\nn =",
                                                               names(LIST_DATA$gene_info))]]$info,
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
      show('genelists1table')
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
      mytab <- "Inclusive Gene Lists"
    }
    if (any(grep("Gene_List_inclusive\nn =", names(LIST_DATA$gene_info)) >
            0)) {
      newnames2 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_inclusive\nn =",
               names(LIST_DATA$gene_info),
               value = TRUE
             ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$genelists2table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Gene_List_inclusive\nn =",
                                                     names(LIST_DATA$gene_info))]]$full,
                           rownames = FALSE,
                           colnames = newnames2,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Gene_List_inclusive\nn =",
                                                               names(LIST_DATA$gene_info))]]$info,
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
      show('genelists2table')
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
      if (mytab == "Inclusive Gene Lists") {
        mytab <- "Exclusive Gene Lists"
      }
    }
    if (any(grep("Gene_List_exclusive\nn =", names(LIST_DATA$gene_info)) >
            0)) {
      newnames3 <-
        gsub("\n",
             " ",
             grep(
               "Gene_List_exclusive\nn =",
               names(LIST_DATA$gene_info),
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
                                                     names(LIST_DATA$gene_info))]]$full,
                           rownames = FALSE,
                           colnames = newnames3,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Gene_List_exclusive\nn =",
                                                               names(LIST_DATA$gene_info))]]$info,
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
      show('genelists3table')
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
  observeEvent(input$pickersortfile,
               ignoreNULL = FALSE,
               ignoreInit = TRUE,
               {
                 if (is.null(input$pickersortfile)) {
                   hide('sorttable')
                 } else if (input$pickersortfile[1] == "") {
                   hide('sorttable')
                 }
               })
  
  # sort tool action ----
  observeEvent(input$actionsorttool, {
    print("sort tool")
    hide('actionsortdatatable')
    hide('sorttable')
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
      ol <- input$selectsortfile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- grep("Sort\nn", names(LIST_DATA$gene_file), value = TRUE)
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
    } else {
      return()
    }
  })
  
  # Sort gene list show data table ----
  observeEvent(input$actionsortdatatable, ignoreInit = TRUE, {
    print("show data table")
    if (any(grep("Sort\nn", names(LIST_DATA$gene_info)) > 0)) {
      newnames <-
        gsub("(.{20})", "\\1... ", names(LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$full))
      dt <- datatable(
        LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$full,
        rownames = FALSE,
        colnames = strtrim(newnames, 24),
        class = 'cell-border stripe compact',
        filter = 'top',
        caption = LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$info,
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
                "return type === 'display' && data.length > 24 ?",
                "'<span title=\"' + data + '\">' + data.substr(0, 19) + '...</span>' : data;",
                "}"
              )
            )
          )
        )
      ) %>% formatPercentage(names(LIST_DATA$gene_file[[grep("Sort\nn", names(LIST_DATA$gene_info))]]$full)[-1])
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
    hide('actionsortdatatable')
    show('sorttable')
  })
  
  # sort tool gene list $use ----
  observeEvent(input$sorttable_rows_all, ignoreInit = TRUE, {
    newname <-
      strtrim(gsub("(.{30})", "\\1... ", paste0(
        sub(
          "([0-9]+)",
          length(input$sorttable_rows_all),
          LIST_DATA$STATE[2]
        )
      )), 33)
    oldname <- grep("Sort\nn =", names(LIST_DATA$gene_file))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("sort filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$sorttable_rows_all])
      print("sort $use picker")
      glo <- input$selectgenelistoptions
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_file),
        selected = glo
      )
      ol <- input$selectsortfile
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- newname
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
  observeEvent(
    c(input$pickerratio1file, input$pickerratio2file),
    ignoreInit = TRUE,
    ignoreNULL = FALSE,
    {
      if (input$pickerratio1file != "" | input$pickerratio2file != "") {
        
      } else {
        hide('ratio1table')
        hide('ratio2table')
        hide('ratio3table')
      }
    }
  )
  
  # ratio tool gene lists $use ----
  observeEvent(input$ratio1table_rows_all, ignoreInit = TRUE, {
    newname <-
      paste("Ratio_Up_file1\nn =",
            length(input$ratio1table_rows_all))
    oldname <-
      grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("ratio1 filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio1table_rows_all])
      glo <- input$selectgenelistoptions
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    }
  })
  
  observeEvent(input$ratio2table_rows_all, ignoreInit = TRUE, {
    newname <-
      paste("Ratio_Up_file2\nn =",
            length(input$ratio2table_rows_all))
    oldname <-
      grep("Ratio_Up_file2\nn =", names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("ratio2 filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio2table_rows_all])
      glo <- input$selectgenelistoptions
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    }
  })
  
  observeEvent(input$ratio3table_rows_all, ignoreInit = TRUE, {
    newname <-
      paste("Ratio_No_Diff\nn =", length(input$ratio3table_rows_all))
    oldname <-
      grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("no ratio filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$ratio3table_rows_all])
      glo <- input$selectgenelistoptions
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    }
  })
  
  # Ratio tool action ----
  observeEvent(input$actionratiotool, ignoreInit = TRUE, {
    print("ratio tool action")
    hide('ratio1table')
    hide('ratio2table')
    hide('ratio3table')
    if (is.numeric(input$numericratio)) {
      if (input$numericratio < 0) {
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    } else {
      return()
    }
  })
  
  # Ratio show gene list ----
  observeEvent(input$actionratiodatatable, ignoreInit = TRUE, {
    print("generiate ratio table")
    hide('actionratiodatatable')
    if (any(grep("Ratio_Up_file1\nn =", names(LIST_DATA$gene_info)) > 0)) {
      newnames1 <-
        gsub("\n", " ", grep(
          "Ratio_Up_file1\nn =",
          names(LIST_DATA$gene_info),
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
                                                     names(LIST_DATA$gene_info))]]$full,
                           rownames = FALSE,
                           colnames = newnames1,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Ratio_Up_file1\nn", names(LIST_DATA$gene_info))]]$info,
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
      show('ratio1table')
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
    if (any(grep("Ratio_Up_file2\nn =", names(LIST_DATA$gene_info)) > 0)) {
      newnames2 <-
        gsub("\n", " ", grep(
          "Ratio_Up_file2\nn =",
          names(LIST_DATA$gene_info),
          value = TRUE
        ))
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
                     output$ratio2table <-
                       DT::renderDataTable(
                         datatable(
                           LIST_DATA$gene_file[[grep("Ratio_Up_file2\nn =",
                                                     names(LIST_DATA$gene_info))]]$full,
                           rownames = FALSE,
                           colnames = newnames2,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Ratio_Up_file2\nn", names(LIST_DATA$gene_info))]]$info,
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
      show('ratio2table')
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
    if (any(grep("Ratio_No_Diff\nn =", names(LIST_DATA$gene_info)) > 0)) {
      newnames3 <-
        gsub("\n", " ", grep(
          "Ratio_No_Diff\nn =",
          names(LIST_DATA$gene_info),
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
                                                     names(LIST_DATA$gene_info))]]$full,
                           rownames = FALSE,
                           colnames = newnames3,
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep("Ratio_No_Diff\nn", names(LIST_DATA$gene_info))]]$info,
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
      show('ratio3table')
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
  observeEvent(input$pickerclusterfile,
               ignoreInit = TRUE,
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
    newname <-
      paste0(reactive_values$clustergroups,
             "1\nn = ",
             length(input$cluster1table_rows_all))
    oldname <-
      grep(paste0(reactive_values$clustergroups, "1\nn ="),
           names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("cluster1 filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster1table_rows_all])
      LD <- LIST_DATA$gene_info
      
      sapply(names(LD), function(i)
        sapply(names(LD[[i]]), function(j)
          if (i %in% grep(reactive_values$clustergroups, names(LD), value = T) &
              j == input$pickerclusterfile) {
            LD[[i]][[j]][5] <<- input$pickerclusterfile
          } else{
            LD[[i]][[j]][5] <<- 0
          }))
      LIST_DATA$gene_info <- LD
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
      reactive_values$Apply_Cluster_Math <- ApplyMath(
        LIST_DATA,
        input$myMath,
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
      )
                   })
      if (!is.null(reactive_values$Apply_Cluster_Math)) {
        reactive_values$Plot_Cluster_Options <-
          MakePlotOptionFrame(LIST_DATA)
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    }
  })
  
  observeEvent(input$cluster2table_rows_all, ignoreInit = TRUE, {
    newname <-
      paste0(reactive_values$clustergroups,
             "2\nn = ",
             length(input$cluster2table_rows_all))
    oldname <-
      grep(paste0(reactive_values$clustergroups, "2\nn ="),
           names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("cluster2 filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster2table_rows_all])
      LD <- LIST_DATA$gene_info
      sapply(names(LD), function(i)
        sapply(names(LD[[i]]), function(j)
          if (i %in% grep(reactive_values$clustergroups, names(LD), value = T) &
              j == input$pickerclusterfile) {
            LD[[i]][[j]][5] <<- input$pickerclusterfile
          } else{
            LD[[i]][[j]][5] <<- 0
          }))
      LIST_DATA$gene_info <- LD
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
      reactive_values$Apply_Cluster_Math <- ApplyMath(
        LIST_DATA,
        input$myMath,
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
      )
                   })
      if (!is.null(reactive_values$Apply_Cluster_Math)) {
        reactive_values$Plot_Cluster_Options <-
          MakePlotOptionFrame(LIST_DATA)
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    }
  })
  
  observeEvent(input$cluster3table_rows_all, ignoreInit = TRUE, {
    newname <-
      paste0(reactive_values$clustergroups,
             "3\nn = ",
             length(input$cluster3table_rows_all))
    oldname <-
      grep(paste0(reactive_values$clustergroups, "3\nn ="),
           names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("cluster3 filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster3table_rows_all])
      LD <- LIST_DATA$gene_info
      sapply(names(LD), function(i)
        sapply(names(LD[[i]]), function(j)
          if (i %in% grep(reactive_values$clustergroups, names(LD), value = T) &
              j == input$pickerclusterfile) {
            LD[[i]][[j]][5] <<- input$pickerclusterfile
          } else{
            LD[[i]][[j]][5] <<- 0
          }))
      LIST_DATA$gene_info <- LD
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
      reactive_values$Apply_Cluster_Math <- ApplyMath(
        LIST_DATA,
        input$myMath,
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
      )
                   })
      if (!is.null(reactive_values$Apply_Cluster_Math)) {
        reactive_values$Plot_Cluster_Options <-
          MakePlotOptionFrame(LIST_DATA)
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    }
  })
  
  observeEvent(input$cluster4table_rows_all, ignoreInit = TRUE, {
    newname <-
      paste0(reactive_values$clustergroups,
             "4\nn = ",
             length(input$cluster4table_rows_all))
    oldname <-
      grep(paste0(reactive_values$clustergroups, "4\nn ="),
           names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_file)[oldname]) {
      print("cluster4 filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cluster4table_rows_all])
      LD <- LIST_DATA$gene_info
      
      sapply(names(LD), function(i)
        sapply(names(LD[[i]]), function(j)
          if (i %in% grep(reactive_values$clustergroups, names(LD), value = T) &
              j == input$pickerclusterfile) {
            LD[[i]][[j]][5] <<- input$pickerclusterfile
          } else{
            LD[[i]][[j]][5] <<- 0
          }))
      LIST_DATA$gene_info <- LD
      withProgress(message = 'Calculation in progress',
                   detail = 'This may take a while...',
                   value = 0,
                   {
      reactive_values$Apply_Cluster_Math <- ApplyMath(
        LIST_DATA,
        input$myMath,
        input$radioplotnrom,
        as.numeric(input$sliderplotBinNorm)
      )
                   })
      if (!is.null(reactive_values$Apply_Cluster_Math)) {
        reactive_values$Plot_Cluster_Options <-
          MakePlotOptionFrame(LIST_DATA)
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
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
    if (n_distinct(LIST_DATA$gene_file[[input$selectclusterfile]]$use) < 4) {
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
    if (n_distinct(LIST_DATA$gene_file[[input$selectclusterfile]]$use) < 4) {
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
  observeEvent(c(input$selectclusternumber, reactive_values$clustergroups),
               ignoreInit = TRUE,
               {
                 print("cluster tool number")
                 if (is.null(reactive_values$clustergroups)) {
                   return()
                 }
                 hide('actionclusterdatatable')
                 hide('actionclusterplot')
                 hide('plotcluster')
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
                   show('actionclusterdatatable')
                   show('actionclusterplot')
                   glo <- input$selectgenelistoptions
                   if (!glo %in% names(LIST_DATA$gene_file)) {
                     glo <- names(LIST_DATA$gene_file)[1]
                   }
                   updateSelectInput(
                     session,
                     "selectgenelistoptions",
                     choices = names(LIST_DATA$gene_info),
                     selected = glo
                   )
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
                 } else {
                   return()
                 }
               })
  
  # Create and Show cluster data table ----
  observeEvent(input$actionclusterdatatable, ignoreInit = TRUE, {
    updateTabItems(session, "clustertooltab", "Cluster 1")
    newnames <- gsub("(.{20})", "\\1... ", input$pickerclusterfile)
    if (any(grep(
      paste0(reactive_values$clustergroups, "1\nn ="),
      names(LIST_DATA$gene_info)
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
                             names(LIST_DATA$gene_info)
                           )]]$full,
                           rownames = FALSE,
                           colnames = strtrim(newnames, 24),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "1\nn ="),
                             names(LIST_DATA$gene_info)
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
      show("cluster1table")
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
    
    if (any(grep(
      paste0(reactive_values$clustergroups, "2\nn ="),
      names(LIST_DATA$gene_info)
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
                             names(LIST_DATA$gene_info)
                           )]]$full,
                           rownames = FALSE,
                           colnames = strtrim(newnames, 24),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "2\nn ="),
                             names(LIST_DATA$gene_info)
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
      show("cluster2table")
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
    
    if (any(grep(
      paste0(reactive_values$clustergroups, "3\nn ="),
      names(LIST_DATA$gene_info)
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
                             names(LIST_DATA$gene_info)
                           )]]$full,
                           rownames = FALSE,
                           colnames = strtrim(newnames, 24),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "3\nn ="),
                             names(LIST_DATA$gene_info)
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
      show("cluster3table")
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
    
    if (any(grep(
      paste0(reactive_values$clustergroups, "4\nn ="),
      names(LIST_DATA$gene_info)
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
                             names(LIST_DATA$gene_info)
                           )]]$full,
                           rownames = FALSE,
                           colnames = strtrim(newnames, 24),
                           class = 'cell-border stripe compact',
                           filter = 'top',
                           caption = LIST_DATA$gene_file[[grep(
                             paste0(reactive_values$clustergroups, "4\nn ="),
                             names(LIST_DATA$gene_info)
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
      show("cluster4table")
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
    hide('actionclusterdatatable')
  })
  
  # creat and show cluster plot ----
  observeEvent(input$actionclusterplot, {
    print("cluster plot button")
    show('plotcluster')
    LD <- LIST_DATA$gene_info
    sapply(names(LIST_DATA$gene_info), function(i)
      sapply(names(LIST_DATA$gene_info[[i]]), function(j)
        if (i %in% grep(reactive_values$clustergroups,
                        names(LIST_DATA$gene_info),
                        value = T) & j == input$pickerclusterfile) {
          LIST_DATA$gene_info[[i]][[j]][5] <<- input$pickerclusterfile
        } else{
          LIST_DATA$gene_info[[i]][[j]][5] <<- 0
        }))
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...',
                 value = 0,
                 {
    reactive_values$Apply_Cluster_Math <- ApplyMath(
      LIST_DATA,
      input$myMath,
      input$radioplotnrom,
      as.numeric(input$sliderplotBinNorm)
    )
                 })
    if (!is.null(reactive_values$Apply_Cluster_Math)) {
      reactive_values$Plot_Cluster_Options <-
        MakePlotOptionFrame(LIST_DATA)
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
  observeEvent(c(input$pickercdffile1),
               ignoreInit = TRUE,
               ignoreNULL = FALSE,
               {
                 if (is.null(input$pickercdffile1)) {
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
    if (is_empty(oldname)) {
      return()
    }
    gene_count <-
      n_distinct(LIST_DATA$gene_file[[oldname]]$full$gene)
    num <- c(
      ceiling(gene_count * input$slidercdfper[1] / 100),
      ceiling(gene_count * input$slidercdfper[2] / 100)
    )
    gene_list <-
      group_by(LIST_DATA$gene_file[[oldname]]$full, gene) %>%
      filter(all(between(bin, num[1], num[2]))) %>%
      distinct(gene) %>%
      ungroup()
    newname <- paste("CDF\nn =", n_distinct(gene_list$gene))
    if (newname != names(LIST_DATA$gene_info)[oldname]) {
      hide('cdftable')
      show('actioncdfdatatable')
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<- gene_list
      df_options <-
        semi_join(
          bind_rows(LIST_DATA$gene_info[[newname]]),
          distinct(LIST_DATA$gene_file[[newname]]$full, set),
          by = "set"
        ) %>%
        mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{17})", "\\1\n", set),
          sep = '\n'
        ))
      df <- inner_join(LIST_DATA$gene_file[[newname]]$full,
                       LIST_DATA$gene_file[[newname]]$use,
                       by = "gene") %>%
        mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{17})", "\\1\n", set),
          sep = '\n'
        ))
      
      use_header <- pull(distinct(df_options, myheader))
      if (n_groups(group_by(df_options, set)) == 2 &
          n_distinct(df$gene) > 1) {
        tt_name <- pull(distinct(df_options, set))
        tt <-
          suppressWarnings(ks.test(pull(filter(
            df, set == tt_name[1]
          ), value),
          pull(filter(
            df, set == tt_name[2]
          ), value)))
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
  observeEvent(input$actioncdfdatatable, ignoreInit = TRUE, {
    if (any(grep("CDF\nn", names(LIST_DATA$gene_info)) > 0)) {
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
    hide('actioncdfdatatable')
    show('cdftable')
  })
  
  
  # cdf tool gene lists $use ----
  observeEvent(input$cdftable_rows_all, ignoreInit = TRUE, {
    newname <- paste("CDF\nn =", length(input$cdftable_rows_all))
    oldname <- grep("CDF\nn =", names(LIST_DATA$gene_info))
    if (newname != names(LIST_DATA$gene_info)[oldname]) {
      print("cdf filter $use")
      LIST_DATA$STATE[2] <<- newname
      names(LIST_DATA$gene_file)[oldname] <<- LIST_DATA$STATE[2]
      names(LIST_DATA$gene_info)[oldname] <<- LIST_DATA$STATE[2]
      LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$use <<-
        tibble(gene = LIST_DATA$gene_file[[LIST_DATA$STATE[2]]]$full$gene[input$cdftable_rows_all])
      df_options <-
        semi_join(
          bind_rows(LIST_DATA$gene_info[[newname]]),
          distinct(LIST_DATA$gene_file[[newname]]$full, set),
          by = "set"
        ) %>%
        mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{17})", "\\1\n", set),
          sep = '\n'
        ))
      df <- inner_join(LIST_DATA$gene_file[[newname]]$full,
                       LIST_DATA$gene_file[[newname]]$use,
                       by = "gene") %>%
        mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{17})", "\\1\n", set),
          sep = '\n'
        ))
      
      use_header <- pull(distinct(df_options, myheader))
      if (n_groups(group_by(df_options, set)) == 2 &
          n_distinct(df$gene) > 1) {
        tt_name <- pull(distinct(df_options, set))
        tt <-
          suppressWarnings(ks.test(pull(filter(
            df, set == tt_name[1]
          ), value),
          pull(filter(
            df, set == tt_name[2]
          ), value)))
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
      ol1 <- input$selectcdffile1
      if (!ol1 %in% names(LIST_DATA$gene_file)) {
        ol1 <- newname
        reactive_values$pickerfile_controler <- input$pickercdffile1
      } else {
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
      newname <-
        grep("CDF\nn =", names(LIST_DATA$gene_info), value = TRUE)
      df_options <-
        semi_join(
          bind_rows(LIST_DATA$gene_info[[newname]]),
          distinct(LIST_DATA$gene_file[[newname]]$full, set),
          by = "set"
        ) %>%
        mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{17})", "\\1\n", set),
          sep = '\n'
        ))
      df <- inner_join(LIST_DATA$gene_file[[newname]]$full,
                       LIST_DATA$gene_file[[newname]]$use,
                       by = "gene") %>%
        mutate(set = paste(
          sub("\n", " ", newname),
          gsub("(.{17})", "\\1\n", set),
          sep = '\n'
        ))
      
      use_header <- pull(distinct(df_options, myheader))
      if (n_groups(group_by(df_options, set)) == 2 &
          n_distinct(df$gene) > 1) {
        tt_name <- pull(distinct(df_options, set))
        tt <-
          suppressWarnings(ks.test(pull(filter(
            df, set == tt_name[1]
          ), value),
          pull(filter(
            df, set == tt_name[2]
          ), value)))
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
      if (!glo %in% names(LIST_DATA$gene_file)) {
        glo <- names(LIST_DATA$gene_file)[1]
      }
      updateSelectInput(
        session,
        "selectgenelistoptions",
        choices = names(LIST_DATA$gene_info),
        selected = glo
      )
      ol <- input$selectcdffile1
      if (!ol %in% names(LIST_DATA$gene_file)) {
        ol <- grep("CDF\nn", names(LIST_DATA$gene_file), value = TRUE)
        reactive_values$pickerfile_controler <- input$pickercdffile1
      } else {
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
                    width = 8,
                    status = "primary",
                    solidHeader = T,
                    textInput("textnickname", "Update Nickname"),
                    actionButton("actionoptions", "Set Nickname"),
                    fluidRow(
                      column(4,
                    numericInput("normfactor", "Set norm factor, score/rpm", value = 1)
                    ),
                    column(4, style = "padding-top:10%;",
                    actionButton("actionnormfactor", "Apply norm factor")
                    )
                    ),
                    helpText("Need to update Nickname and/or nrom factor"),
                    checkboxInput("checkboxremovefile",
                                  "remove all files and restart", value = FALSE),
                    actionButton("actionremovefile", "Remove File(s)")
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
                    selectInput("kbrewer",
                                "color brewer",
                                c(choices = "select", kBrewerList))
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
                pickerInput(
                  "pickerdenominator",
                  label = "denominator",
                  width = "90%",
                  choices = "Load data file",
                  multiple = F,
                  options = list(title = "Select second file")
                )
              )),
          actionButton("actionnorm", label = "create norm file"),
          checkboxInput("checkboxnormmean", label = "gene by gene", value = TRUE),
          checkboxInput("checkboxnormzero", label = "denom 0 -> min/2", value = TRUE),
          helpText(
            "if 0's are not converted genes containing will be removed from all gene lists"
          )
          
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
                      min = -20,
                      max = 120,
                      post = "%",
                      value = c(0, 100)
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
                    awesomeRadio(
                      "radioplotnrom",
                      label = "Set Y Normalization",
                      choices = c("none", "relative frequency", "rel gene frequency"),
                      selected = "none"
                    ),
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
                  actionButton("actiongenelists", "Compare Gene lists"),
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
                    tabPanel(
                      "Intersected Gene Lists",
                      DT::dataTableOutput('genelists1table')
                    ),
                    tabPanel(
                      "Inclusive Gene Lists",
                      DT::dataTableOutput('genelists2table')
                    ),
                    tabPanel(
                      "Exclusive Gene Lists",
                      DT::dataTableOutput('genelists3table')
                    )
                  )
                )
              )),
      # main sort tab ----
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
                      6,
                      sliderInput(
                        "slidersortpercent",
                        label = "% select:",
                        post = "%",
                        min = 1,
                        max = 100,
                        value = 75
                      )
                    ),
                    column(
                      6,
                      sliderInput(
                        "slidersortbinrange",
                        label = "Select Bin Range:",
                        min = 0,
                        max = 80,
                        value = c(0, 80)
                      )
                    )
                  ),
                  fluidRow(
                    column(
                      3,
                      pickerInput(
                        "selectsorttop",
                        "Sort Options",
                        choices = c("Top%", "Middle%", "Bottom%"),
                        selected = "Middle%"
                      )
                    )
                  ),
                  actionButton("actionsorttool", "Sort")
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
              )),
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
                  helpText(
                    "if 0's are not converted gene's containing region with sum(0) will be removed from results"
                  )
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
                  fluidRow(column(
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
                  )),
                  actionButton("actionclustertool", "Get clusters"),
                  actionButton("actiongroupstool", "Get groups")
                  
                  
                ),
                box(
                  title = "Cluster Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  width = 12,
                  collapsible = TRUE,
                  collapsed = TRUE,
                  actionButton("actionclusterplot", "plot"),
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
              div(
                id = "enablemaincdf",
                box(
                  title = "CDF tool",
                  status = "primary",
                  solidHeader = T,
                  width = 12,
                  fluidRow(column(
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
                  )),
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
                  helpText(
                    "if 0's are not converted gene's containing region with sum(0) will be removed from results"
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
                    width = 12,
                    actionButton("actioncdfdatatable", "Show gene list(s)"),
                    DT::dataTableOutput('cdftable')
                  )
                )
              ))
    )
  )
      )

# exicute ----
shinyApp(ui = ui, server = server)