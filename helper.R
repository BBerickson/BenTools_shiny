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
    "ggplot2",
    "dplyr",
    "tidyr",
    "readr",
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
  STATE = c(0, "common")
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
kListColorSet <- brewer.pal(8, kBrewerList[3])

# math options avalible
kMathOptions <- c("mean", "sum", "median", "var")

# functions ----
RgbToHex <- function(my_hex = NULL, my_rgb = NULL){
  if(!is.null(my_hex)){
    return(paste(col2rgb(c(my_hex)),collapse = ","))
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
# reads in file, tests, fills out info and returns list_data
LoadTableFile <- function(file_path, file_name, list_data) {
  # load remote files TODO
  # if(length(file_path) == 1 && length(grep(".url",file_path))==1){
  #   file_path <- read_lines(file_path)
  # }
  file_count <- length(list_data$table_file)
  for (x in file_path) {
    legend_nickname <-
      strsplit(as.character(file_name), '.tab')[[1]][1]
    if (any(legend_nickname == names(list_data$table_file))) {
      # tkmessageBox(message = "This file has already been loaded", icon = "info") TODO
      next()
    } else {
      # TODO numbins test
      num_bins <-
        count_fields(x,
                     n_max = 1,
                     skip = 1,
                     tokenizer = tokenizer_tsv()) - 1
      
      tablefile <- suppressMessages(read_tsv(x,
                                             col_names = c("gene", 1:(num_bins)),
                                             skip = 1) %>%
                                      gather(., bin, score, 2:(num_bins + 1))) %>%
        mutate(
          set = legend_nickname,
          bin = as.numeric(bin),
          score = as.numeric(score)
        ) %>%
        na_if(Inf) %>%
        replace_na(list(score = 0))
    }
    zero_genes <-
      group_by(tablefile, gene) %>% summarise(test = sum(score, na.rm = T)) %>% filter(test !=
                                                                                         0)
    
    tablefile <- semi_join(tablefile, zero_genes, by = "gene")
    gene_names <- collapse(distinct(tablefile, gene))[[1]]
    if (file_count > 0) {
      gene_names <-
        c(list_data$gene_file$common$use, gene_names)
      gene_names <- gene_names[duplicated(gene_names)]
      if (length(gene_names) == 0) {
        # tkmessageBox("no genes in common") TODO
        break()
      }
    } else {
      list_data$x_plot_range <- c(1, num_bins)
    }
    color_safe <-
      (length(list_data$table_file) + 1) %% length(kListColorSet)
    if (color_safe == 0) {
      color_safe <- 1
    }
    color_select <- kListColorSet[color_safe]
    
    list_data$table_file[[legend_nickname]] <- tablefile
    list_data$gene_file$common$use <- gene_names
    list_data$gene_info$common[[legend_nickname]] <-
      # don't change the order of postions
      tibble(
        set = legend_nickname,
        mydot = kDotOptions[1],
        myline = kLineOptions[1],
        mycol = color_select,
        onoff = legend_nickname,
        rnorm = "none"
      )
    
    # generate info for new file for loaded gene list(s)
    sapply(names(list_data$gene_file), function(g) {
      if (g != "common" & length(list_data$gene_file[[g]]$full) > 0) {
        enesg <- c(gene_names, list_data$gene_file[[g]]$full)
        enesg <- enesg[duplicated(enesg)]
        list_data$gene_file[[g]]$use <<- enesg
        list_data$gene_info[[g]][[legend_nickname]] <<-
          tibble(
            set = legend_nickname,
            mydot = kDotOptions[1],
            myline = kLineOptions[1],
            mycol = color_select,
            onoff = 0,
            rnorm = "none"
          )
        
        # if 0 length message and remove file TODO
      }
    })
    file_count <- 1
  }
  list_data
}

# records check box on/off
CheckBoxOnOff <- function(gene_set, check_box, list_data) {
  for (i in names(list_data[[gene_set]])) {
    # make function
    if (!i %in% check_box) {
      list_data[[gene_set]][[i]]["onoff"] <- 0
    } else {
      list_data[[gene_set]][[i]]["onoff"] <- i
    }
  }
  list_data
}

# Applys math to on data in each gene list
ApplyMath <-
  function(list_data, 
           use_math, 
           r_checkbox_gene_relative_frequency,
           sel_list = NULL) {
    
    table_file = list_data$table_file
    gene_file = list_data$gene_file
    gene_info = list_data$gene_info
    list_data_frame <- NULL
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
            enesg <- data_frame(gene = gene_file[[i]]$use)
          }
          truefalse <-
            c(sapply(
              gene_info[[i]], "[[", 5
            ) != 0)
          list_data_frame[[i]] <-
            bind_rows(table_file[truefalse]) %>%
            semi_join(., enesg, by = "gene") %>%
            mutate(., set = paste(gsub("(.{17})", "\\1\n", i), gsub("(.{17})", "\\1\n", set), sep = '-\n'))
        }
      }
      
      if (is.null(names(list_data_frame))) {
        print("nothing to plot")
        return(NULL)
      }
    # applys math to pared down data file
    if (r_checkbox_gene_relative_frequency == 1) {
      list_long_data_frame <- bind_rows(list_data_frame) %>%
        group_by(set, gene) %>%
        mutate(score = score / sum(score, na.rm = TRUE)) %>%
        ungroup() %>%
        group_by(set, bin) %>%
        summarise(value = get(use_math)(score, na.rm = T)) %>%
        ungroup()
      
    } else {
      list_long_data_frame <- bind_rows(list_data_frame) %>%
        group_by(set, bin) %>%
        summarise(value = get(use_math)(score, na.rm = T)) %>%
        ungroup()
    }
    
    return(list_long_data_frame)
  }

# gather relavent plot option data
MakePlotOptionFrame <- function(list_data){
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
               set = paste(gsub("(.{17})", "\\1\n", i), gsub("(.{17})", "\\1\n", set), sep = '-\n'))
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
YAxisLable <- function(use_math = "mean", use_log2 = 0, norm_bin = 0){
  use_y_label <- paste(use_math, "of bin counts")
  if ("r_checkbox_gene_relative_frequency" == 1) {
    use_y_label <- paste("RF per gene :", use_y_label)
  } else if ("r_checkbox_relative_frequency" == 1) {
    use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                         "bins : RF")
  }
  if (norm_bin > 0) {
    if (r_checkbox_gene_relative_frequency == 1) {
      use_y_label <- paste(use_y_label, " : Norm bin ", norm_bin)
    } else {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : Normalize to bin ",
                           norm_bin)
    }
  }
  if (use_log2 == 1) { # fix
    use_y_label <- paste("log2(", use_y_label, ")", sep = "")
  }
  # if(tclvalue(tcl_smooth) == 1){
  #   use_y_label <- paste0("smoothed(", use_y_label, ")")
  # }
  use_y_label
}

# Sets plot lines and lables fix
LinesLablesList <- function(use_pos_plot_ticks,
                        use_label_plot_ticks,
                        use_virtical_line_type = c(1,1,1,1)
                                               ){
  
  
  if (length(use_label_plot_ticks) != length(use_pos_plot_ticks)) {
    # fix tkmessageBox(message = "The number of Positions and labels are unequal", icon = "error")
    return ()
  }
  
  use_plot_breaks <- c(
    # Destring(tclvalue(tcl_pos_one_line)),
    # Destring(tclvalue(tcl_pos_two_line)),
    # Destring(tclvalue(tcl_pos_three_line)),
    # Destring(tclvalue(tcl_pos_four_line)),
    use_pos_plot_ticks
  )
  use_plot_breaks[use_plot_breaks == 0] <- NA
  
  use_plot_breaks_labels <- c(
    # tclvalue(tcl_one_tss_tts_option),
    # tclvalue(tcl_two_tss_tts_option),
    # tclvalue(tcl_three_tss_tts_option),
    # tclvalue(tcl_four_tss_tts_option),
    use_label_plot_ticks
  )
  
  use_virtical_line_color <- c("green", "red", "black", "black")
  use_plot_breaks_labels <-
    use_plot_breaks_labels[!is.na(use_plot_breaks)]
  use_virtical_line_type <-
    use_virtical_line_type[!is.na(use_plot_breaks[1:4])]
  use_virtical_line_color <-
    use_virtical_line_color[!is.na(use_plot_breaks[1:4])]
  use_virtical_line <-
    use_plot_breaks[1:4][!is.na(use_plot_breaks[1:4])]
  use_plot_breaks <- use_plot_breaks[!is.na(use_plot_breaks)]
  
  virtical_line_data_frame <- data.frame(
    use_virtical_line,
    use_virtical_line_type,
    use_virtical_line_color,
    use_plot_breaks,
    stringsAsFactors = FALSE
  )
}

# main ggplot function
GGplotLineDot <-
  function(list_long_data_frame, xBinRange, plot_options, use_smooth = 0, legend_space = 1) {
    print("ggplot")
    gp <-
        ggplot(
          list_long_data_frame,
          aes(
            x = as.numeric(bin),
            y = value,
            color = set,
            shape = set,
            #size = set
            linetype = set
            
          )
        )
    if (use_smooth == 1) { #fix
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
      # scale_size_manual(name = "Sample", values = use_size) +
      scale_color_manual(name = "Sample", values = plot_options$mycol) + # do I need to add names to keep things straight? 
      scale_shape_manual(name = "Sample", values = plot_options$mydot) +
      scale_linetype_manual(name = "Sample", values = plot_options$myline) +
      # xlab(use_x_label) + ylab(use_y_label) +  # Set axis labels
      # scale_x_continuous(breaks = use_plot_breaks,
      #                    labels = use_plot_breaks_labels) +
      # 
      # geom_vline(
      #   data = virtical_line_data_frame,
      #   aes(xintercept = use_virtical_line),
      #   size = 2,
      #   linetype = virtical_line_data_frame$use_virtical_line_type,
      #   color = virtical_line_data_frame$use_virtical_line_color
      # ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank()) +
      theme(axis.title.y = element_text(size =  15)) +
      theme(axis.title.x = element_text(size =  10, vjust = .5)) +
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
        legend.text = element_text(size = 9)
      )  +
      coord_cartesian(xlim = xBinRange)
      #coord_cartesian(xlim = use_x_limits, ylim = unlist(use_y_limits))
    #suppressMessages(print(gp))
    return(gp)
  }

