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
    "shinythemes",
    "ggplot2",
    "tcltk2",
    "dplyr",
    "tidyr",
    "readr",
    "fastcluster",
    "RColorBrewer"
  )
))

# TODO do I need all theses lists
LIST_DATA <- list(
  table_file = list(),
  # [[]] gene X1 X2 ...
  gene_file = list(),
  # holds $common genes from files and $gene file(s)
  gene_info = list(),
  # for holding gene file info in a list of lists, a set for $common and each $gene file(s) [c("dot", "line", "color", plot?, NickName,nrom)]
  clust = list()
)      # Cluster holder

# values for comboboxs ----

kplotBinRange <- c(0, 0, 0, 0)

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
          set = gsub("(.{17})", "\\1\n", legend_nickname),
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
      kplotBinRange <<- c(1, num_bins, 1, num_bins)
    }
    color_safe <-
      (length(LIST_DATA$table_file) + 1) %% length(kListColorSet)
    if (color_safe == 0) {
      color_safe <- 1
    }
    color_select <- kListColorSet[color_safe]
    
    list_data$table_file[[legend_nickname]] <- tablefile
    list_data$gene_file$common$use <- gene_names
    list_data$gene_info$common[[legend_nickname]] <-
      c(
        kDotOptions[1],
        kLineOptions[1],
        color_select,
        legend_nickname,
        legend_nickname,
        "none"
      )
    
    # generate info for new file for loaded gene list(s)
    sapply(names(list_data$gene_file), function(g) {
      if (g != "common" & length(list_data$gene_file[[g]]$full) > 0) {
        enesg <- c(gene_names, list_data$gene_file[[g]]$full)
        enesg <- enesg[duplicated(enesg)]
        list_data$gene_file[[g]]$use <<- enesg
        list_data$gene_info[[g]][[legend_nickname]] <<-
          c(
            kDotOptions[1],
            kLineOptions[1],
            color_select,
            legend_nickname,
            legend_nickname,
            "none"
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
      list_data[[gene_set]][[i]][4] <- 0
    } else {
      list_data[[gene_set]][[i]][4] <- i
    }
  }
  list_data
}

# Makes data frame and gathers plot settings for plotting active samples
MakeDataFrame <-
  function(sel_list = NULL,
           table_file = LIST_DATA$table_file,
           gene_file = LIST_DATA$gene_file,
           gene_info = LIST_DATA$gene_info,
           nickname = NULL) {
    if (is.null(names(table_file))) {
      return ()
    } else {
      use_col <- NULL
      use_dot <- NULL
      use_line <- NULL
      use_size <- NULL
      use_nickname <- NULL
      use_x_label <- NULL
      legend_space <- 1
      list_data_frame <- NULL
      
      for (i in names(gene_file)) {
        # checks to see if at least one file in list is acitve
        if (sum(sapply(gene_info[[i]], "[[", 4) == 0) != 0) {
          next ()
        } else {
          if (!is.null(sel_list)) {
            enesg <- c(sel_list, gene_file[[i]]$use)
            enesg <- data_frame(gene = enesg[duplicated(enesg)])
            use_x_label <- paste(use_x_label, paste(i, "n = ",
                                                    length(enesg[[1]])), sep = '  ')
            if (length(enesg[[1]]) == 0) {
              break()
            }
          } else {
            enesg <- data_frame(gene = gene_file[[i]]$use)
            use_x_label <- paste(use_x_label, paste(i, "n = ",
                                                    length(enesg[[1]])), sep = '  ')
          }
          tf <-
            c(sapply(
              gene_info[[i]], "[[", 4
            ) != 0)
          list_data_frame[[i]] <-
            bind_rows(table_file[tf]) %>%
            semi_join(., enesg, by = "gene") %>%
            mutate(., set = paste(gsub("(.{17})", "\\1\n", i), set, sep = '-\n'))
          
          dots <-
            match(sapply(gene_info[[i]][tf], "[[", 1), kDotOptions)
          use_dot <- c(use_dot,
                       if_else(dots == 1, 0, dots + 13))
          use_size <- c(use_size, if_else(use_dot == 0, 0.01, 4.5))
          my_lines <-
            match(sapply(gene_info[[i]][tf], "[[", 2), kLineOptions)
          use_line <- c(use_line,
                        if_else(my_lines > 6, 0, as.double(my_lines)))
          use_col <- c(use_col, sapply(gene_info[[i]][tf], "[[", 3))
          current_nickname <- select(list_data_frame[[i]], set) %>%
            distinct()
          use_nickname <-
            c(use_nickname, current_nickname$set)
          legend_space <-
            max(legend_space, (lengths(
              strsplit(current_nickname$set, "\n")
            ) - 0.75))
        }
      }
      if (!is.null(nickname)) {
        use_x_label <- paste(nickname, use_x_label)
      }
      if (!is.null(names(list_data_frame))) {
        print("tt")
        ApplyMath(
          list_data_frame,
          use_col,
          use_dot,
          use_line,
          use_size,
          use_nickname,
          use_x_label,
          legend_space
        )
      } 
    }
  }

# Applys math to long list
ApplyMath <-
  function(list_data_frame,
           use_col,
           use_dot,
           use_line,
           use_size,
           use_nickname,
           use_x_label,
           legend_space) {
    # math set and y label
    
    use_math <- "mean" # fix
    
    use_y_label <- paste(use_math, "of bin counts")
    # set normilization to relative frequency or bin number or 1 for none,
    # and update y label
    r_checkbox_relative_frequency <- 0 # fix
    r_checkbox_gene_relative_frequency <- 0 # fix
    norm_bin <- 0 # fix
    if (r_checkbox_gene_relative_frequency == 1) {
      if (r_checkbox_relative_frequency == 1) {
        # tktoggle(checkbox_relative_frequency) # fix
      }
      list_long_data_frame <- bind_rows(list_data_frame) %>%
        group_by(set, gene) %>%
        mutate(score = score / sum(score, na.rm = TRUE)) %>%
        ungroup() %>%
        group_by(set, bin) %>%
        summarise(value = get(use_math)(score, na.rm = T)) %>%
        ungroup()
      
      use_y_label <- paste("RF per gene :", use_y_label)
    } else {
      list_long_data_frame <- bind_rows(list_data_frame) %>%
        group_by(set, bin) %>%
        summarise(value = get(use_math)(score, na.rm = T)) %>%
        ungroup()
    }
    
    if (r_checkbox_relative_frequency == 1 && norm_bin == 0) {
      use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                           "bins : RF")
      list_long_data_frame <-
        group_by(list_long_data_frame, set) %>%
        mutate(value = value / sum(value)) %>%
        ungroup()
    } else if (norm_bin > 0) {
      if (r_checkbox_relative_frequency == 1) {
        # tktoggle(checkbox_relative_frequency) # fix
      }
      if (r_checkbox_gene_relative_frequency == 1) {
        use_y_label <- paste(use_y_label, " : Norm bin ", norm_bin)
      } else {
        use_y_label <- paste(strsplit(use_y_label, split = " ")[[1]][1],
                             "bins : Normalize to bin ",
                             norm_bin)
      }
      list_long_data_frame <-
        group_by(list_long_data_frame, set) %>%
        mutate(value = value / nth(value, norm_bin)) %>%
        ungroup()
    }
    
    use_plot_limits <-
      c(kplotBinRange[3], kplotBinRange[4])
    
    # update y label if log2
    if (0 == 1) { # fix
      use_y_label <- paste("log2(", use_y_label, ")", sep = "")
      if (summarise(list_long_data_frame, sum(value, na.rm = TRUE)) != 0) {
        list_long_data_frame$value[list_long_data_frame$value == 0] <-
          min(list_long_data_frame$value[list_long_data_frame$value > 0]) / 2
      } else {
        list_long_data_frame$value <- 1
      }
      use_y_limits <- group_by(list_long_data_frame, set) %>%
        filter(bin %in% use_plot_limits[1]:use_plot_limits[2]) %>%
        ungroup() %>%
        summarise(min(log2(value), na.rm = T), max(log2(value), na.rm =
                                                     T))
    } else {
      use_y_limits <- group_by(list_long_data_frame, set) %>%
        filter(bin %in% use_plot_limits[1]:use_plot_limits[2]) %>%
        ungroup() %>%
        summarise(min(value, na.rm = T), max(value, na.rm = T))
    }
    if (0 == 1 & # fix
        !suppressWarnings(is.na(as.numeric(NA))) & # fix
        !suppressWarnings(is.na(as.numeric(NA)))) { # fix
      use_y_limits <-
        c(as.numeric(0), as.numeric(0)) # fix
    } else { # fix
      # tkdelete(entry_Y_axis_min, 0, "end")
      # tkdelete(entry_Y_axis_max, 0, "end")
      # tkinsert(entry_Y_axis_min, 0, round(unlist(use_y_limits)[1], digits = 3))
      # tkinsert(entry_Y_axis_max, 0,  round(unlist(use_y_limits)[2], digits = 3))
    }
    use_pos_plot_ticks <- c(1:16) # fix 
    # Destring(unlist(strsplit(tclvalue(tcl_pos_plot_ticks), " ")))
    use_label_plot_ticks <- c(1:16) # fix
      # unlist(strsplit(tclvalue(tcl_label_plot_ticks),
      #                 " "))
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
    use_virtical_line_type <-
      c(1,1,1,1
        # as.character(tclvalue(tcl_tss_tts_line_type)),
        # as.character(tclvalue(tcl_tss_tts_line_type)),
        # as.character(tclvalue(tcl_body_line_type)),
        # as.character(tclvalue(tcl_body_line_type))
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
    
    # TODO need controls?
    virtical_line_data_frame <- data.frame(
      use_virtical_line,
      use_virtical_line_type,
      use_virtical_line_color,
      stringsAsFactors = FALSE
    )
    names(use_col) <- use_nickname
    names(use_dot) <- use_nickname
    names(use_line) <- use_nickname
    # if(tclvalue(tcl_smooth) == 1){
    #   use_y_label <- paste0("smoothed(", use_y_label, ")")
    # }
    
    GGplotF(
      list_long_data_frame,
      use_col,
      use_dot,
      use_line,
      use_size,
      use_y_label,
      use_x_label,
      use_plot_breaks,
      virtical_line_data_frame,
      use_plot_breaks_labels,
      use_plot_limits,
      use_y_limits,
      legend_space
    )
  }

# main ggplot function
GGplotF <-
  function(list_long_data_frame,
           use_col,
           use_dot,
           use_line,
           use_size,
           use_y_label,
           use_x_label,
           use_plot_breaks,
           virtical_line_data_frame,
           use_plot_breaks_labels,
           use_plot_limits,
           use_y_limits,
           legend_space) {
    if (0 == 1) { # fix
      gp <- ggplot(
        list_long_data_frame,
        aes(
          x = as.numeric(bin),
          y = log2(value),
          color = set,
          shape = set,
          linetype = set,
          size = set
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
            linetype = set,
            size = set
          )
        )
    }
    if (0 == 1) { #fix
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
      scale_size_manual(name = "Sample", values = use_size) +
      scale_color_manual(name = "Sample", values = use_col) +
      scale_shape_manual(name = "Sample", values = use_dot) +
      scale_linetype_manual(name = "Sample", values = use_line) +
      xlab(use_x_label) + ylab(use_y_label) +  # Set axis labels
      scale_x_continuous(breaks = use_plot_breaks,
                         labels = use_plot_breaks_labels) +
      
      geom_vline(
        data = virtical_line_data_frame,
        aes(xintercept = use_virtical_line),
        size = 2,
        linetype = virtical_line_data_frame$use_virtical_line_type,
        color = virtical_line_data_frame$use_virtical_line_color
      ) +
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
      ) +
      coord_cartesian(xlim = use_plot_limits, ylim = unlist(use_y_limits))
    suppressMessages(print(gp))
  }