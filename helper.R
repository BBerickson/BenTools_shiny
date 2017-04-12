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
  STATE = c(0, "common") # flow control, gene list flow control
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
      num_bins <-
        count_fields(file_path[x],
                     n_max = 1,
                     skip = 1,
                     tokenizer = tokenizer_tsv()) - 1
      # checks if last row is an empty and fixes
      tablefile <- suppressMessages(read_tsv(file_path[x],
                                             col_names = c("gene", 1:(num_bins)),
                                             skip = 1, n_max = 5))
      if(sum(is.na(tablefile[,num_bins+1])) == 5){
        test <- TRUE
        num_bins2 <- num_bins - 1
      } else {
        test <- FALSE
        num_bins2 <- num_bins
      }
      
      if(file_count > 0 & num_bins2 != list_data$x_plot_range[2]){
        showModal(modalDialog(
          title = "Information message",
          "Can't load file, different number of bins", size = "s",
          easyClose = TRUE
        ))
        next()
      }
      
      
      tablefile <- suppressMessages(read_tsv(file_path[x],
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
      
      if(test){
        tablefile <- tablefile %>% filter(., bin != num_bins)
      }
      
      num_bins <- num_bins2

    zero_genes <-
      group_by(tablefile, gene) %>% summarise(test = sum(score, na.rm = T)) %>% filter(test !=
                                                                                         0)
    
    tablefile <- semi_join(tablefile, zero_genes, by = "gene")
    gene_names <- collapse(distinct(tablefile, gene))[[1]]
    if (file_count > 0) {
      gene_names <-
        c(list_data$gene_file[[1]]$use, gene_names)
      gene_names <- gene_names[duplicated(gene_names)]
      if (length(gene_names) == 0) {
        showModal(modalDialog(
          title = "Information message",
          " No genes in common ", size = "s",
          easyClose = TRUE
        ))
        break()
      }
    } else {
      list_data$x_plot_range <- c(1, num_bins)

    }
    my_name <- paste("common\nn =", length(gene_names))
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
    
    # generate info for new file for loaded gene list(s)
    sapply(names(list_data$gene_file), function(g) {
      if (length(list_data$gene_file[[g]]$full) > 0) {
        enesg <- c(list_data$gene_file[[g]]$full, gene_names)
        enesg <- enesg[duplicated(enesg)]
        if(length(enesg) < 1){
          showModal(modalDialog(
            title = "Information message",
            " No genes in common, need to remove gene file", size = "s",
            easyClose = TRUE
          ))
        }
        my_name_g <- paste(strsplit(g, "\nn =")[[1]][1], "\nn =", length(enesg))
        names(list_data$gene_file)[which(names(list_data$gene_file) == g)] <<- my_name_g
        names(list_data$gene_info)[which(names(list_data$gene_info) == g)] <<- my_name_g
        list_data$gene_file[[my_name_g]]$use <<- enesg
        list_data$gene_info[[my_name_g]][[legend_nickname]] <<-
          tibble(
            set = legend_nickname,
            mydot = kDotOptions[1],
            myline = kLineOptions[1],
            mycol = color_select,
            onoff = 0,
            rnorm = "1"
          )
      }
    })
    file_count <- 1
  }
  list_data
}

# reads in gene list files
LoadGeneFile <- function(file_path, file_name, list_data) {
  
    genefile <-
      suppressMessages(read_tsv(
        file_path,
        col_names = "gene",
        comment = "#",
        cols(gene = col_character())
      ))
    
    enesg <- c(collapse(distinct(genefile, gene))[[1]],
               list_data$gene_file[[1]]$use)
    enesg <- enesg[duplicated(enesg)]
    if (length(enesg) == 0) {
      
      showModal(modalDialog(
        title = "Information message",
        " No genes in common, might need to reformat gene name style", size = "s",
        easyClose = TRUE
      ))
      return()
      
    }
    
    legend_nickname <-
      paste(strsplit(as.character(file_name), '.txt')[[1]][1], "\nn =", length(enesg))
    
    list_data$gene_file[[legend_nickname]]$full <-
      collapse(distinct(genefile, gene))[[1]]
    list_data$gene_file[[legend_nickname]]$use <- enesg
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
          mycol = list_data$gene_info[[1]][[i]]$mycol,
          onoff = 0,
          rnorm = "1"
        ))
    list_data$STATE[2] <- legend_nickname
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
            enesg <- data_frame(gene = gene_file[[i]]$use)
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
    use_plot_breaks <- c(1,rev(seq(totbins, by = -everybin, length.out = (totbins/everybin))))
    use_plot_breaks[near(use_plot_breaks, tssbin, tol = everybin -1)] <- tssbin + .5
    use_plot_breaks_labels <- rev(seq((totbins-tssbin)*binbp, by = -everybp, length.out = (totbins/everybin) + 1))
    use_plot_breaks_labels[near(use_plot_breaks, tssbin, tol = everybin -1)] <- "TSS"
    use_virtical_line <- c(tssbin, NA, NA, NA) + .5
  } else if(mytype == "3'"){
    use_plot_breaks <- c(1,rev(seq(totbins, by = -everybin, length.out = (totbins/everybin))))
    use_plot_breaks[near(use_plot_breaks, tesbin, tol = everybin -1)] <- tesbin + .5
    use_plot_breaks_labels <- abs(seq(-(tesbin-totbins)*binbp, by = binbp * everybin, length.out = (totbins/everybin) + 1))
    use_plot_breaks_labels[near(use_plot_breaks, tesbin, tol = everybin -1)] <- "TES"
    use_virtical_line <- c(NA, tesbin, NA, NA) + .5
  } else if(mytype == "none"){
    use_plot_breaks <- c(1,rev(seq(totbins, by = -everybin, length.out = (totbins/everybin))))
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
  if(mytype == "543 bins 20,20,40"){
    tt <- c(20,40,15,45,100,5)
  } else if (mytype == "543 bins 10,10,10"){ 
    tt <- c(10,20,5,25,100,5)
  } else if(mytype == "5' 1k 1k 80bins"){
    tt <- c(0,0,40,0,25,20)
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
    unlist(.,use.names=FALSE) #%>% round(.,4)
  tt <- round(c(tt, tt[1]-tt[1]*.1, tt[2]+tt[2]*.1),4) 
}

# main ggplot function
GGplotLineDot <-
  function(list_long_data_frame, xBinRange, plot_options, yBinRange, line_list, use_smooth, use_y_label, legend_space = 1) {
    print("ggplot")
    
    use_col <- plot_options$mycol
    use_dot <- plot_options$mydot
    use_line <- plot_options$myline
    names(use_col) <- plot_options$set
    names(use_dot) <- plot_options$set
    names(use_line) <- plot_options$set
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

