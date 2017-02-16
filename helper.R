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
        na_if(Inf)
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
      c(kDotOptions[1],
        kLineOptions[1],
        color_select,
        legend_nickname,
        legend_nickname,
        "none")
    
    # generate info for new file for loaded gene list(s)
    sapply(names(list_data$gene_file), function(g) {
      if (g != "common" & length(list_data$gene_file[[g]]$full) > 0 ) {
        enesg <- c(gene_names, list_data$gene_file[[g]]$full)
        enesg <- enesg[duplicated(enesg)]
        list_data$gene_file[[g]]$use <<- enesg
        list_data$gene_info[[g]][[legend_nickname]] <<-
          c(kDotOptions[1],
            kLineOptions[1],
            color_select,
            legend_nickname,
            legend_nickname,
            "none")
        
        # if 0 length message and remove file TODO
      }
    })
  }
  list_data
}
