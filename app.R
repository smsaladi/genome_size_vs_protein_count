#
# Genome size vs gene count
# This is a Shiny web application that plots genome size vs gene count,
# staying up to date with NCBI's statistics
#
# Find out more about building applications with Shiny here:
#    http://shiny.rstudio.com/
#
#
# Author:  Shyam Saladi (saladi@caltech.edu)
# Date:    October 2018
# License: MIT
#

library(shiny)

library(tidyverse)
library(magrittr)
library(dplyrExtras)
library(svglite)
library(cowplot)

ui <- fluidPage(
  includeHTML("github_corner.html"),
  titlePanel("Genome size vs. protein count across NCBI genomes"),
  p("A log-log plot of the total number of annotated genes in genomes submitted to GenBank
    as a function of genome size (based on data provided by NCBI genome reports). The app stays
    up-to-date since it retrieves data from the NCBI FTP server upon loading."),
  # Select archaea/bacteria data source
  sidebarLayout(
    sidebarPanel(
      p("Please be patient. A few large files must first be retrieved from
         NCBI Genome Reports, so it may take a moment for the initial plot to appear."),
      br(),
      radioButtons("genomes_to_plot",
                   "Genomes to plot:",
                   c("All" = "all",
                     "Distinct names" = "distinct"),
                   selected = "distinct"
      ),
      radioButtons("bact_arch_filt",
                   "Bacterial/Archaeal filter:",
                   c("All" = "all",
                     "Reference" = "ref",
                     "Representative" = "rep")
      ),
      checkboxGroupInput("seq_status",
                         "Sequence type:",
                         c("Chromosome",
                           "Complete Genome",
                           "Contig",
                           "Scaffold"),
                         selected = c("Chromosome",
                                      "Complete Genome",
                                      "Contig",
                                      "Scaffold")
      ),
      br(),
      br(),
      downloadButton("downloadSVG", "Download SVG")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot", height = "700px")
    )
  )
)

base_url <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/"

# Download and parse table from NCBI's FTP site
read_genome_report <- function(fn) {
  df <- fn %>%
    paste0(base_url, ., ".txt") %>%
    read_tsv(
      col_types = cols(
        TaxID = col_integer(),
        `BioProject ID` = col_integer(),
        `Size (Mb)` = col_double(),
        `Size (Kb)` = col_double(),
        `GC%` = col_double(),
        Scaffolds = col_integer(),
        Genes = col_integer(),
        Proteins = col_integer(),
        `Release Date` = col_date(format = ""),
        `Modify Date` = col_date(format = ""),
        .default = col_character()
      ),
      na = c("-"),
      quoted_na = FALSE,
      quote = ""
    )
  
  if("Size (Kb)" %in% colnames(df))
    df %<>%
    mutate(`Size (Mb)` = `Size (Kb)` / 1000)
  
  df %>%
    # Remove if missing the numbers we want
    drop_na(`Size (Mb)`, Genes) %>%
    # Just keep necessary columns around
    transmute(
      `Organism name` = `#Organism/Name`,
      `BioProject ID`,
      `Group`,
      `SubGroup`,
      `Size (bp)` = `Size (Mb)` * 1000 * 1000,
      Genes,
      Proteins,
      Status,
      fn = fn
    ) %>%
    # remove identical Group/Size/Protein (very likely duplicates)
    distinct(SubGroup, `Size (bp)`, Proteins, .keep_all = TRUE)
}

#
# Variables to hold data, to avoid downloading dataset everytime the plot is updated
#

df_taxa <- "overview.txt" %>%
  paste0(base_url, .) %>%
  read_tsv(
    col_types = cols(.default = col_character()),
    na = c("-"),
    quoted_na = FALSE,
    quote = ""
  ) %>%
  rename(Superkingdom = Kingdom) %>%
  distinct(Group, SubGroup, Superkingdom) %>%
  # there are so few viroids (and for this plot, they can be considered as viruses)
  mutate_rows(Superkingdom == "Viroids", Superkingdom = "Viruses")

df_rep <- "prok_representative_genomes.txt" %>%
  paste0(base_url, .) %>%
  read_tsv(
    col_types = cols(.default = col_character()),
    na = c("-"),
    quoted_na = FALSE,
    quote = ""
  ) %>%
  transmute(`Organism name`,
            rep = TRUE) %>%
  distinct()

df_ref <- "prok_reference_genomes.txt" %>%
  paste0(base_url, .) %>%
  read_tsv(
    col_types = cols(.default = col_character()),
    na = c("-"),
    quoted_na = FALSE,
    quote = ""
  ) %>%
  transmute(`Organism name`,
            ref = TRUE) %>%
  distinct()

#
## Final/full dataframe used for plotting ##
#
df_genomes <- c("eukaryotes", "prokaryotes", "viruses") %>%
  map_dfr(read_genome_report) %>%
  left_join(df_taxa, by = c("Group", "SubGroup")) %>%
  left_join(df_rep, by = "Organism name") %>%
  left_join(df_ref, by = "Organism name") %>%
  # replace NA with FALSE
  mutate(rep = !is.na(rep),
         ref = !is.na(ref)) %>%
  # fix superkingdoms where Other/Other
  mutate_rows(Group == "Other" & fn == "eukaryotes", Superkingdom = "Eukaryota") %>%
  mutate_rows(Group == "Other" & fn == "viruses", Superkingdom = "Viruses") %>%
  # no good way of assigning Bacteria or Archaea, so remove
  filter(!(Group == "Other" & fn == "prokaryotes"))

# remove to reduce memory usage
rm(df_taxa, df_rep, df_ref)

# global variable that holds plot object
# so we don't need to recreate it to download svg
plot <- NULL

# Define server logic required to draw the plot
server <- function(input, output) {
  output$downloadSVG <- downloadHandler(
    filename = function() {
      paste0("genomes_vs_gene_count_", Sys.Date(), ".svg")
    },
    content = function(con) {
      ggsave(file = con, plot = plot,
             width = 10, height = 10, device = "svg")
    }
  )
  
  output$distPlot <- renderPlot({
    
    df <- df_genomes %>%
      filter(Status %in% input$seq_status)
    
    if(input$genomes_to_plot == "distinct")
      df %<>%
      group_by(`Superkingdom`, `Organism name`) %>%
      filter(`Size (bp)` == max(`Size (bp)`))
    
    if(input$bact_arch_filt == "ref")
      df %<>%
      filter(ref | Superkingdom %in% c("Eukaryota", "Viruses"))
    else if(input$bact_arch_filt == "rep")
      df %<>%
      filter(rep | Superkingdom %in% c("Eukaryota", "Viruses"))
    
    plot <<- df %>%
      filter(Superkingdom != "Bacteria", Superkingdom != "Viruses") %>%
      ggplot(aes(x = `Size (bp)`, y = Proteins, color = Superkingdom)) +
      # Do it this way so all bacteria and viruses are in back
      geom_point(shape = 1, alpha = 0.7,
                 data = df %>% filter(Superkingdom %in% c("Bacteria", "Viruses"))) +
      geom_point(shape=1) +
      scale_x_log10(
        breaks = 10^c(0:10),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
      ) +
      scale_y_log10(
        breaks = 10^c(0:10),
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
      ) +
      # from iwanthue.com
      scale_color_manual(values = c("Bacteria" = "#6ed1b4",
                                    "Archaea" = "#bf90d7",
                                    "Eukaryota" = "#e19775",
                                    "Viruses" = "#8fcd6b" #, "Viroids" = "#bf90d7"
      )) +
      annotation_logticks() +
      ggtitle("Genome size vs. protein count across NCBI genomes") +
      xlab('Genome size (bp)') +
      ylab("Protein coding genes (count)") +
      theme_bw(base_size = 18) +
      theme(
        aspect.ratio = 1,
        legend.title = element_blank(),
        legend.position = c(0.15, 0.80)
      )
    plot
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

