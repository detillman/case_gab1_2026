## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(clusterProfiler)
library(Homo.sapiens)
library(readr)
library(DBI)
library(RSQLite)
library(RColorBrewer)
library(shiny)
library(plotly)
library(pheatmap)
library(readxl)
library(ggrepel)
library(UpSetR)
library(tidyverse)
go_reduce <- function (pathway_df, orgdb = "org.Hs.eg.db", threshold = 0.7,
                       scores = NULL, measure = "Wang")
{
  if (!measure %in% c("Resnik", "Lin", "Rel", "Jiang", "Wang")) {
    stop("Chosen measure is not one of the recognised measures, c(\"Resnik\", \"Lin\", \"Rel\", \"Jiang\", \"Wang\").")
  }
  if (measure == "Wang") {
    computeIC <- FALSE
  }
  else {
    computeIC <- TRUE
  }
  ont <- pathway_df %>% .[["go_type"]] %>% unique()
  if (any(!ont %in% c("BP", "CC", "MF"))) {
    stop("Column go_type does not contain the recognised sub-ontologies, c(\"BP\", \"CC\", \"MF\")")
  }
  go_similarity <- setNames(object = vector(mode = "list",
                                            length = length(ont)), nm = ont)
  for (i in 1:length(ont)) {
    print(stringr::str_c("Reducing sub-ontology: ", ont[i]))
    hsGO <- GOSemSim::godata(OrgDb = orgdb, ont = ont[i],
                             computeIC = computeIC)
    terms <- pathway_df %>% dplyr::filter(go_type ==
                                            ont[i]) %>% .[["go_id"]] %>% unique()
    sim <- GOSemSim::mgoSim(GO1 = terms, GO2 = terms, semData = hsGO,
                            measure = measure, combine = NULL)
    go_similarity[[i]] <- rrvgo::reduceSimMatrix(simMatrix = sim,
                                                 threshold = threshold, 
                                                 orgdb = orgdb, scores = scores) %>%
      tibble::as_tibble() %>% dplyr::rename(parent_id = parent,
                                            parent_term = parentTerm, 
                                            parent_sim_score = termDispensability)
  }
  go_sim_df <- go_similarity %>% qdapTools::list_df2df(col1 = "go_type")
  pathway_go_sim_df <- pathway_df %>% dplyr::inner_join(go_sim_df %>%
                                                          dplyr::select(go_type, 
                                                                        go_id = go, 
                                                                        contains("parent")),
                                                        by = c("go_type", "go_id")) %>% dplyr::arrange(go_type,
                                                                                                       parent_id, -parent_sim_score)
  return(pathway_go_sim_df)
}


## ----idIDRs-----------------------------------------------------------------------------------------------------------------
# defining idr parsing function
idr_parse <- function(idr_path) {

  idr_lines <- readLines(idr_path)

  # find header lines
  header_indices <- which(str_detect(idr_lines, "^#>"))
  header_lastIndex <- c(header_indices[-1] - 1, length(idr_lines))

  # parsing each fasta block
  map2_dfr(header_indices, header_lastIndex, function(i, j) {

    # defining each fasta block
    fasta_header <- str_remove(idr_lines[i], "^#>")
    fasta_block_temp <- idr_lines[(i + 1):j]
    fasta_block <- fasta_block_temp[nzchar(fasta_block_temp)]

    # storing information in temporary dataframe
    temp_df <- read_table2(I(fasta_block),
                           col_names = c("position", "residue",
                           "idr_probability"),
                           col_types = "icc")

    # extracting additional fasta information
    fasta_entry <- str_trim(fasta_header)
    uniprot_id <- str_match(fasta_entry, "^[^|]*\\|([^|]+)\\|")[,2]

    # making final dataframe
    temp_df %>%
      mutate(fasta_entry = fasta_entry,
             uniprot_id = uniprot_id) %>%
      select(uniprot_id, fasta_entry, position, residue, idr_probability)

  })

}

# running idr parse function
idr_parsed <- idr_parse("UP000005640_canonical_aiupred.txt")

# saving parsed idr file
write.csv(idr_parsed, "idr_parsed.csv", row.names = F)


## ----idPhosphoY-------------------------------------------------------------------------------------------------------------
# defining phosphoY parsing function
phosphoY_parse <- function(phosphoY_path) {

  phosphoY_lines <- if (length(phosphoY_path) == 1 &&
                        file.exists(phosphoY_path)) readLines(phosphoY_path,
                                                              warn = F) else phosphoY_path
  stopifnot(is.character(phosphoY_lines))

  # drop blank lines
  phosphoY_lines <- phosphoY_lines[nzchar(trimws(phosphoY_lines))]

  # find header lines
  header_indices <- which(str_detect(phosphoY_lines, "^>"))
  if (length(header_indices) == 0) stop("No '>' header lines found in file")
  header_lastIndex <- c(header_indices[-1] - 1, length(phosphoY_lines))

  # parsing each fasta block
  map2_dfr(header_indices, header_lastIndex, function(i, j) {

    # extracting additional fasta information
    fasta_entry <- str_trim(str_remove(phosphoY_lines[i], "^>"))
    uniprot_id <- str_match(fasta_entry, "^[^|]*\\|([^|]+)\\|")[,2]

    # defining each fasta block
    fasta_block_temp <- phosphoY_lines[(i + 1):j]

    # removing column header line if present
    fasta_block_clean <- fasta_block_temp[!str_detect(fasta_block_temp,
                                                      "^Position\\s+Residue\\s+PTMscores")]

    # final fasta block
    fasta_block <- fasta_block_clean[nzchar(trimws(fasta_block_clean))]

    # ignoring empty blocks
    if(length(fasta_block) == 0) return(tibble())

    # parse each line individually
    fasta_match <- str_match(fasta_block,
                             "^\\s*(\\d+)\\s+([A-Z])\\s+(Phosphotyrosine:[0-9.]+)\\s+(.*)\\s*$")

    # catching incorrectly parsed lines
    fasta_bad <- which(is.na(fasta_match[,1]))
    if (length(fasta_bad) > 0) {
      stop("Failed to parse some lines like: ",
           fasta_block[fasta_bad[1]])
    }

    # making final dataframe
    tibble(uniprot_id = uniprot_id,
           position = as.integer(fasta_match[,2]),
           residue = fasta_match[,3],
           ptm_score = fasta_match[,4],
           called = fasta_match[,5]) %>%
      transmute(uniprot_id,
                position,
                residue,
                phosphoY_probability = as.numeric(str_match(ptm_score,
                                                            "Phosphotyrosine:([0-9.]+)")[,2]))
  })

}

# running phosphoY parse function
phosphoY_parsed <- phosphoY_parse("UP000005640_canonical_musite_results.txt")

# saving parsed idr file
write.csv(phosphoY_parsed, "phosphoY_parsed.csv", row.names = F)


## ----idTopologies-----------------------------------------------------------------------------------------------------------
# defining topology parsing function
topology_parse <- function(topology_path) {

  # read in file with flexible number of columns
  topology_raw <- read_tsv(file = topology_path,
                           col_names = F,
                           comment = "#",
                           na = c("", "NA"),
                           show_col_types = F,
                           progress = F,
                           trim_ws = T)

  # drop separator lines
  topology_clean <- topology_raw %>%
    filter(.data$X1 != "//")

  # only keep first four columns
  if (ncol(topology_clean) < 4) stop("Topology file has fewer than four columns")

  # making final dataframe
  topology_final <- topology_clean %>%
    transmute(fasta_entry = as.character(.data$X1),
              topology_type = as.character(.data$X2),
              topology_start = as.character(.data$X3),
              topology_end = as.character(.data$X4),
              uniprot_id = str_match(fasta_entry, "^[^|]*\\|([^|]+)\\|")[, 2]) %>%
    select(-fasta_entry) %>%
    relocate(uniprot_id, .before = topology_type)
}

# running topology parse function
topology_parsed <- topology_parse("TMRs.gff3") %>%
  mutate(topology_type = case_when(topology_type == "outside" ~ "Outside",
                                   topology_type == "TMhelix" ~ "TM Helix",
                                   topology_type == "inside" ~ "Inside",
                                   topology_type == "signal" ~ "Signal",
                                   topology_type == "periplasm" ~ "Periplasm",
                                   topology_type == "Beta sheet" ~ "Beta Sheet"))

# saving parsed topology file
write.csv(topology_parsed, "topology_parsed.csv", row.names = F)


## ----tableMaking------------------------------------------------------------------------------------------------------------
# reading in tables
idr_parsed <- read.csv("idr_parsed.csv")
phosphoY_parsed <- read.csv("phosphoY_parsed.csv")
topology_parsed <- read.csv("topology_parsed.csv")

# finding itms
find_itms <- function(protein_sequence) {

  # empty tibble for no hits
  empty_tibble <- tibble(motif_type = character(),
                         motif_start = integer(),
                         motif_end = integer(),
                         motif_sequence = character())

  # define protein parameters
  upper_sequence <- toupper(protein_sequence)
  protein_length <- nchar(protein_sequence)

  # scan function
  motif_scanner <- function(motif_type, motif_width, motif_pattern) {

    map_dfr(motif_width, function(w) {

      # ignore proteins with small sequences
      if(protein_length < w) return(empty_tibble[0, ])

      # define motif parameters
      motif_starts <- 1:(protein_length - w + 1L)
      motif_kmers <- substring(upper_sequence, motif_starts, motif_starts + w - 1L)
      motif_hits <- str_detect(motif_kmers, motif_pattern)

      # ignore proteins with no hits
      if (!any(motif_hits)) return(empty_tibble[0, ])

      # save scan results
      tibble(motif_type = motif_type,
             motif_start = motif_starts[motif_hits],
             motif_end = motif_starts[motif_hits] + w - 1L,
             motif_sequence = motif_kmers[motif_hits])
    })

  }

  # scan sequences
  itim_df <- motif_scanner("ITIM", 6L,  "^[SIVL].Y..[IVL]$")
  itsm_df <- motif_scanner("ITSM", 6L,  "^T.Y..[VI]$")
  itam_df <- motif_scanner("ITAM", 14:20, "^Y..[IL].{6,12}Y..[IL]$")

  # save all results
  bind_rows(itim_df, itsm_df, itam_df) %>%
    arrange(motif_start, motif_end, motif_type)

}

# merging per-residue tables
residue_scores <- idr_parsed %>%
  left_join(phosphoY_parsed,
            by = c("uniprot_id", "position", "residue"))

# making protein-level table
protein_table <- residue_scores %>%
  group_by(uniprot_id, fasta_entry) %>%
  summarize(protein_length = max(position),
            protein_sequence = paste0(residue[order(position)], collapse = ""),
            .groups = "drop") %>%
  mutate(gene_id = str_match(fasta_entry, "\\bGN=([^ ]+)")[,2])

# finding itms
itm_hits <- protein_table %>%
  select(uniprot_id, protein_sequence) %>%
  mutate(motif_hits = map(protein_sequence, find_itms)) %>%
  select(-protein_sequence) %>%
  unnest(motif_hits)


## ----databaseBuilding-------------------------------------------------------------------------------------------------------
# define sql connection location
sql_connect <- dbConnect(SQLite(), "protein_features.sqlite")

# writing tables to sql
dbWriteTable(sql_connect, "proteins", protein_table, overwrite = TRUE)
dbWriteTable(sql_connect, "residue_scores", residue_scores, overwrite = TRUE)
dbWriteTable(sql_connect, "topologies", topology_parsed, overwrite = TRUE)
dbWriteTable(sql_connect, "itm_hits", itm_hits, overwrite = TRUE)

# creating sql indices
dbExecute(sql_connect, "CREATE INDEX IF NOT EXISTS idx_residue_uniprot ON residue_scores(uniprot_id)")
dbExecute(sql_connect, "CREATE INDEX IF NOT EXISTS idx_topology_uniprot ON topologies(uniprot_id)")
dbExecute(sql_connect, "CREATE INDEX IF NOT EXISTS idx_itm_uniprot ON itm_hits(uniprot_id)")

# disconnecting from database
dbDisconnect(sql_connect)


## ----shinyRunning-----------------------------------------------------------------------------------------------------------
# define sql connection location
sql_connect <- dbConnect(SQLite(), "protein_features.sqlite")

# making function to get per-protein data
protein_data <- function(sql_connection, uniprot_id) {
  
  # get protein id
  protein_query <- dbGetQuery(sql_connection, 
                              "SELECT * FROM proteins WHERE uniprot_id = ?", 
                              params = list(uniprot_id))
  
  # ignore non-existent protein ids
  if (nrow(protein_query) == 0) return(NULL)
  
  # get protein information
  protein_residues <- dbGetQuery(sql_connection, 
                                 "SELECT position, residue, 
                                 idr_probability, phosphoY_probability FROM 
                                 residue_scores WHERE 
                                 uniprot_id = ? ORDER BY position",
                         params = list(uniprot_id))
  protein_topology <- dbGetQuery(sql_connection, 
                                 "SELECT topology_type, 
                                 topology_start, 
                                 topology_end FROM 
                                 topologies WHERE 
                                 uniprot_id = ?",
                     params = list(uniprot_id))
  protein_itms <- dbGetQuery(sql_connection, 
                               "SELECT motif_type, 
                               motif_start, 
                               motif_end, 
                               motif_sequence FROM 
                               itm_hits WHERE 
                               uniprot_id = ?",
                       params = list(uniprot_id))
  
  # summarize protein information
  list(protein = protein_query, residues = protein_residues, 
       topology = protein_topology, itms = protein_itms)
}

# making custom fill palette
fill_df <- data.frame(fill_type = c("Outside", "TM Helix", "Inside", "Signal",
                                    "Periplasm", "Beta Sheet", 
                                    "ITAM", "ITIM", "ITSM")) %>%
  mutate(fill_type = as.factor(fill_type))
fill_colors <- c("#CCEBC5", "#80B1D3", "#8DD3C7", "#BC80BD", "#BEBADA",
                 "#FDB462", "#FCCDE5", "#FB8072", "#B3DE69")
names(fill_colors) <- levels(fill_df$fill_type)
fill_palette <- scale_fill_manual(name = "", values = fill_colors)

# making function to fix track names
fix_tracks <- function(old_plot) {
  for (i in seq_along(old_plot$x$data)) {
    
    # get names and ignore tracks with no names
    nm <- old_plot$x$data[[i]]$name
    if (is.null(nm)) next

    # clean names
    nm2 <- str_replace(nm, "^\\(([^,]+),\\s*\\d+\\)$", "\\1")
    old_plot$x$data[[i]]$name <- nm2

    # clean hover
    if (!is.null(old_plot$x$data[[i]]$legendgroup)) {
      old_plot$x$data[[i]]$legendgroup <- nm2
    }
  }
  old_plot
}

# making function to plot protein information
protein_plot <- function(residues, topology, itms, title_text) {
  
  # annotating phosphotyrosines
  plot_residues <- residues %>%
    mutate(is_tyrosine = residue == "Y",
           idr_state = if_else(idr_probability >= 0.5, TRUE, FALSE),
           idr_segment = cumsum(idr_state != lag(idr_state, 
                                                 default = first(idr_state))),
           is_phosphoY = case_when(phosphoY_probability < 0.5 ~ F,
                                   phosphoY_probability >= 0.5 ~ T),
           hover = paste0("Pos: ", position, 
                          "<br>AA: ", residue,
                          "<br>IDR: ", round(idr_probability, 2),
                          ifelse(is.na(phosphoY_probability), "",
                                 paste0("<br>pY: ", 
                                        round(phosphoY_probability, 2))))) %>%
    mutate(idr_state = fct_relevel(as.factor(idr_state), "TRUE", "FALSE"),
           is_phosphoY = fct_relevel(as.factor(is_phosphoY), "TRUE", "FALSE"))
  
  # making plot
  to_plot <- ggplot() +
    
    # idr track
    geom_line(data = plot_residues, aes(x = position, 
                                   y = idr_probability + 2.5,
                                   group = idr_segment,
                                   text = hover,
                                   color = idr_state),
              linewidth = 0.5) +
    geom_hline(yintercept = 3, linewidth = 0.125, linetype = "dashed") +
    
    # phosphoY track
    geom_point(data = plot_residues %>%
                 filter(is_tyrosine, !is.na(phosphoY_probability)),
               aes(x = position, y = phosphoY_probability + 1.5, text = hover,
                   color = is_phosphoY),
               size = 1.5) +
    scale_color_manual(values = c("red", "grey50")) +
    geom_hline(yintercept = 2, linewidth = 0.125, linetype = "dashed")
    
  # itm track
  if(!is.null(itms) && nrow(itms) > 0) {
    to_plot <- to_plot + 
      geom_rect(data = itms,
                aes(xmin = motif_start, xmax = motif_end, 
                    ymin = 0.85, ymax = 1.15, fill = motif_type),
                alpha = 0.7)
  }
  
  # topology_track
  if(!is.null(topology) && nrow(topology) > 0) {
    to_plot <- to_plot +
      geom_rect(data = topology,
                aes(xmin = topology_start, xmax = topology_end,
                    ymin = -0.15, ymax = 0.15, fill = topology_type),
                alpha = 0.7)
  }
  
  # plotting parameters
  to_plot <- to_plot +
    fill_palette +
    scale_y_continuous(breaks = c(0, 1, 2, 3),
                       labels = c("Topology", "ITMs", "pY", "IDR")) +
    labs(x = "Position", y = "", fill = "", color = "", title = title_text) +
    geom_hline(yintercept = 0.5, color = "black", linewidth = 0.125) +
    geom_hline(yintercept = 1.5, color = "black", linewidth = 0.125) +
    geom_hline(yintercept = 2.5, color = "black", linewidth = 0.125) +
    geom_hline(yintercept = 3.5, color = "black", linewidth = 0.125) +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
  
  # plot
  final_plot <- ggplotly(to_plot, tooltip = "text") %>%
    layout(hovermode = "closest"); fix_tracks(final_plot)
  
}

# making shiny ui
shiny_ui <- fluidPage(
  textInput("uniprot", "UniProt ID", value = "Q15116"),
  actionButton("go", "Load"),
  plotlyOutput("tracks", height = "600px"),
  verbatimTextOutput("msg")
)

# making shiny server
shiny_server <- function(input, output, session) {
  
  # connect to sql database
  sql_connect <- dbConnect(RSQLite::SQLite(), "protein_features.sqlite")
  onStop(function() dbDisconnect(sql_connect))

  # get protein data
  shiny_data <- eventReactive(input$go, {
    req(input$uniprot)
    protein_data(sql_connect, toupper(trimws(input$uniprot)))
  })

  # plot output
  output$msg <- renderPrint({
    if (is.null(shiny_data())) "No protein found for that UniProt ID." else shiny_data()$protein
  })
  output$tracks <- renderPlotly({
    req(shiny_data())
    used_protein <- shiny_data()$protein[1, ]
    protein_name <- used_protein$uniprot_id
    gene_name <- used_protein$gene_id
    title_text <- if(!is.na(gene_name) && nzchar(gene_name)) {
      paste0("UniProt ", protein_name, " (", gene_name, ")")
    } else {
      protein_name
    }
    
    
    protein_plot(shiny_data()$residues, 
                 shiny_data()$topology, 
                 shiny_data()$itms,
                 title_text)
  })
}

# run shiny
shinyApp(shiny_ui, shiny_server)


## ----egfrGrb2Analysis-------------------------------------------------------------------------------------------------------
# load interaction dataframes
egfr_hpa <- read_excel("egfr_grb2_interactors.xlsx", 
                       sheet = "egfr_hpa_uniprot") %>%
  mutate(comp = "egfr_hpa")
grb2_hpa <- read_excel("egfr_grb2_interactors.xlsx", 
                       sheet = "grb2_hpa_uniprot") %>%
  mutate(comp = "grb2_hpa")

# combine interaction dataframes
egfr_grb2_interactions <- bind_rows(egfr_hpa %>%
                                      select(uniprot_id, comp), 
                                    grb2_hpa %>%
                                      select(uniprot_id, comp))

# define sql connection location
sql_connect <- dbConnect(SQLite(), "protein_features.sqlite")

# function to extract residue information
get_res_info <- function(con, uniprot_ids) {

  # write temporary table with uniprot ids
  dbWriteTable(
    con, "input_ids",
    data.frame(uniprot_id = uniprot_ids),
    temporary = TRUE, overwrite = TRUE
  )

  # define sql query
  q <- "
SELECT
  r.uniprot_id,
  r.position,
  t.topology_type AS topology,
  r.idr_probability,
  r.phosphoY_probability,
  r.residue
FROM residue_scores r
JOIN input_ids ids
  ON ids.uniprot_id = r.uniprot_id
LEFT JOIN topologies t
  ON t.uniprot_id = r.uniprot_id
 AND r.position BETWEEN t.topology_start AND t.topology_end
ORDER BY r.uniprot_id, r.position;
  "

  # get sql query
  dbGetQuery(con, q)
}

# run sql function
egfr_grb2_residues <- get_res_info(sql_connect, 
                                   unique(egfr_grb2_interactions$uniprot_id))

# summarize sql data
egfr_grb2_summary <- egfr_grb2_residues %>%
mutate(is_cyto = case_when(topology == "Inside" ~ T, T ~ F),
       is_cytoIDR = case_when(is_cyto == T & idr_probability > 0.5 ~ T,
                              T ~ F),
       is_cytoY = case_when(is_cyto == T & residue == "Y" ~ T,
                            T ~ F),
       is_cytoPY = case_when(is_cyto == T & phosphoY_probability > 0.5 ~ T,
                             T ~ F),
       is_cytoIDRy = case_when(is_cytoIDR == T & is_cytoY == T ~ T, T ~ F),
       is_cytoIDRpY = case_when(is_cytoIDR == T & is_cytoPY == T ~ T, T ~ F)) %>%
  group_by(uniprot_id) %>%
  summarize(cyto = sum(is_cyto) + 0.5,
            cytoIDR = sum(is_cytoIDR) + 0.5,
            cytoY = sum(is_cytoY),
            cytoPY = sum(is_cytoPY),
            cytoIDRy = sum(is_cytoIDRy),
            cytoIDRpY = sum(is_cytoIDRpY)) %>%
  ungroup()

# making name dataframe
egfr_grb2_summary_hpa_or <- egfr_grb2_summary %>%
  mutate(plot_color = case_when(uniprot_id %in% egfr_hpa$uniprot_id &
                                  uniprot_id %in% grb2_hpa$uniprot_id ~ "EGFR/GRB2",
                                uniprot_id %in% egfr_hpa$uniprot_id & 
                                  !(uniprot_id %in% grb2_hpa$uniprot_id) ~ "EGFR-Only",
                                !(uniprot_id %in% egfr_hpa$uniprot_id) &
                                  uniprot_id %in% grb2_hpa$uniprot_id ~ "GRB2-Only",
                                T ~ "Neither")) %>%
  filter(plot_color != "Neither") %>%
  left_join(bind_rows(egfr_hpa) %>% distinct(uniprot_id, gene_name),
            by = "uniprot_id") %>%
  mutate(gene_name = case_when(is.na(gene_name) == T ~ uniprot_id,
                               T ~ gene_name),
         plot_label = case_when(cytoIDRpY > 6 ~ gene_name, T ~ NA)) %>%
  mutate(plot_label = case_when(uniprot_id == "P21860" ~ "EGFR3",
                                uniprot_id == "Q15303" ~ "EGFR4",
                                uniprot_id == "Q9UQC2" ~ "GAB2",
                                uniprot_id == "P04626" ~ "EGFR2",
                                uniprot_id == "Q05397" ~ "FAK",
                                uniprot_id == "P62993" ~ "GRB2",
                                uniprot_id == "O00401" ~ "N-WASP",
                                uniprot_id == "Q07889" ~ "SOS1",
                                uniprot_id == "Q07890" ~ "SOS2",
                                uniprot_id == "P42336" ~ "PIK3CA",
                                uniprot_id == "P27986" ~ "PIK3R1",
                                uniprot_id == "P29353" ~ "SHC1",
                                uniprot_id == "Q06124" ~ "PTPN11",
                                uniprot_id == "Q13905" ~ "RAPGEF1",
                                uniprot_id == "P22681" ~ "CBL",
                                uniprot_id == "Q9UJM3" ~ "ERRFI1",
                                uniprot_id == "Q13094" ~ "LCP2",
                                uniprot_id == "O00750" ~ "PIK3C2B",
                                uniprot_id == "Q07666" ~ NA, #KHDRBS1
                                uniprot_id == "Q14247" ~ NA, #CTTN
                                T ~ plot_label))

# subsetting for only egfr/grb2 interactors
egfr_grb2_summary_hpa_or_subset <- egfr_grb2_summary_hpa_or %>%
  filter(plot_color == "EGFR/GRB2" | plot_label == "PIK3CA" | plot_label == "GRB2" |
           plot_label == "PTPN11" | plot_label == "SHC1") %>%
  mutate(plot_color = case_when(plot_color == "EGFR-Only" ~ "gray",
                                plot_color == "GRB2-Only" ~ "gray",
                                T ~ plot_color))

## supplemental figure 1f
# plotting subset of only egfr/grb2 interactors
gab1_selection <- ggplot(egfr_grb2_summary_hpa_or_subset, aes(x = cytoIDR - 0.5, y = cytoIDRpY,
                                     color = plot_color, label = plot_label)) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  geom_text_repel(min.segment.length = 0.1, max.overlaps = 10,
                   size = 2, segment.size = 0.25) +
  #scale_x_continuous(trans = "log10") +
  labs(x = "Cytoplasmic IDR Length", y = "# pY in Cytoplasmic IDR",
       color = "") +
  scale_color_manual(values = c("black", "gray50")) +
  theme_classic(base_size = 7, base_family = "sans") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); gab1_selection
ggsave("figure_gab1Selection.pdf", gab1_selection, units = "in",
       width = 3.33, height = 2.5)


## ----ncprCalc---------------------------------------------------------------------------------------------------------------
# make dataframe for ncpr calculations
ncpr_df <- data.frame(gene_name = c("EGFR", "GAB1", "LAT", "PD1"),
                      uniprot_id = c("P00533", "Q13480", "O43561", "Q15116"),
                      cyto_idr_seq = c("LQRYSSDPTGALTEDSIDDTFLPVPEYINQSVPKRPAGSVQNPVYHNQPLNPAPSRDPHYQDPHSTAVGNPEYLNTVQPTCVNSTFDSPAHWAQKGSHQISLDNPDYQQDFFPKEAKPNGIFKGSTAENAEYLRVAPQSSEFIGA", "GFNPTEEDPVKPPGSSLQAPADLPLAINTAPPSTQADSSSATLPPPYQLINVPPHLETLGIQEDPQDYLLLINCQSKKPEPTRTHADSAKSTSSETDCNDNVPSHKNPASSQSKHGMNGFFQQQMIYDSPPSRAPSASVDSSLYNLPRSYSHDVLPKVSPSSTEADGELYVFNTPSGTSSVETQMRHVSISYDIPPTPGNTYQIPRTFPEGTLGQTSKLDTIPDIPPPRPPKPHPAHDRSPVETCSIPRTASDTDSSYCIPTAGMSPSRSNTISTVDLNKLRKDASSQDCYDIPRAFPSDRSSSLEGFHNHFKVKNVLTVGSVSSEELDENYVPMNPNSPPRQHSSSFTEPIQEANYVPMTPGTFDFSSFGMQVPPPAHMGFRSSPKTPPRRPVPVADCEPPPVDRNLKPDRKVKPAPLEIKPLPEWEELQAPVRSPITRSFARDSSRFPMSPRPDSVHSTTSSSDSHDSEENYVPMNPNLSSEDPNLFGSNSLDGGSSPMIKPKGDKQVEYLDLDLDSGKSTPPRKQKSSGSGSSVADERVDYVVVDQQKTLALKSTREAWTDGRQSTESETPAKSVK", "HCHRLPGSYDSTSSDSLYPRGIQFKRPHTVAPWPPAYPPVTSYPPLSQPDLLPIPRSPQPLGGSHRTPSSRRDSDGANSVASYENEGASGIRGAQAGWGVWGPSWTRLTPVSLPPEPACEDADEDEDDYHNPGYLVVLPDSTPATSTAAPSAPALSTPGIRDSAFSMESIDDYVNVPESGESAEASLDGSREYVNVSQELHPGAAKTEPAALSSQEAEEVEEEGAPDYENLQELN", "CSRAARGTIGARRTGQPLKEDPSAVPVFSVDYGELDFQWREKTPEPPVPCVPEQTEYATIVFPSGMGTSSPARRGSADGPRSAQPLRPEDGHCSWPL"))

# function to calculate ncpr for a single window
calc_ncpr <- function(window_seq) {
  chars <- strsplit(window_seq, "")[[1]]
  pos <- sum(chars %in% c("K", "R"))
  neg <- sum(chars %in% c("D", "E"))
  (pos - neg) / length(chars)
}

# function to calculate per-residue mean ncpr for one sequence
sliding_ncpr_per_residue <- function(seq, window_size = 5) {
  n <- nchar(seq)
  chars <- strsplit(seq, "")[[1]]

  # if sequence shorter than window, just return overall ncpr for every residue
  if (n < window_size) {
    val <- calc_ncpr(seq)
    return(tibble(
      position = seq_len(n),
      residue  = chars,
      mean_NCPR = val
    ))
  }

  # generate start indices for each window
  starts <- 1:(n - window_size + 1)

  # data frame: one row per window, listing which residue positions it covers
  window_df <- map_dfr(starts, function(s) {
    e <- s + window_size - 1
    window_seq <- substr(seq, s, e)
    tibble(
      window_id = s,
      ncpr = calc_ncpr(window_seq),
      position = s:e
    )
  })

  # average ncpr across all windows that cover each residue position
  window_df %>%
    group_by(position) %>%
    summarise(mean_NCPR = mean(ncpr), .groups = "drop") %>%
    mutate(residue = chars[position]) %>%
    select(position, residue, mean_NCPR)
}

# calculate ncpr
ncpr_result <- ncpr_df %>%
  mutate(ncpr_data = map(cyto_idr_seq, sliding_ncpr_per_residue, window_size = 5)) %>%
  select(gene_name, uniprot_id, ncpr_data) %>%
  unnest(ncpr_data) %>%
  mutate(fl_position = case_when(gene_name == "EGFR" ~ position + 1065,
                                 gene_name == "GAB1" ~ position + 115,
                                 gene_name == "LAT" ~ position + 27,
                                 gene_name == "PD1" ~ position + 191))

# making combined pd1 and gab1 matrix
combo_ncpr_matrix <- ncpr_result %>%
  filter(gene_name == "PD1" | gene_name == "GAB1") %>%
  mutate(keep_residue = case_when(gene_name == "PD1" & fl_position > 207 ~ T,
                                  gene_name == "GAB1" & fl_position > 613 ~ T,
                                  T ~ F)) %>%
  filter(keep_residue == T) %>%
  mutate(plot_position = case_when(gene_name == "PD1" ~ fl_position - 288,
                                   gene_name == "GAB1" ~ fl_position - 694)) %>%
  select(gene_name, mean_NCPR, plot_position) %>%
  pivot_wider(names_from = plot_position, values_from = mean_NCPR) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
combo_ncpr_labels <- ifelse(as.numeric(colnames(combo_ncpr_matrix)) %% 10 == 0,
                        colnames(combo_ncpr_matrix), "")

## supplemental figure 2c
# plotting combo ncpr heatmap
pheatmap(combo_ncpr_matrix, cluster_rows = F, cluster_cols = F,
         breaks = seq(-0.6, 0.6, length.out = 101),
         border_color = NA,
         color = colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(100),
         labels_col = combo_ncpr_labels,
         fontsize = 7,
         show_colnames = F,
         filename = "pd1_ncpr_heatmap.pdf",
         width = 3, height = 0.75)


## ----idCandidates-----------------------------------------------------------------------------------------------------------
# define sql connection location
sql_connect <- dbConnect(SQLite(), "protein_features.sqlite")

# make upsest plot sql query
pd1_upset_query <- "
SELECT i.uniprot_id, i.motif_type, i.motif_start, i.motif_end,
i.motif_start + 2 AS y_position,
r.residue, r.idr_probability, r.phosphoY_probability,
t.topology_type, t.topology_start, t.topology_end
FROM itm_hits i
JOIN residue_scores r
ON r.uniprot_id = i.uniprot_id
AND r.position = (i.motif_start + 2)
JOIN topologies t
ON t.uniprot_id = i.uniprot_id
AND (i.motif_start + 2) BETWEEN t.topology_start AND topology_end
WHERE i.motif_type IN ('ITIM', 'ITSM')
ORDER BY i.uniprot_id, y_position
"

# submit upset plot sql query
pd1_upset_hits <- dbGetQuery(sql_connect, pd1_upset_query)

# making pd1 upset plot dataframe
pd1_upset_df <- pd1_upset_hits %>%
  mutate(Cytoplasmic = case_when(topology_type == "Inside" ~ 1, T ~ 0),
         Disordered = case_when(idr_probability > 0.5 ~ 1, T ~ 0),
         Phosphorylated = case_when(phosphoY_probability > 0.5 ~ 1, T ~ 0),
         is_empty = case_when((Cytoplasmic == F & Disordered == F & Phosphorylated == F) ~ 1,
                              T ~ 0)) %>%
  unite(protein_itm, uniprot_id, motif_start, remove = F) %>%
  select(protein_itm, Cytoplasmic, Disordered, Phosphorylated, is_empty)

## figure 2g
# plotting pd1 upsets
pdf("pd1_upset_empty.pdf", height = 3.41, width = 3.41)
print(upset(pd1_upset_df,
      nintersects = NA, main.bar.color = "black", sets.bar.color = "black",
      order.by = "degree", matrix.color = "black", decreasing = T))
dev.off()

pdf("pd1_upset_noEmpty.pdf", height = 2.25, width = 3)
print(upset(pd1_upset_df %>% select(-is_empty),
      nintersects = NA, main.bar.color = "black", sets.bar.color = "black",
      order.by = "degree", matrix.color = "black", decreasing = T))
dev.off()

# make intersection sql query
pd1_query <- "
SELECT i.uniprot_id, i.motif_type, i.motif_start, i.motif_end,
i.motif_start + 2 AS y_position,
r.residue, r.idr_probability, r.phosphoY_probability,
t.topology_type, t.topology_start, t.topology_end
FROM itm_hits i
JOIN residue_scores r
ON r.uniprot_id = i.uniprot_id
AND r.position = (i.motif_start + 2)
JOIN topologies t
ON t.uniprot_id = i.uniprot_id
AND (i.motif_start + 2) BETWEEN t.topology_start AND topology_end
WHERE i.motif_type IN ('ITIM', 'ITSM')
AND r.idr_probability > 0.5
AND r.phosphoY_probability > 0.5
AND t.topology_type = 'Inside'
ORDER BY i.uniprot_id, y_position
"

# submitting intersection query
pd1_hits <- dbGetQuery(sql_connect, pd1_query)

# calculating number of valid itms for each candidate protein
pd1_itmTotals <- pd1_hits %>%
  group_by(uniprot_id) %>%
  summarize(valid_itms = n()) %>%
  arrange(-valid_itms) %>%
  mutate(itmTotal_score = case_when(valid_itms > 1 ~ 1,
                                    valid_itms == 1 ~ 0))


## ----goCandidates-----------------------------------------------------------------------------------------------------------
# reading in human universe of proteins (20,416 proteins)
human_universe <- read.csv("topology_parsed.csv") %>%
  distinct(uniprot_id)

# reading in itm universe of proteins (10,028 proteins)
itm_universe <- pd1_upset_hits %>%
  distinct(uniprot_id)

# reading in all valid itm hits (588 proteins)
pd1_allHits <- pd1_itmTotals %>%
  distinct(uniprot_id)

# reading in multiple valid itm hits (93 proteins)
pd1_multiHits <- pd1_itmTotals %>%
  filter(valid_itms > 1)

# running go analysis
go_allHits_hUniverse <- enrichGO(pd1_allHits$uniprot_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "UNIPROT",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       readable = T,
                       pool = T,
                       universe = human_universe$uniprot_id) %>%
  as.data.frame() %>%
  mutate(comp = "allHits_hUniverse")
go_multiHits_hUniverse <- enrichGO(pd1_multiHits$uniprot_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "UNIPROT",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       readable = T,
                       pool = T,
                       universe = human_universe$uniprot_id) %>%
  as.data.frame() %>%
  mutate(comp = "multiHits_hUniverse")
go_allHits_iUniverse <- enrichGO(pd1_allHits$uniprot_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "UNIPROT",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       readable = T,
                       pool = T,
                       universe = itm_universe$uniprot_id) %>%
  as.data.frame() %>%
  mutate(comp = "allHits_iUniverse")
go_multiHits_iUniverse <- enrichGO(pd1_multiHits$uniprot_id,
                       OrgDb = org.Hs.eg.db,
                       keyType = "UNIPROT",
                       ont = "ALL",
                       pAdjustMethod = "BH",
                       readable = T,
                       pool = T,
                       universe = itm_universe$uniprot_id) %>%
  as.data.frame() %>%
  mutate(comp = "multiHits_iUniverse")


# combing and saving go dataframe for later retrieval
pd1_go_df <- bind_rows(go_allHits_hUniverse, go_multiHits_hUniverse,
                       go_allHits_iUniverse, go_multiHits_iUniverse) %>%
  separate(GeneRatio, c("Num", "Denom")) %>%
  mutate(GeneRatio = as.numeric(Num) / as.numeric(Denom))
write.csv(pd1_go_df, file = "pd1_go_df.csv", row.names = F)

# loading and reducing go dataframe
pd1_go_df <- read.csv("pd1_go_df.csv")
pd1_go_df_comps <- pd1_go_df %>%
  distinct(comp)
pd1_go_df_reduced <- data.frame()
for (i in 1:nrow(pd1_go_df_comps)) {

  # get appropriate comparison
  plot_comp <- pd1_go_df_comps[i, ]
  plot_go_df <- pd1_go_df %>%
    filter(comp == plot_comp) %>%
    arrange(qvalue) %>%
    mutate(trans_qvalue = -log10(qvalue))

  # ensure more than one row per ontology
  go_df_reduceTest <- plot_go_df %>%
    group_by(ONTOLOGY) %>%
    summarize(ont_count = n()) %>%
    filter(ont_count == 1)
  rows_to_bind <- plot_go_df %>%
    filter(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY)
  plot_go_df <- plot_go_df %>%
    filter(!(ONTOLOGY %in% go_df_reduceTest$ONTOLOGY))

  # format for reduce function
  go_df_preReduce <- plot_go_df %>%
    select(ONTOLOGY, ID) %>%
    rename(go_type = ONTOLOGY,
           go_id = ID)
  go_df_scores <- plot_go_df$trans_qvalue %>%
    setNames(plot_go_df$ID)

  # run reduce function
  go_df_reduce <- go_reduce(go_df_preReduce,
                            orgdb = "org.Hs.eg.db",
                            threshold = 0.7,
                            scores = go_df_scores,
                            measure = "Wang")
  plot_go_df_reduce <- plot_go_df %>%
    filter(ID %in% go_df_reduce$parent_id) %>%
    bind_rows(rows_to_bind)

  # saving reduced go dataframe
  pd1_go_df_reduced <- bind_rows(pd1_go_df_reduced, plot_go_df_reduce)
}
write.csv(pd1_go_df_reduced, file = "pd1_go_df_reduced.csv", row.names = F)

# plotting reduced go dataframes
pd1_go_df_reduced <- read.csv("pd1_go_df_reduced.csv")
pd1_go_df_comps <- pd1_go_df_reduced %>%
  distinct(comp)
for (i in 1:nrow(pd1_go_df_comps)) {
  
  # get appropriate comparison
  plot_comp <- pd1_go_df_comps[i, ]
  plot_go_df <- pd1_go_df_reduced %>%
    filter(comp == plot_comp) %>%
    arrange(qvalue) %>%
    group_by(ONTOLOGY) %>%
    slice(1:10) %>%
    ungroup() %>%
    mutate(Description = paste0(ONTOLOGY,
                                ": ",
                                Description),
           Description = fct_reorder2(Description,
                                      -qvalue,
                                      qvalue))
  
  print(ggplot(plot_go_df, aes(x = GeneRatio, y = Description, 
                               fill = qvalue, size = Count)) +
          geom_point(color = "black", pch = 21) +
          facet_wrap(~ ONTOLOGY) + 
          labs(y = "", title = plot_comp) +
          theme_classic(base_size = 10, base_family = "sans") +
          scale_fill_distiller(palette = "RdBu", direction = 1) +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black")))
}

## figure 2h
pd1_go_allValidItims_iUniverse <- ggplot(pd1_go_df_reduced %>%
         filter(comp == "allHits_iUniverse") %>%
         arrange(qvalue) %>%
           mutate(Description = case_when(Description == "signaling adaptor activity" ~
                                            "signaling adaptor",
                                          Description == "adaptive immune response" ~
                                            "adaptive immunity",
                                          Description == "actin filament-based process" ~
                                            "actin-based process",
                                          Description == "negative regulation of antigen receptor-mediated signaling pathway" ~ "antigen signaling",
                                          Description == "actin cytoskeleton organization" ~
                                            "actin organization",
                                          T ~ Description)) %>%
           slice(1:10) %>%
         mutate(Description = paste0(ONTOLOGY,
                                ": ",
                                Description),
           Description = fct_reorder2(Description,
                                      -qvalue,
                                      qvalue)), 
       aes(x = GeneRatio, y = Description, 
           fill = qvalue, size = Count)) +
          geom_point(color = "black", pch = 21) +
          labs(y = "", title = "GO Analysis: 1+ pITIM/pITSM in Cytoplasmic IDR") +
          theme_classic(base_size = 7, base_family = "sans") +
          scale_fill_distiller(palette = "RdBu", direction = 1) +
  scale_y_discrete(labels = scales::label_wrap(16)) +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black"),
                legend.position = "inside",
                legend.position.inside = c(0.9, 0.9)); pd1_go_allValidItims_iUniverse
ggsave("pd1_go_allValidItims_iUniverse.pdf", pd1_go_allValidItims_iUniverse,
       units = "in", height = 1.91*1.25, width = 2.16)


## ----sessionInfo------------------------------------------------------------------------------------------------------------
sessionInfo()

