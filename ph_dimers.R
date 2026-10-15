## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
library(ggpubr)
library(readxl)
library(biomaRt)
library(ggrepel)
library(tidyverse)


## ----uniprotLoad------------------------------------------------------------------------------------------------------------
# reading in uniprot data for original alphaPulldown run with ph-containing human proteins
original_ph <- read_excel("alphaServer/phDimers/phProteins.xlsx") %>%

  # focus on first gene name
  mutate(gene_name = gsub(" .*", "", gene_names), .keep = "unused") %>%

  # extract ph domains
  mutate(ph_matches = str_match_all(domains,
                                    '(\\d+)\\.\\.(\\d+);\\s+/note="PH(?:\\s+\\d+)?"')) %>%

  # extract ph domain positions
  mutate(ph_vals = map(ph_matches, \(m) {

    # get start and end positions
    ph_starts <- if (!is.null(nrow(m)) && nrow(m) > 0) as.integer(m[, 2]) else integer(0)
    ph_ends   <- if (!is.null(nrow(m)) && nrow(m) > 0) as.integer(m[, 3]) else integer(0)

    # set a maximum of 5 ph domains (from manual inspection)
    ph_starts <- c(ph_starts, rep(NA_integer_, 5))[1:5]
    ph_ends   <- c(ph_ends,   rep(NA_integer_, 5))[1:5]

    # set names of ph columns
    setNames(as.list(c(rbind(ph_starts, ph_ends))),
             c(rbind(paste0("PHstart_", 1:5),
                     paste0("PHend_", 1:5))))})) %>%
  unnest_wider(ph_vals) %>%
  select(-ph_matches) %>%

  # get ph domain sequences
  mutate(PHseq_1 = str_sub(full_sequence, PHstart_1, PHend_1),
         PHseq_2 = str_sub(full_sequence, PHstart_2, PHend_2),
         PHseq_3 = str_sub(full_sequence, PHstart_3, PHend_3),
         PHseq_4 = str_sub(full_sequence, PHstart_4, PHend_4),
         PHseq_5 = str_sub(full_sequence, PHstart_5, PHend_5))

# making final ph dataframe
final_ph <- original_ph %>%
  select(uniprot_id, gene_name, full_sequence, contains("PH")) %>%
  pivot_longer(
    cols = matches("^PH(start|end|seq)_\\d+$"),
    names_to = c(".value", "PH_number"),
    names_pattern = "^PH(start|end|seq)_(\\d+)$"
  ) %>%
  rename(
    PH_start = start,
    PH_end = end,
    PH_seq = seq
  ) %>%
  mutate(PH_number = as.integer(PH_number)) %>%
  filter(!is.na(PH_start))

# making alphaPulldown dataframe
alphaPulldown_fullLength <- final_ph %>%
  distinct(uniprot_id) %>%
  mutate(uniprot_id = paste0(uniprot_id, "+", uniprot_id))

alphaPulldown_phOnly <- final_ph %>%
  select(uniprot_id, contains("PH")) %>%
  mutate(uniprot_id = paste0(uniprot_id, ":", PH_start, "-", PH_end)) %>%
  mutate(uniprot_id = paste0(uniprot_id, "+", uniprot_id)) %>%
  distinct(uniprot_id)

alphaPulldown_FLplusPH <- final_ph %>%
  select(uniprot_id, contains("PH")) %>%
  mutate(ph_id = paste0(uniprot_id, ":", PH_start, "-", PH_end)) %>%
  mutate(uniprot_id = paste0(uniprot_id, "+", ph_id)) %>%
  distinct(uniprot_id)

alphaPulldown_final <- bind_rows(alphaPulldown_fullLength,
                                 alphaPulldown_phOnly,
                                 alphaPulldown_FLplusPH)

write.csv(alphaPulldown_final,
          "alphaServer/phDimers/alphaPulldown_sampleSheet.csv",
          quote = F, row.names = F)

## reading in uniprot data for additional alphaPulldown run with ph-containing human proteins that were missed
# loading old uniprot ph spreadhseet
original_ph <- read_excel("alphaServer/phDimers/phProteins.xlsx")

# loading new uniprot ph spreadsheet
new_ph <- read_excel("alphaServer/phDimers/uniprotkb_taxonomy_id_9606_AND_reviewedDomains_2026_08_06.xlsx") %>%
  mutate(has_ph = grepl("note=\"PH", `Domain [FT]`)) %>%
  filter(has_ph == T) %>%
  filter(!(Entry %in% original_ph$uniprot_id)) %>%
  rename("uniprot_id" = "Entry",
         "uniprot_name" = "Entry Name",
         "gene_names" = "Gene Names",
         "length" = "Length",
         "domains" = "Domain [FT]",
         "full_sequence" = "Sequence") %>%
  select(-has_ph)

# processing new uniprot ph spreadsheet
new_ph_df <- new_ph %>%

  # focus on first gene name
  mutate(gene_name = gsub(" .*", "", gene_names), .keep = "unused") %>%

  # extract ph domains
  mutate(ph_matches = str_match_all(domains,
                                    '(\\d+)\\.\\.(\\d+);\\s+/note="PH(?:\\s+\\d+)?"')) %>%

  # extract ph domain positions
  mutate(ph_vals = map(ph_matches, \(m) {

    # get start and end positions
    ph_starts <- if (!is.null(nrow(m)) && nrow(m) > 0) as.integer(m[, 2]) else integer(0)
    ph_ends   <- if (!is.null(nrow(m)) && nrow(m) > 0) as.integer(m[, 3]) else integer(0)

    # set a maximum of 5 ph domains (from manual inspection)
    ph_starts <- c(ph_starts, rep(NA_integer_, 5))[1:5]
    ph_ends   <- c(ph_ends,   rep(NA_integer_, 5))[1:5]

    # set names of ph columns
    setNames(as.list(c(rbind(ph_starts, ph_ends))),
             c(rbind(paste0("PHstart_", 1:5),
                     paste0("PHend_", 1:5))))})) %>%
  unnest_wider(ph_vals) %>%
  select(-ph_matches) %>%

  # get ph domain sequences
  mutate(PHseq_1 = str_sub(full_sequence, PHstart_1, PHend_1),
         PHseq_2 = str_sub(full_sequence, PHstart_2, PHend_2),
         PHseq_3 = str_sub(full_sequence, PHstart_3, PHend_3),
         PHseq_4 = str_sub(full_sequence, PHstart_4, PHend_4),
         PHseq_5 = str_sub(full_sequence, PHstart_5, PHend_5))

# making final ph dataframe for new ph proteins
new_ph_df_final <- new_ph_df %>%
  select(uniprot_id, gene_name, full_sequence, contains("PH")) %>%
  pivot_longer(
    cols = matches("^PH(start|end|seq)_\\d+$"),
    names_to = c(".value", "PH_number"),
    names_pattern = "^PH(start|end|seq)_(\\d+)$"
  ) %>%
  rename(
    PH_start = start,
    PH_end = end,
    PH_seq = seq
  ) %>%
  mutate(PH_number = as.integer(PH_number)) %>%
  filter(!is.na(PH_start))

# making new alphaPulldown dataframe
alphaPulldown_fullLength_new <- new_ph_df_final %>%
  distinct(uniprot_id) %>%
  mutate(uniprot_id = paste0(uniprot_id, "+", uniprot_id))

alphaPulldown_phOnly_new <- new_ph_df_final %>%
  select(uniprot_id, contains("PH")) %>%
  mutate(uniprot_id = paste0(uniprot_id, ":", PH_start, "-", PH_end)) %>%
  mutate(uniprot_id = paste0(uniprot_id, "+", uniprot_id)) %>%
  distinct(uniprot_id)

alphaPulldown_FLplusPH_new <- new_ph_df_final %>%
  select(uniprot_id, contains("PH")) %>%
  mutate(ph_id = paste0(uniprot_id, ":", PH_start, "-", PH_end)) %>%
  mutate(uniprot_id = paste0(uniprot_id, "+", ph_id)) %>%
  distinct(uniprot_id)

alphaPulldown_final_new <- bind_rows(alphaPulldown_fullLength_new,
                                 alphaPulldown_phOnly_new,
                                 alphaPulldown_FLplusPH_new)

# saving new alphaPulldown dataframe
write.csv(alphaPulldown_final_new,
          "alphaServer/phDimers/alphaPulldown_new_sampleSheet.csv",
          quote = F, row.names = F)


## ----metricLoad-------------------------------------------------------------------------------------------------------------
# define predictions directory
predictions_dir <- "alphaServer/phDimers/temp_fold/predictions"

# define job subdirectories
job_subDirs <- list.dirs(predictions_dir, recursive = FALSE, full.names = TRUE)

# function to process one job subdirectory at a time
get_best_model <- function(dir_path) {

  # define file paths
  ranking_file <- file.path(dir_path, "ranking_scores.csv")
  interface_file <- file.path(dir_path, "interfaces.csv")

  # skip if either file is missing
  if (!file.exists(ranking_file) || !file.exists(interface_file)) {
    return(NULL)
  }

  # read files
  interface_df <- read_csv(interface_file, show_col_types = FALSE)
  ranking_df <- read_csv(ranking_file, show_col_types = FALSE) %>%
    mutate(model_used = paste0("seed-", seed, "_sample-", sample)) %>%
    filter(model_used %in% interface_df$model_used)

  # skip if no models in interface dataframe
  if (nrow(ranking_df) == 0) {
    return(NULL)
  }

  # find best model via ranking dataframe
  best_model <- ranking_df %>%
    filter(ranking_score == max(ranking_score))

  # extract matching row from interfaces.csv
  best_metrics <- interface_df %>%
    filter(model_used %in% best_model$model_used)

  return(best_metrics)
}

# extract and save prediction metrics for best models
metric_df <- map_dfr(job_subDirs, get_best_model)
write.csv(metric_df, "alphaServer/phDimers/best_models.csv", row.names = F)


## ----idInterface------------------------------------------------------------------------------------------------------------
## calculate distances between Cbeta carbons (except glycine which is Calpha) using a cutoff of < 8 angstroms
# from geometry.py in AlphaJudge GitHub

# function to read in .cif files
read_cif <- function(cif_file) {

  # read in cif file
  lines_cif <- readLines(cif_file, warn = FALSE)

  # find all loop_ lines
  loop_index <- which(trimws(lines_cif) == "loop_")
  if (length(loop_index) == 0) {
    stop("No loop_ blocks in cif file: ", cif_file)
  }

  # make empty vectors
  atom_headers <- NULL
  atom_data_start <- NULL
  atom_data_end <- NULL

  # find appropriate _loop lines
  for (i in seq_along(loop_index)) {
    read_start <- loop_index[i]
    read_end <- if (i < length(loop_index)) loop_index[i + 1] - 1 else length(lines_cif)
    read_block <- lines_cif[read_start:read_end]
    read_block_trim <- trimws(read_block)

    # retrieve header lines from loop block
    header_index <- which(startsWith(read_block_trim, "_atom_site."))
    if (length(header_index) == 0) next

    # find atom_site loop
    atom_headers <- read_block_trim[header_index]
    last_header_local <- max(header_index)
    atom_data_start_local <- last_header_local + 1

    # continue until blank line or #
    data_lines_local <- character()
    for (j in atom_data_start_local:length(read_block_trim)) {
      x <- read_block_trim[j]
      if (x == "" || x == "#") break
      data_lines_local <- c(data_lines_local, read_block[j])
    }

    # break if no atom data in file
    if (length(data_lines_local) == 0) {
      stop("no atom data in file: ", cif_file)
    }

    # retrieve atom data
    atom_data_start <- read_start + atom_data_start_local - 1
    atom_data_end <- atom_data_start + length(data_lines_local) - 1
    atom_data <- lines_cif[atom_data_start:atom_data_end]
    break
  }

  # break if no atom data in file
  if (is.null(atom_headers)) {
    stop("no _atom_site loop in file: ", cif_file)
  }

  # clean column names
  col_names <- sub("^_atom_site\\.", "", atom_headers)

  # make atom dataframe
  atom_df <- read.table(text = atom_data, header = FALSE,
                        stringsAsFactors = FALSE, fill = TRUE, quote = "",
                        comment.char = "")

  # break if error while making atom dataframe
  if (ncol(atom_df) != length(col_names)) {
    stop(
      "column count doesn't match header count: ",
      cif_file, "\n",
      "found ", ncol(atom_df), " columns but expected ", length(col_names)
    )
  }

  colnames(atom_df) <- col_names
  atom_df
}

# function to filter for relavent atoms (Cbeta for non-GLY, Calpha for GLY)
get_rel_atoms <- function(cif_file) {

  # get atomic dataframe from cif file
  atom_df <- read_cif(cif_file)

  # clean up atom dataframe
  atom_df_clean <- atom_df %>%
    filter(group_PDB == "ATOM") %>%
    mutate(atom_id = trimws(label_atom_id),
           residue_id = trimws(label_comp_id),
           chain_id = trimws(label_asym_id),
           sequence_id = as.integer(label_seq_id),
           cart_x = as.numeric(Cartn_x),
           cart_y = as.numeric(Cartn_y),
           cart_z = as.numeric(Cartn_z))

  # filter for relavent atoms
  rel_atoms <- atom_df_clean %>%
    filter((residue_id == "GLY" & atom_id == "CA") |
             (residue_id != "GLY" & atom_id == "CB")) %>%
    transmute(chain = chain_id, resNumber = sequence_id, resID = residue_id,
              atom = atom_id, x = cart_x, y = cart_y, z = cart_z) %>%
    arrange(chain, resNumber)

  rel_atoms
}

# function to generate a pairwise distance matrix
pairwise_distances <- function(df_1, df_2) {

  # getting relavent atom coordinates
  coords_1 <- as.matrix(df_1[, c("x", "y", "z")])
  coords_2 <- as.matrix(df_2[, c("x", "y", "z")])

  # making distance matrix
  distance_matrix <- outer(seq_len(nrow(coords_1)), seq_len(nrow(coords_2)),
                           Vectorize(function(i, j) {
                             sqrt(sum((coords_1[i, ] - coords_2[j, ])^2))}))

  distance_matrix
}

# function to find interacting residues
find_interactors <- function(cif_file, distance_cutoff = 8) {

  # get relavent atoms from cif file
  rel_atoms <- get_rel_atoms(cif_file)

  chains <- unique(rel_atoms$chain)

  # break if non-dimer
  if (length(chains) != 2) {
    stop("expected 2 chains in ", basename(cif_file),
      " but found: ", paste(chains, collapse = ", "))
  }

  # separate atom data by chain
  chain_1 <- chains[1]
  chain_2 <- chains[2]
  df_1 <- rel_atoms %>% filter(chain == chain_1)
  df_2 <- rel_atoms %>% filter(chain == chain_2)

  # obtain distance matrices and apply distance cutoff
  distance_matrix <- pairwise_distances(df_1, df_2)
  interface_hits <- which(distance_matrix <= distance_cutoff, arr.ind = TRUE)

  # generate interface table
  interface_table <- tibble(file = basename(cif_file),
                            chain_1 = chain_1,
                            resNumber_1 = df_1$resNumber[interface_hits[, 1]],
                            resID_1 = df_1$resID[interface_hits[, 1]],
                            atom_1  = df_1$atom[interface_hits[, 1]],
                            chain_2 = chain_2,
                            resNumber_2 = df_2$resNumber[interface_hits[, 2]],
                            resID_2 = df_2$resID[interface_hits[, 2]],
                            atom_2  = df_2$atom[interface_hits[, 2]],
                            distance = distance_matrix[interface_hits]) %>%
    arrange(distance, resNumber_1, resNumber_2)

  # identify interface residues
  interface_residues <- bind_rows(interface_table %>%
                                    transmute(file, chain = chain_1,
                                              resNumber = resNumber_1,
                                              resID = resID_1),
                                  interface_table %>%
                                    transmute(file, chain = chain_2,
                                              resNumber = resNumber_2,
                                              resID = resID_2)) %>%
    distinct() %>%
    arrange(chain, resNumber)

  list(file = basename(cif_file),
       rel_atoms = rel_atoms,
       interface_table = interface_table,
       interface_residues = interface_residues)
}

# function to find interface residues from .cif files in a directory
find_interfaces_inDir <- function(dir_path, file_pattern = "\\.cif$",
                                           distance_cutoff = 8) {

  # get list of .cif files from directory
  cif_files <- list.files(dir_path, pattern = file_pattern, full.names = TRUE)

  # break if no .cif files in direcotry
  if (length(cif_files) == 0) {
    stop("no .cif files in directory: ", dir_path)
  }

  # run functions
  interface_results <- map(cif_files, ~find_interactors(.x,
                                              distance_cutoff = distance_cutoff))

  # return function outputs
  interface_table <- map_dfr(interface_results, "interface_table")
  interface_residues <- map_dfr(interface_results, "interface_residues")

  list(files = cif_files, interface_results = interface_results,
       interface_table = interface_table, interface_residues = interface_residues)
}

# get and save all interface data
all_interfaces <- find_interfaces_inDir(dir_path =
                                          "alphaServer/phDimers/temp_fold/best_cifs",
                                        file_pattern = "\\.cif$",
                                        distance_cutoff = 8)
interface_table_all <- all_interfaces$interface_table
interface_residues_all <- all_interfaces$interface_residues
write.csv(interface_table_all, "alphaServer/phDimers/interface_table.csv",
          row.names = F)
write.csv(interface_residues_all, "alphaServer/phDimers/interface_residues.csv",
          row.names = F)


## ----phID-------------------------------------------------------------------------------------------------------------------
# loading data
metric_df_all <- read.csv("alphaServer/phDimers/best_models.csv") %>%
  select(-interface) %>%
  mutate(colon_count = str_count(jobs, ":")) %>%
  mutate(job_type = case_when(colon_count == 0 ~ "FL_FL",
                              colon_count == 1 ~ "FL_PH",
                              colon_count == 2 ~ "PH_PH"))
interface_residues_all <- read.csv("alphaServer/phDimers/interface_residues.csv") %>%
  mutate(hypen_count = str_count(file, "-"),
         and_count = str_count(file, "and"),
         protein = str_extract(file, "^[a-zA-Z0-9]*")) %>%
  mutate(job_type = case_when(hypen_count == 0 & and_count == 0 ~ "FL_FL",
                              hypen_count == 1 & and_count == 1 ~ "FL_PH",
                              hypen_count == 1 & and_count == 0 ~ "PH_PH"))
interface_table_all <- read.csv("alphaServer/phDimers/interface_table.csv")

# extracting ph domain residues
phDomains <- metric_df_all %>%
  mutate(protein = str_extract(jobs, "^\\w*"),
         phDomain = str_extract(jobs, "[0-9]*-[0-9]*")) %>%
  separate(phDomain, into = c("phStart", "phEnd")) %>%
  distinct(protein, phStart, phEnd) %>%
  filter(is.na(phStart) == F) %>%
  mutate(phStart = as.numeric(phStart),
         phEnd = as.numeric(phEnd))

# classifying ph domain residues
interface_residues_ph <- data.frame()
for (i in 1:nrow(interface_residues_all)) {

  # get residue of interest
  plot_row <- interface_residues_all %>%
    slice(i)

  # get protein of interest
  plot_protein <- plot_row$protein

  # find ph domain(s) of interest
  plot_ph <- phDomains %>%
    filter(protein == plot_protein)

  # make temporary dataframe to store ph test results
  ph_df <- data.frame(ph_domain = 1:nrow(plot_ph),
                      in_ph = rep(F, nrow(plot_ph)))

  # FL + FL -> ph if in any ph domain for either chain
  if (plot_row$job_type == "FL_FL") {
    for (j in 1:nrow(plot_ph)) {

      # get each ph domain
      plot_each_ph <- plot_ph %>%
        slice(j)

      # test for ph domain presence
      if (plot_row$resNumber >= plot_each_ph$phStart &
          plot_row$resNumber <= plot_each_ph$phEnd) {
        ph_df[j, 2] <- T
      }
    }
  }

  # FL + PH -> ph if in any ph domain for chain a or for any chain b residue
  if (plot_row$job_type == "FL_PH") {
    if (plot_row$chain == "A") {
      for (j in 1:nrow(plot_ph)) {

        # get each ph domain
        plot_each_ph <- plot_ph %>%
          slice(j)

        # test for ph domain presence
        if (plot_row$resNumber >= plot_each_ph$phStart &
            plot_row$resNumber <= plot_each_ph$phEnd) {
          ph_df[j, 2] <- T
        }
      }
    }
    if (plot_row$chain == "B") {
      ph_df$in_ph <- 1
    }
  }

  # PH + PH -> ph for all residues
  if (plot_row$job_type == "PH_PH") {
    ph_df$in_ph <- 1
  }

  # make dataframe and add to previous
  plot_df <- plot_row %>%
    mutate(in_ph = as.logical(max(ph_df$in_ph)))
  interface_residues_ph <- bind_rows(interface_residues_ph, plot_df)
}
write.csv(interface_residues_ph, "alphaServer/phDimers/interface_residues_ph.csv",
          row.names = F)


## ----pipID------------------------------------------------------------------------------------------------------------------
# reading in pymol sasa files for ph homodimers
sasa_files <- list.files(path = "alphaServer/phDimers/temp_fold/best_cifs/sasa_csv_outputs",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE)
sasa_df_original <- map_dfr(sasa_files, function(csv_file) {
  read_csv(csv_file,
           show_col_types = F,
           col_names = c("chain", "residue_number", "residue_id",
                         "secondary_structure", "sasa_monomer", "sasa_dimer"),
           skip = 1) %>%
    mutate(sasa_delta = sasa_dimer - sasa_monomer,
           sasa_plot = -sasa_delta,
           sample = basename(csv_file))}) %>%
  mutate(sample = gsub("__x2_model_residue_sasa.csv", "", sample)) %>%
  separate(sample, into = c("uniprot_id", "ph_start", "ph_end"), remove = F)

# identifying beta strands and connecting loops
sasa_df_strands <- sasa_df_original %>%
  group_by(sample, chain) %>%
  arrange(residue_number, .by_group = T) %>%
  mutate(is_strand = secondary_structure == "strand",
         run_id = cumsum(is_strand != lag(is_strand, default = first(is_strand))),
         strand_start = is_strand & !lag(is_strand, default = F),
         strand_number = cumsum(strand_start)) %>%
  group_by(sample, chain, run_id) %>%
  mutate(run_is_strand = first(is_strand),
         run_strand_number = if_else(run_is_strand,
                                     first(strand_number), NA_integer_)) %>%
  ungroup() %>%
  group_by(sample, chain) %>%
  mutate(previous_strand = if_else(is_strand, strand_number, strand_number),
         strand_region = if_else(is_strand,
                                 paste0("strand_", strand_number),
                                 paste0("connect_", strand_number, "_",
                                        strand_number + 1))) %>%
  ungroup() %>%
  select(-is_strand, -run_id, -strand_start, -strand_number, -run_is_strand,
         -run_strand_number,-previous_strand) %>%
  mutate(pip_moderate = case_when(strand_region == "connect_1_2" &
                                 residue_id == "ARG" ~ T,
                               strand_region == "connect_1_2" &
                                 residue_id == "LYS" ~ T,
                               strand_region == "strand_1" &
                                 residue_id == "ARG" ~ T,
                               strand_region == "strand_1" &
                                 residue_id == "LYS" ~ T,
                               strand_region == "strand_2" &
                                 residue_id == "ARG" ~ T,
                               strand_region == "strand_2" &
                                 residue_id == "LYS" ~ T,
                               T ~ F))

# save pip-binding residue results
write.csv(sasa_df_strands,
          "alphaServer/phDimers/sasa_strands.csv",
          row.names = F)


## ----databaseLoad-----------------------------------------------------------------------------------------------------------
## loading human genome databases
entrez_uniprot <- as.data.frame(org.Hs.egUNIPROT)
entrez_symbol <- as.data.frame(org.Hs.egSYMBOL)
entrez_combo <- full_join(entrez_uniprot, entrez_symbol, by = "gene_id") %>%
  filter(!(uniprot_id == "Q5VW22" & symbol == "AGAP4"))


## ----phAnalysis-------------------------------------------------------------------------------------------------------------
# loading data
metric_df_all <- read.csv("alphaServer/phDimers/best_models.csv") %>%
  select(-interface) %>%
  mutate(uniprot_id = str_extract(jobs, "^[a-zA-Z0-9]*")) %>%
  mutate(colon_count = str_count(jobs, ":")) %>%
  mutate(job_type = case_when(colon_count == 0 ~ "FL_FL",
                              colon_count == 1 ~ "FL_PH",
                              colon_count == 2 ~ "PH_PH")) %>%
  filter(confidence_score > 0) %>%
  left_join(entrez_combo %>% distinct(uniprot_id, symbol), by = "uniprot_id")
interface_residues_ph <- read.csv("alphaServer/phDimers/interface_residues_ph.csv")
interface_table_all <- read.csv("alphaServer/phDimers/interface_table.csv")
sasa_strands <- read.csv("alphaServer/phDimers/sasa_strands.csv") %>%
  
  # removing non-human gab1 ph domains from human proteome screen
  filter(uniprot_id != "Q9QYY0") %>% # mouse gab1
  filter(uniprot_id != "G5EDJ4") %>% # worm soc1
  filter(uniprot_id != "Q9VZZ9") # fly dos

# run t test on all complexes for difference in pip sasa
sasa_test <- sasa_strands %>%
  filter(pip_moderate == T) %>%
  group_by(sample) %>%
  summarize(sasa_p = t.test(sasa_monomer, sasa_dimer)$p.value) %>%
  ungroup() %>%
  mutate(sasa_q = p.adjust(sasa_p, method = "fdr"),
         job_id = gsub("_", ":", sample),
         jobs = paste0(job_id, "+", job_id),
         uniprot_id = gsub("_.*", "", sample))

# summarizing delta sasa for pip residues per complex
sasa_summary <- sasa_strands %>%
  filter(pip_moderate == T) %>%
  group_by(sample) %>%
  summarize(pip_sasa_monomer = sum(sasa_monomer),
            pip_sasa_dimer = sum(sasa_dimer)) %>%
  ungroup() %>%
  mutate(percent_sasa_dimer = pip_sasa_dimer / pip_sasa_monomer * 100,
         sasa_log2FC = log2(pip_sasa_dimer / pip_sasa_monomer),
         job_id = gsub("_", ":", sample),
         jobs = paste0(job_id, "+", job_id)) %>%
  left_join(sasa_test, by = c("jobs", "job_id", "sample"))

# getting ph domain information
phDomains <- metric_df_all %>%
  mutate(protein = str_extract(jobs, "^\\w*"),
         phDomain = str_extract(jobs, "[0-9]*-[0-9]*")) %>%
  separate(phDomain, into = c("phStart", "phEnd")) %>%
  distinct(protein, phStart, phEnd) %>%
  filter(is.na(phStart) == F) %>%
  mutate(phStart = as.numeric(phStart),
         phEnd = as.numeric(phEnd),
         phLength = phEnd - phStart + 1) %>%
  group_by(protein) %>%
  summarize(ph_residues = sum(phLength))

# summarzing ph residues per job
interface_ph_summary <- interface_residues_ph %>%
  group_by(file, protein, job_type) %>%
  summarize(interface_residues = n(),
            interface_ph = sum(in_ph),
            .groups = "drop_last") %>%
  left_join(phDomains, by = "protein") %>%
  mutate(percent_interface_with_ph = interface_ph / interface_residues * 100,
         percent_ph_in_interface = interface_ph / ph_residues * 100 / 2,
         jobs_1 = gsub("_model.cif", "", file),
         jobs_1 = gsub("__x2", "", jobs_1),
         jobs_2 = jobs_1,
         jobs_1 = gsub("_and.*", "", jobs_1),
         jobs_2 = gsub(".*_and_", "", jobs_2),
         jobs = paste0(jobs_1, "+", jobs_2),
         jobs = gsub("_", ":", jobs),
         phDomain = str_extract(jobs, "[0-9]*-[0-9]*")) %>%
  select(-contains("jobs_"))

# creating a dataframe with filtering results
interface_ph_summary_wide <- interface_ph_summary %>%
  ungroup() %>%
  select(jobs, protein, phDomain, job_type, percent_interface_with_ph) %>%
  pivot_wider(names_from = job_type, values_from = percent_interface_with_ph) %>%
  rename(uniprot_id = protein) %>%
  distinct(uniprot_id, phDomain, FL_FL, FL_PH)
interface_ph_summary_fl_fl <- interface_ph_summary_wide %>%
  distinct(uniprot_id, FL_FL) %>%
  filter(is.na(FL_FL) == F)
interface_ph_summary_fl_ph <- interface_ph_summary_wide %>%
  distinct(uniprot_id, phDomain, FL_PH) %>%
  filter(is.na(FL_PH) == F)

# creating a datframe with filtering cutoffs
metric_df_cutoffs <- metric_df_all %>%
  filter(job_type == "PH_PH") %>%
  mutate(phDomain = str_extract(jobs, "[0-9]*-[0-9]*")) %>%
  left_join(interface_ph_summary_fl_fl, by = "uniprot_id") %>%
  left_join(interface_ph_summary_fl_ph, by = c("uniprot_id", "phDomain")) %>%
  mutate(FL_FL = case_when(is.na(FL_FL) == T ~ 0,
                           T ~ FL_FL),
         FL_PH = case_when(is.na(FL_PH) == T ~ 0,
                           T ~ FL_PH))

# applying filter cutoffs
metric_df_upset <- metric_df_cutoffs %>%
  distinct(jobs, uniprot_id, phDomain, symbol, iptm, pDockQ.mpDockQ,
           interface_score, FL_FL, FL_PH) %>%
  mutate(pass_iptm = case_when(iptm > 0.3 ~ 1, T ~ 0),
         pass_dock = case_when(pDockQ.mpDockQ > 0.23 ~ 1, T ~ 0),
         pass_interfaceScore = case_when(interface_score > 130 ~ 1, T ~ 0),
         pass_flFLinterface = case_when(FL_FL > 10 ~ 1, T ~ 0),
         pass_flPHinterface = case_when(FL_PH > 50 ~ 1, T ~ 0))

## supplemental figure 7c
# applying filtered cutoffs to dataframe
metric_df_filtered <- metric_df_cutoffs %>%
  filter(job_type == "PH_PH") %>%
  filter(iptm > 0.3 & pDockQ.mpDockQ > 0.23 & interface_score > 130 &
           FL_FL > 10 & FL_PH > 50) %>%
  left_join(sasa_summary %>%
              select(jobs, contains("sasa")), by = "jobs")

# making dataframe of filter cutoffs for facet plot
plot_metric_cutoffs <- data.frame(filter_metric = c("ipTM", "pDockQ", "Interface Score", 
                               "FL-FL PH Interface", "FL-PH PH Interface"),
                               filter_cutoff = c(0.3, 0.23, 130, 10, 50))

## supplemental figure 7a
# plotting filter cutoffs
af_cutoffs_plot <- ggplot(metric_df_cutoffs %>%
         distinct(jobs, uniprot_id, phDomain, symbol, iptm, pDockQ.mpDockQ,
                  interface_score, FL_FL, FL_PH) %>%
         rename("ipTM" = "iptm",
                "pDockQ" = "pDockQ.mpDockQ",
                "Interface Score" = "interface_score",
                "FL-FL PH Interface" = "FL_FL",
                "FL-PH PH Interface" = "FL_PH") %>%
         pivot_longer(cols = c("ipTM", "pDockQ", "Interface Score", 
                               "FL-FL PH Interface", "FL-PH PH Interface"),
                      names_to = "filter_metric", values_to = "filter_value"),
       aes(x = filter_value)) +
  geom_histogram(bins = 50, fill = "black") +
  facet_wrap(~ factor(filter_metric, levels = c("ipTM", "pDockQ", 
                                                "Interface Score", 
                                                "FL-FL PH Interface", 
                                                "FL-PH PH Interface")), 
             scales = "free", nrow = 5) +
  geom_vline(data = plot_metric_cutoffs %>% filter(filter_metric == "ipTM"),
             aes(xintercept = filter_cutoff), 
             color = "red", linetype = "dashed", linewidth = 0.1) +
  geom_vline(data = plot_metric_cutoffs %>% filter(filter_metric == "pDockQ"),
             aes(xintercept = filter_cutoff), 
             color = "red", linetype = "dashed", linewidth = 0.1) +
  geom_vline(data = plot_metric_cutoffs %>% filter(filter_metric == "Interface Score"),
             aes(xintercept = filter_cutoff), 
             color = "red", linetype = "dashed", linewidth = 0.1) +
  geom_vline(data = plot_metric_cutoffs %>% filter(filter_metric == "FL-FL PH Interface"),
             aes(xintercept = filter_cutoff), 
             color = "red", linetype = "dashed", linewidth = 0.1) +
  geom_vline(data = plot_metric_cutoffs %>% filter(filter_metric == "FL-PH PH Interface"),
             aes(xintercept = filter_cutoff), 
             color = "red", linetype = "dashed", linewidth = 0.1) +
  theme_classic(base_size = 7) +
  labs(x = "Score", y = "Count") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); af_cutoffs_plot
ggsave("manuscript/non_maria_figures/phDimers/af_cutoffs_plot.pdf", af_cutoffs_plot,
       units = "in", height = 3.41, width = 2.22)

## figure 6b
# volcano plot with sasa results
sasa_volcano <- ggplot(metric_df_filtered %>%
         mutate(plot_label = case_when(sasa_q < 0.1 ~ symbol,
                                       symbol == "BTK" ~ "BTK",
                                       symbol == "PLEK2" ~ "PLEK2",
                                       symbol == "SWAP70" ~ "SWAP70",
                                       T ~ ""),
                pip_hiding = case_when(sasa_q < 0.1 ~ T,
                                       T ~ F)),
       aes(x = sasa_log2FC, y = -log10(sasa_q), 
           label = plot_label, color = pip_hiding)) +
  geom_point(size = 0.000000000000000000000000000000000000000000000000000001) +
  theme_classic(base_size = 7) +
  geom_text_repel(min.segment.length = 0.05, max.overlaps = 20,
                   size = 2, segment.size = 0.25) +
  scale_color_manual(values = c("gray50", "black")) +
  labs(x = "log2(PIP SASA in Dimer / PIP SASA in Monomer)", y = "-log10(q-value)") +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 1, linewidth = 0.1, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); sasa_volcano
ggsave("manuscript/non_maria_figures/phDimers/sasa_volcano.pdf", sasa_volcano,
       units = "in", height = 1, width = 2.75)


## ----pipBury----------------------------------------------------------------------------------------------------------------
## figure 4i
# plotting delta pip sasa for gab1 ph-ph complex
gab1PH_plot <- ggplot(sasa_strands %>%
         filter(sample == "Q13480_5-116") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "GAB1 PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); gab1PH_plot
ggsave("manuscript/non_maria_figures/phDimers/gab1_plot.pdf", gab1PH_plot,
       units = "in", height = 2.36, width = 2.06)

## figure 6c
# plotting delta pip sasa for btk ph-ph complex
btkPH_plot <- ggplot(sasa_strands %>%
         filter(sample == "Q06187_3-133") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "BTK PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); btkPH_plot
ggsave("manuscript/non_maria_figures/phDimers/btkPH_plot.pdf", btkPH_plot,
       units = "in", height = 1.55, width = 2.06/2)

## figure 6d
# plotting delta pip sasa for plek2 ph-ph complex
plek2PH_plot <- ggplot(sasa_strands %>%
         filter(sample == "Q9NYT0_4-104") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "PLEK2 PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); plek2PH_plot
ggsave("manuscript/non_maria_figures/phDimers/plek2PH_plot.pdf", plek2PH_plot,
       units = "in", height = 1.55, width = 2.06/2)

## figure 6e
# plotting delta pip sasa for swap70 ph-ph complex
swap70PH_plot <- ggplot(sasa_strands %>%
         filter(sample == "Q9UH65_210-306") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "SWAP70 PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); swap70PH_plot
ggsave("manuscript/non_maria_figures/phDimers/swap70PH_plot.pdf", swap70PH_plot,
       units = "in", height = 1.55, width = 2.06/2)

## figure 6f
# plotting delta pip sasa for rasa4b ph-ph complex
rasa4bPH_plot <- ggplot(sasa_strands %>%
         filter(sample == "C9J798_566-673") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "RASA4B PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); rasa4bPH_plot
ggsave("manuscript/non_maria_figures/phDimers/rasa4bPH_plot.pdf", rasa4bPH_plot,
       units = "in", height = 1.55, width = 2.06/2)

## figure 6g
# plotting delta pip sasa for sh2b2 ph-ph complex
sh2b2PH_plot <- ggplot(sasa_strands %>%
         filter(sample == "O14492_193-306") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "SH2B2 PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); sh2b2PH_plot
ggsave("manuscript/non_maria_figures/phDimers/sh2b2PH_plot.pdf", sh2b2PH_plot,
       units = "in", height = 1.55, width = 2.06/2)

## figure 6h
# plotting delta pip sasa for veph1 ph-ph complex
veph1PH_plot <- ggplot(sasa_strands %>%
         filter(sample == "Q14D04_716-819") %>%
         filter(pip_moderate == T) %>%
         rename(Monomer = sasa_monomer,
                Dimer = sasa_dimer) %>%
         pivot_longer(cols = c("Monomer", "Dimer"),
                      names_to = "sasa_state", values_to = "sasa_value") %>%
         mutate(sasa_state = fct_relevel(as.factor(sasa_state), 
                                         "Monomer", "Dimer")),
       aes(x = sasa_state, y = sasa_value)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Monomer", "Dimer")),
                     label = "p.signif", method = "t.test",
                     size = 2) + 
  theme_classic(base_size = 7) +
  labs(x = "State", y = "VEPH1 PH/PH PIP-Binding SASA (A2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")); veph1PH_plot
ggsave("manuscript/non_maria_figures/phDimers/veph1PH_plot.pdf", veph1PH_plot,
       units = "in", height = 1.55, width = 2.06/2)


## ----sessionInfo------------------------------------------------------------------------------------------------------------
sessionInfo()

