## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(ggpubr)
library(tidyverse)


## ----slbEGFRgab1Cor---------------------------------------------------------------------------------------------------------
# load in data
egfr_cor_files <- list.files(path = "maria_data/maria_figures/slb_gab1FL_egfrCorrelation/",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE)
egfr_cor_df_original <- map_dfr(egfr_cor_files, function(csv_file) {
  read_csv(csv_file,
           show_col_types = F,
           col_names = c("roi_number", "roi_area", "roi_mean", "roi_sd",
                         "roi_min", "roi_max"),
           skip = 1) %>%
    mutate(exp = basename(dirname(csv_file)),
      sample = basename(csv_file))}) %>%
  mutate(exp = gsub("2024.", "", exp),
         sample = gsub(".csv", "", sample)) %>%
  separate(sample, into = c("image", "condition", "channel", "rois"), remove = F) %>%
  unite("full_comp", "exp", "sample", remove = F) %>%
  mutate(roi_totalIntensity = roi_area*roi_mean,
         gab1_level = case_when(rois == "egfrROIs" ~ "Low",
                                rois == "gab1ROIs" ~ "High")) %>%
  mutate(gab1_level = fct_relevel(gab1_level, "Low", "High")) 

## figure 1d
# plotting per-condensate correlation between egfr and gab1
slb_egfr_gab1_correlation <- ggplot(egfr_cor_df_original %>% filter(channel == "EGFR"), 
       aes(x = gab1_level, y = roi_mean)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("Low", "High")), 
                     label = "p.signif", method = "t.test", # wilcox is same
                     size = 2) + 
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GAB1 Level", y = "EGFR Intensity"); slb_egfr_gab1_correlation
ggsave("maria_data/maria_figures/slb_egfr_gab1_correlation.pdf",
       slb_egfr_gab1_correlation, units = "in",
       width = 0.86, height = 2.16)


## ----cellEGFRgab1Cor--------------------------------------------------------------------------------------------------------
# load in baseline data
cell_egfr_gab1_correlation_baseline <- read.csv("maria_data/maria_figures/cell_egfr_gab1_correlation_baseline.csv") %>%
  rename("Medium" = "Med") %>%
  pivot_longer(cols = c("Low", "Medium", "High"), 
               names_to = "gab1_level", values_to = "egfr_intensity") %>%
  mutate(gab1_level = fct_relevel(as.factor(gab1_level), "Low", "Medium", "High"))

## figure 1g
# plotting baseline data
cell_egfr_gab1_correlation_baseline_plot <- ggplot(cell_egfr_gab1_correlation_baseline, 
       aes(x = gab1_level, y = egfr_intensity)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("High", "Medium"),
                                        c("High", "Low"),
                                        c("Medium", "Low")), 
                     label = "p.signif", method = "t.test", # wilcox is same
                     size = 2) + 
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GAB1 Level", y = "EGFR Intensity"); cell_egfr_gab1_correlation_baseline_plot
ggsave("maria_data/maria_figures/cell_egfr_gab1_correlation_baseline.pdf",
       cell_egfr_gab1_correlation_baseline_plot, units = "in",
       width = 1.20, height = 2.16)


# load in serum-free data
cell_egfr_gab1_correlation_serumFree <- read.csv("maria_data/maria_figures/cell_egfr_gab1_correlation_serumFree.csv") %>%
  rename("Medium" = "Med") %>%
  pivot_longer(cols = c("Low", "Medium", "High"), 
               names_to = "gab1_level", values_to = "egfr_intensity") %>%
  filter(is.na(egfr_intensity) == F) %>%
  mutate(gab1_level = fct_relevel(as.factor(gab1_level), "Low", "Medium", "High"))

## supplemental figure 1i: left
# plotting serumFree data
cell_egfr_gab1_correlation_serumFree_plot <- ggplot(cell_egfr_gab1_correlation_serumFree, 
       aes(x = gab1_level, y = egfr_intensity)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("High", "Medium"),
                                        c("High", "Low"),
                                        c("Medium", "Low")), 
                     label = "p.signif", method = "t.test", # wilcox is same
                     size = 2) + 
  coord_cartesian(ylim = c(100, 1400)) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "GAB1 Level", y = "EGFR Intensity",
       title = "Serum-Starved"); cell_egfr_gab1_correlation_serumFree_plot
ggsave("maria_data/maria_figures/cell_egfr_gab1_correlation_serumFree.pdf",
       cell_egfr_gab1_correlation_serumFree_plot, units = "in",
       width = 2.36, height = 2.16)

# load in 5 minutes post-egf data
cell_egfr_gab1_correlation_EGF5min <- read.csv("maria_data/maria_figures/cell_egfr_gab1_correlation_EGF5min.csv") %>%
  rename("Medium" = "Med") %>%
  pivot_longer(cols = c("Low", "Medium", "High"), 
               names_to = "gab1_level", values_to = "egfr_intensity") %>%
  filter(is.na(egfr_intensity) == F) %>%
  mutate(gab1_level = fct_relevel(as.factor(gab1_level), "Low", "Medium", "High"))

## supplemental figure 1i: middle
# plotting EGF5min data
cell_egfr_gab1_correlation_EGF5min_plot <- ggplot(cell_egfr_gab1_correlation_EGF5min, 
       aes(x = gab1_level, y = egfr_intensity)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("High", "Medium"),
                                        c("High", "Low"),
                                        c("Medium", "Low")), 
                     label = "p.signif", method = "t.test", # wilcox is same
                     size = 2) + 
  coord_cartesian(ylim = c(100, 1400)) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "GAB1 Level", y = "EGFR Intensity",
       title = "5 Minutes Post-EGF"); cell_egfr_gab1_correlation_EGF5min_plot
ggsave("maria_data/maria_figures/cell_egfr_gab1_correlation_EGF5min.pdf",
       cell_egfr_gab1_correlation_EGF5min_plot, units = "in",
       width = 2.36, height = 2.16)


# load in 30 minutes post-egf data
cell_egfr_gab1_correlation_EGF30min <- read.csv("maria_data/maria_figures/cell_egfr_gab1_correlation_EGF30min.csv") %>%
  rename("Medium" = "Med") %>%
  pivot_longer(cols = c("Low", "Medium", "High"), 
               names_to = "gab1_level", values_to = "egfr_intensity") %>%
  filter(is.na(egfr_intensity) == F) %>%
  mutate(gab1_level = fct_relevel(as.factor(gab1_level), "Low", "Medium", "High"))

## supplemental figure 1i: right
# plotting EGF30min data
cell_egfr_gab1_correlation_EGF30min_plot <- ggplot(cell_egfr_gab1_correlation_EGF30min, 
       aes(x = gab1_level, y = egfr_intensity)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("High", "Medium"),
                                        c("High", "Low"),
                                        c("Medium", "Low")), 
                     label = "p.signif", method = "t.test", # wilcox is same
                     size = 2) + 
  coord_cartesian(ylim = c(100, 1400)) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5)) +
  labs(x = "GAB1 Level", y = "EGFR Intensity",
       title = "30 Minutes Post-EGF"); cell_egfr_gab1_correlation_EGF30min_plot
ggsave("maria_data/maria_figures/cell_egfr_gab1_correlation_EGF30min.pdf",
       cell_egfr_gab1_correlation_EGF30min_plot, units = "in",
       width = 2.36, height = 2.16)


## ----gab1DomainsSLBs--------------------------------------------------------------------------------------------------------
# load in data
slb_domain_files <- list.files(path = "maria_data/maria_figures/slb_gab1Domains/",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE)
slb_domain_df_original <- map_dfr(slb_domain_files, function(csv_file) {
  read_csv(csv_file,
           show_col_types = F,
           col_names = c("roi_number", "roi_area", "roi_mean", "roi_sd",
                         "roi_min", "roi_max", "roi_dir", "roi_file"),
           skip = 1) %>%
    mutate(exp = basename(dirname(csv_file)),
      sample = basename(csv_file))}) %>%
  mutate(exp = gsub("2024.", "", exp),
         sample = gsub("_results.csv", "", sample)) %>%
  separate(sample, into = c("image", "condition", "channel"), remove = F) %>%
  unite("full_comp", "exp", "sample", remove = F) %>%
  mutate(roi_totalIntensity = roi_area*roi_mean)

# load in metadata
slb_domain_df <- slb_domain_df_original %>%
  left_join(read_excel("maria_data/maria_figures/slb_gab1Domains/slb_gab1Domains_metadata.xlsx") %>%
              mutate(exp = gsub("delete.", "", exp),
                     image = gsub("delete.", "", image)),
            by = c("exp", "image", "condition", "channel")) %>%
  mutate(gab1 = case_when(gab1 == "IDR" ~ "IDR493", T ~ gab1),
         gab1 = case_when(condition == "postGrb2" ~ "preGAB1", T ~ gab1))

# isolate whole-fov measurements
slb_domain_fov <- slb_domain_df %>%
  filter(roi_number == 1) %>%
  distinct(full_comp, exp, image, condition, channel, gab1, roi_area, roi_mean, 
           roi_totalIntensity) %>%
  rename(fov_area = roi_area,
         fov_mean = roi_mean,
         fov_intensity = roi_totalIntensity)
  
# summarizing condensate metrics for each fov
slb_domain_condense <- slb_domain_df %>%
  filter(roi_number != 1) %>%
  group_by(full_comp, exp, image, condition, channel, gab1) %>%
  summarize(condense_area = sum(roi_area),
            condense_intensity = sum(roi_totalIntensity),
            condense_count = max(roi_number),
            .groups = "drop") %>%
  mutate(condense_mean = condense_intensity / condense_area,
         condense_size = condense_area / condense_count) %>%
  full_join(slb_domain_fov, 
            by = c("full_comp", "exp", "image", "condition", "channel", "gab1")) %>%
  mutate(condense_area = case_when(is.na(condense_area) == T ~ 0, 
                                   T ~ condense_area),
         condense_intensity = case_when(is.na(condense_intensity) == T ~ 0, 
                                        T ~ condense_intensity),
         condense_mean = case_when(is.na(condense_mean) == T ~ 0,
                                   T ~ condense_mean),
         condense_count = case_when(is.na(condense_count) == T ~ 0,
                                    T ~ condense_count),
         condense_size = case_when(is.na(condense_size) == T ~ 0,
                                   T ~ condense_size))

# plotting condensate metrics
ggplot(slb_domain_condense, 
       aes(x = gab1, y = condense_area)) +
  geom_boxplot() +
  facet_wrap(~ channel + gab1, scales = "free") +
  theme_classic()

## potential postGab1 outliers
# EGFR dPH: 10.09_037_postGab1_EGFR (fine)
# EGFR FL: 10.09_034_postGab1_EGFR (fine), 09.24_026_postGab1_EGFR (outlier)
# EGFR IDR: 09.24_032_postGab1_EGFR (fine)
# EGFR PRR: 10.09_044_postGab1_EGFR (fine)
# EGFR NONE: 10.09_047_postGab1_EGFR (fine)
# GAB1 dPH: 11.13_043_postGab1_Gab1 (outlier), 11.13_042_postGab1_Gab1 (outlier)
# GAB1 FL: 10.09_034_postGab1_Gab1 (fine), 11.22_007_postGab1_Gab1 (fine)
# GAB1 IDR: 11.13_044_postGab1_Gab1 (outlier)

# EGFR postGrb2: 11.22_001_postGrb2_EGFR (outlier), 
# 11.22_003_postGrb2_EGFR (outlier), 11.22_002_postGrb2_EGFR (outlier), 
# 11.22_011_postGrb2_EGFR (fine), 11.22_007_postGrb2_EGFR (fine)

# removing outliers and adding pseudocount
slb_domain_condense_final <- slb_domain_condense %>%
  filter(full_comp != "09.24_026_postGab1_EGFR") %>%
  filter(full_comp != "11.13_043_postGab1_Gab1") %>%
  filter(full_comp != "11.13_042_postGab1_Gab1") %>%
  filter(full_comp != "11.13_044_postGab1_Gab1") %>%
  filter(full_comp != "11.22_001_postGrb2_EGFR") %>%
  filter(full_comp != "11.22_003_postGrb2_EGFR") %>%
  filter(full_comp != "11.22_002_postGrb2_EGFR") %>%
  mutate(gab1 = case_when(gab1 == "preGAB1" ~ "Pre GAB1",
                          gab1 == "FL" ~ "GAB1-FL",
                          gab1 == "dPH" ~ "GAB1-dPH",
                          gab1 == "PRR" ~ "GAB1-PRR",
                          gab1 == "IDR493" ~ "GAB1-IDR493",
                          gab1 == "NONE" ~ "No GAB1")) %>%
  mutate(gab1 = fct_relevel(as.factor(gab1), "Pre GAB1", "GAB1-FL", "GAB1-dPH", 
                            "GAB1-PRR", "GAB1-IDR493", "No GAB1"),
         dilute_area = fov_area - condense_area,
         dilute_intensity = fov_intensity - condense_intensity,
         dilute_mean = dilute_intensity / dilute_area,
         partition = condense_mean / dilute_mean,
         condense_area_pseudo = condense_area + 0.1,
         condense_count_pseudo = condense_count + 1)

# plotting egfr condensate metrics without outliers
ggplot(slb_domain_condense_final %>% filter(channel == "EGFR"), 
       aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(ref.group = "Pre GAB1", label = "p.signif") +
  facet_wrap(~ channel, scales = "free") +
  theme_classic()

# plotting gab1 condensate metrics without outliers
ggplot(slb_domain_condense_final %>% filter(channel == "Gab1"), 
       aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(ref.group = "GAB1-FL", label = "p.format") +
  facet_wrap(~ channel, scales = "free") +
  theme_classic()

# mean condensate areas
slb_domain_condense_means <- slb_domain_condense_final %>%
  group_by(channel, gab1) %>%
  summarize(mean_pseudoArea = mean(condense_area_pseudo), .groups = "drop_last")

# finding representative images
slb_domain_condense_rep <- slb_domain_condense_final %>%
  left_join(slb_domain_condense_means, by = c("channel", "gab1")) %>%
  mutate(percent_diff = (abs(condense_area_pseudo - mean_pseudoArea)) / mean_pseudoArea * 100,
         condense_area_pseudo_norm = case_when(channel == "EGFR" ~ 
                                                 condense_area_pseudo * 100 / 
                                                 slb_domain_condense_means %>%
                                                 filter(channel == "EGFR" & gab1 == "Pre GAB1") %>%
                                                 pull(mean_pseudoArea),
                                               channel == "Gab1" ~ 
                                                 condense_area_pseudo * 100 /
                                                 slb_domain_condense_means %>%
                                                 filter(channel == "Gab1" & gab1 == "GAB1-FL") %>%
                                                 pull(mean_pseudoArea)))

## representative images
# EGFR Pre: 09.24_014_postGrb2_EGFR
# EGFR FL: 09.24_024_postGab1_EGFR
# EGFR dPH: 09.24_029_postGab1_EGFR
# EGFR PRR: 10.09_045_postGab1_EGFR
# EGFR IDR: 10.09_040_postGab1_EGFR
# EGFR No: 09.24_037_postGab1_EGFR
# GAB1 FL: 09.24_024_postGab1_Gab1
# GAB1 dPH: 09.24_029_postGab1_Gab1
# GAB1 PRR: 10.09_045_postGab1_Gab1
# GAB1 IDR: 10.09_040_postGab1_Gab1
# GAB1 No: 09.24_037_postGab1_Gab1

## figure 2b: left
egfr_condense_area_domains <- ggplot(slb_domain_condense_final %>% 
                                        filter(channel == "EGFR") %>%
                                       filter(gab1 != "Pre GAB1") %>%
                                       filter(gab1 != "No GAB1") %>%
                                        mutate(gab1 = case_when(gab1 == "GAB1-FL" ~ "GAB1-WT",
                                                                T ~ gab1),
                                               gab1 = fct_relevel(as.factor(gab1), "GAB1-WT",
                                                                  "GAB1-dPH",
                                                                  "GAB1-PRR",
                                                                  "GAB1-IDR493")), 
                                      aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(comparisons = list(c("GAB1-WT", "GAB1-dPH"),
                                        c("GAB1-WT", "GAB1-PRR"),
                                        c("GAB1-WT", "GAB1-IDR493"),
                                        c("GAB1-dPH", "GAB1-PRR"),
                                        c("GAB1-dPH", "GAB1-IDR493"),
                                        c("GAB1-PRR", "GAB1-IDR493")),
                     label = "p.signif", method = "wilcox",
                     size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "Total EGFR Condensate Area (um2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); egfr_condense_area_domains
ggsave("maria_data/maria_figures/slb_gab1Domains_egfr.pdf",
       egfr_condense_area_domains, units = "in",
       width = 1.46, height = 2.15)

## figure 2b: right
gab1_condense_area_domains <- ggplot(slb_domain_condense_final %>% 
                                        filter(channel == "Gab1") %>%
                                       filter(gab1 != "Pre GAB1") %>%
                                       filter(gab1 != "No GAB1") %>%
                                        mutate(gab1 = case_when(gab1 == "GAB1-FL" ~ "GAB1-WT",
                                                                T ~ gab1),
                                               gab1 = fct_relevel(as.factor(gab1), "GAB1-WT",
                                                                  "GAB1-dPH",
                                                                  "GAB1-PRR",
                                                                  "GAB1-IDR493")), 
                                      aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(comparisons = list(c("GAB1-WT", "GAB1-dPH"),
                                        c("GAB1-WT", "GAB1-PRR"),
                                        c("GAB1-WT", "GAB1-IDR493"),
                                        c("GAB1-dPH", "GAB1-PRR"),
                                        c("GAB1-dPH", "GAB1-IDR493"),
                                        c("GAB1-PRR", "GAB1-IDR493")),
                     label = "p.signif", method = "wilcox",
                     size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "Total GAB1 Condensate Area (um2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); gab1_condense_area_domains
ggsave("maria_data/maria_figures/slb_gab1Domains_gab1.pdf",
       gab1_condense_area_domains, units = "in",
       width = 1.46, height = 2.15)


## ----gab1ResidueSLBs--------------------------------------------------------------------------------------------------------
# load in data
slb_residue_files <- list.files(path = "maria_data/maria_figures/slb_gab1Residues/",
  pattern = "\\.csv$",
  full.names = TRUE,
  recursive = TRUE)
slb_residue_df_original <- map_dfr(slb_residue_files, function(csv_file) {
  read_csv(csv_file,
           show_col_types = F,
           col_names = c("roi_number", "roi_area", "roi_mean", "roi_sd",
                         "roi_min", "roi_max", "roi_dir", "roi_file"),
           skip = 1) %>%
    mutate(sample = basename(csv_file))}) %>%
  mutate(sample = gsub("_results.csv", "", sample)) %>%
  separate(sample, into = c("exp", "image", "condition", "channel"), 
           sep = "_") %>%
  unite("sample", "image", "condition", "channel", remove = F) %>%
  unite("full_comp", "exp", "sample", remove = F) %>%
  mutate(roi_totalIntensity = roi_area*roi_mean) %>%
  filter(condition != "postGrb2")

# load in metadata
slb_residue_df <- slb_residue_df_original %>%
  inner_join(read_excel("maria_data/maria_figures/slb_gab1Residues/slb_gab1Residues_metadata.xlsx") %>%
              mutate(image = gsub("delete.", "", image)),
            by = c("exp", "image", "condition", "channel")) %>%
  mutate(gab1 = case_when(gab1 == "FL" ~ "GAB1-FL", 
                          gab1 == "p" ~ "pGAB1",
                          gab1 == "E" ~ "GAB1-E",
                          gab1 == "F" ~ "GAB1-F",
                          T ~ gab1))

# isolate whole-fov measurements
slb_residue_fov <- slb_residue_df %>%
  filter(roi_number == 1) %>%
  distinct(full_comp, exp, image, condition, channel, gab1, roi_area, roi_mean, 
           roi_totalIntensity) %>%
  rename(fov_area = roi_area,
         fov_mean = roi_mean,
         fov_intensity = roi_totalIntensity)
  
# summarizing condensate metrics for each fov
slb_residue_condense <- slb_residue_df %>%
  filter(roi_number != 1) %>%
  group_by(full_comp, exp, image, condition, channel, gab1) %>%
  summarize(condense_area = sum(roi_area),
            condense_intensity = sum(roi_totalIntensity),
            condense_count = max(roi_number),
            .groups = "drop") %>%
  mutate(condense_mean = condense_intensity / condense_area,
         condense_size = condense_area / condense_count) %>%
  full_join(slb_residue_fov, 
            by = c("full_comp", "exp", "image", "condition", "channel", "gab1")) %>%
  mutate(condense_area = case_when(is.na(condense_area) == T ~ 0, 
                                   T ~ condense_area),
         condense_intensity = case_when(is.na(condense_intensity) == T ~ 0, 
                                        T ~ condense_intensity),
         condense_mean = case_when(is.na(condense_mean) == T ~ 0,
                                   T ~ condense_mean),
         condense_count = case_when(is.na(condense_count) == T ~ 0,
                                    T ~ condense_count),
         condense_size = case_when(is.na(condense_size) == T ~ 0,
                                   T ~ condense_size))

# plotting condensate metrics
ggplot(slb_residue_condense, 
       aes(x = gab1, y = condense_area)) +
  geom_boxplot() +
  facet_wrap(~ channel + gab1, scales = "free") +
  theme_classic()

## potential postGab1 outliers
# EGFR E: 2025.03.31_030_postGab1_EGFR (fine)
# EGFR FL: 2025.03.31_024_postGab1_EGFR (outlier), 2025.03.31_026_postGab1_EGFR (outlier)
# GAB1 E: ALL SHOULD BE ZERO!!
# GAB1 F: ALL SHOULD BE ZERO!!
# GAB1 FL: 2025.03.31_024_postGab1_GAB1 (outlier)

# removing outliers and adding pseudocount
slb_residue_condense_final <- slb_residue_condense %>%
  filter(full_comp != "2025.03.31_024_postGab1_EGFR") %>%
  filter(full_comp != "2025.03.31_026_postGab1_EGFR") %>%
  filter(full_comp != "2025.03.31_024_postGab1_GAB1") %>%
  mutate(condense_area = case_when((channel == "GAB1" & gab1 == "GAB1-E") ~ 0,
                                   (channel == "GAB1" & gab1 == "GAB1-F") ~ 0,
                                   T ~ condense_area),
         condense_intensity = case_when((channel == "GAB1" & gab1 == "GAB1-E") ~ 0,
                                   (channel == "GAB1" & gab1 == "GAB1-F") ~ 0,
                                   T ~ condense_intensity),
         condense_mean = case_when((channel == "GAB1" & gab1 == "GAB1-E") ~ 0,
                                   (channel == "GAB1" & gab1 == "GAB1-F") ~ 0,
                                   T ~ condense_mean),
         condense_count = case_when((channel == "GAB1" & gab1 == "GAB1-E") ~ 0,
                                   (channel == "GAB1" & gab1 == "GAB1-F") ~ 0,
                                   T ~ condense_count),
         condense_size = case_when((channel == "GAB1" & gab1 == "GAB1-E") ~ 0,
                                   (channel == "GAB1" & gab1 == "GAB1-F") ~ 0,
                                   T ~ condense_size)) %>%
  mutate(gab1 = fct_relevel(as.factor(gab1), "GAB1-FL", "pGAB1", 
                            "GAB1-E", "GAB1-F"),
         dilute_area = fov_area - condense_area,
         dilute_intensity = fov_intensity - condense_intensity,
         dilute_mean = dilute_intensity / dilute_area,
         partition = condense_mean / dilute_mean,
         condense_area_pseudo = condense_area + 0.1,
         condense_count_pseudo = condense_count + 1)

# plotting egfr condensate metrics without outliers
ggplot(slb_residue_condense_final %>% filter(channel == "EGFR"), 
       aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(ref.group = "GAB1-FL", label = "p.format") +
  facet_wrap(~ channel, scales = "free") +
  theme_classic()

# plotting gab1 condensate metrics without outliers
ggplot(slb_residue_condense_final %>% filter(channel == "GAB1"), 
       aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot() +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(ref.group = "GAB1-FL", label = "p.format") +
  facet_wrap(~ channel, scales = "free") +
  theme_classic()

# mean condensate areas
slb_residue_condense_means <- slb_residue_condense_final %>%
  group_by(channel, gab1) %>%
  summarize(mean_pseudoArea = mean(condense_area_pseudo), .groups = "drop_last")

# finding representative images
slb_residue_condense_rep <- slb_residue_condense_final %>%
  left_join(slb_residue_condense_means, by = c("channel", "gab1")) %>%
  mutate(percent_diff = (abs(condense_area_pseudo - mean_pseudoArea)) / mean_pseudoArea * 100,
         condense_area_pseudo_norm = case_when(channel == "EGFR" ~ 
                                                 condense_area_pseudo * 100 / 
                                                 slb_residue_condense_means %>%
                                                 filter(channel == "EGFR" & gab1 == "GAB1-FL") %>%
                                                 pull(mean_pseudoArea),
                                               channel == "GAB1" ~ 
                                                 condense_area_pseudo * 100 /
                                                 slb_residue_condense_means %>%
                                                 filter(channel == "GAB1" & gab1 == "GAB1-FL") %>%
                                                 pull(mean_pseudoArea)))

## representative images
# EGFR FL: 2025.03.10_014_postGab1_EGFR
# EGFR p: 2025.03.31_029_postGab1_EGFR
# EGFR E: 2025.03.31_032_postGab1_EGFR
# EGFR F: 2025.03.31_033_postGab1_EGFR
# GAB1 FL: 2025.03.31_025_postGab1_GAB1
# GAB1 p: 2025.03.10_017_postGab1_GAB1
# GAB1 E: 2025.03.10_022_postGab1_GAB1
# GAB1 F: 2025.03.10_025_postGab1_GAB1

## figure 2d: left
egfr_condense_area_residues <- ggplot(slb_residue_condense_final %>% 
                                        filter(channel == "EGFR") %>%
                                        mutate(gab1 = case_when(gab1 == "GAB1-FL" ~ "GAB1-WT",
                                                                gab1 == "pGAB1" ~ "GAB1-pY",
                                                                T ~ gab1),
                                               gab1 = fct_relevel(as.factor(gab1), "GAB1-WT",
                                                                  "GAB1-pY",
                                                                  "GAB1-E",
                                                                  "GAB1-F")), 
                                      aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(comparisons = list(c("GAB1-WT", "GAB1-pY"),
                                        c("GAB1-WT", "GAB1-E"),
                                        c("GAB1-WT", "GAB1-F"),
                                        c("GAB1-pY", "GAB1-E"),
                                        c("GAB1-pY", "GAB1-F"),
                                        c("GAB1-E", "GAB1-F")),
                     label = "p.signif", method = "wilcox",
                     size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "Total EGFR Condensate Area (um2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); egfr_condense_area_residues
ggsave("maria_data/maria_figures/slb_gab1Residues_egfr.pdf",
       egfr_condense_area_residues, units = "in",
       width = 1.46, height = 2.15)

## figure 2d: right
gab1_condense_area_residues <- ggplot(slb_residue_condense_final %>% 
                               filter(channel == "GAB1") %>%
                                        mutate(gab1 = case_when(gab1 == "GAB1-FL" ~ "GAB1-WT",
                                                                gab1 == "pGAB1" ~ "GAB1-pY",
                                                                T ~ gab1),
                                               gab1 = fct_relevel(as.factor(gab1), "GAB1-WT",
                                                                  "GAB1-pY",
                                                                  "GAB1-E",
                                                                  "GAB1-F")), 
       aes(x = gab1, y = condense_area_pseudo)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  scale_y_continuous(trans = "log10") +
  stat_compare_means(comparisons = list(c("GAB1-WT", "GAB1-pY"),
                                        c("GAB1-WT", "GAB1-E"),
                                        c("GAB1-WT", "GAB1-F"),
                                        c("GAB1-pY", "GAB1-E"),
                                        c("GAB1-pY", "GAB1-F"),
                                        c("GAB1-E", "GAB1-F")),
                     label = "p.signif", method = "wilcox",
                     size = 2) +
  theme_classic(base_size = 7) +
  labs(x = "", y = "Total GAB1 Condensate Area (um2)") +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); gab1_condense_area_residues
ggsave("maria_data/maria_figures/slb_gab1Residues_gab1.pdf",
       gab1_condense_area_residues, units = "in",
       width = 1.46, height = 2.15)


## ----slbScanLines-----------------------------------------------------------------------------------------------------------
# read in pre-gab1 line scan data
linescan_pre <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_preGAB1.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("GAB")) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup()

## supplemental figure 2a: first row, first column
# plotting pre-gab1 line scan data
plot_linescan_pre <- ggplot(linescan_pre, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 0 uM GAB1"); plot_linescan_pre
ggsave("maria_data/maria_figures/lineScan_slb_preGAB1.pdf",
       plot_linescan_pre, units = "in",
       width = 1.49, height = 1.72)


# read in gab1-wt line scan data
linescan_wt <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-WT.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "Gab1.FL.0" = "Gab1.FL") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("Gab1.FL", "GAB1-WT", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "GAB1-WT"))

## supplemental figure 2a: second row, first column
# plotting gab1-wt line scan data
plot_linescan_wt <- ggplot(linescan_wt, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.1, 0.9)) +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-WT"); plot_linescan_wt
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-WT.pdf",
       plot_linescan_wt, units = "in",
       width = 1.49, height = 1.72)


# read in pGab1 line scan data
linescan_pY <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-pY.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "pGab1.0" = "pGab1") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("pGab1", "pGAB1", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "pGAB1"))

## supplemental figure 2a: first row, second column
# plotting pGab1 line scan data
plot_linescan_pY <- ggplot(linescan_pY, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-pY"); plot_linescan_pY
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-pY.pdf",
       plot_linescan_pY, units = "in",
       width = 1.49, height = 1.72)


# read in gab1-e line scan data
linescan_e <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-E.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "Gab1.E.0" = "Gab1.E") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  filter(is.na(microns) == F) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("Gab1.E", "GAB1-E", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "GAB1-E"))

## supplemental figure 2a: first row, third column
# plotting gab1-e line scan data
plot_linescan_e <- ggplot(linescan_e, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-E"); plot_linescan_e
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-E.pdf",
       plot_linescan_e, units = "in",
       width = 1.49, height = 1.72)


# read in gab1-f line scan data
linescan_f <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-F.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "Gab1.F.0" = "Gab1.F") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  filter(is.na(microns) == F) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("Gab1.F", "GAB1-F", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "GAB1-F"))

## supplemental figure 2a: first row, fourth column
# plotting gab1-f line scan data
plot_linescan_f <- ggplot(linescan_f, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-F"); plot_linescan_f
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-F.pdf",
       plot_linescan_f, units = "in",
       width = 1.49, height = 1.72)


# read in gab1-dPH line scan data
linescan_dPH <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-dPH.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "Gab1.ΔPH.0" = "Gab1.ΔPH") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  filter(is.na(microns) == F) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("Gab1.ΔPH", "GAB1-dPH", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "GAB1-dPH"))

## supplemental figure 2a: second row, second column
# plotting gab1-dPH line scan data
plot_linescan_dPH <- ggplot(linescan_dPH, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-dPH"); plot_linescan_dPH
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-dPH.pdf",
       plot_linescan_dPH, units = "in",
       width = 1.49, height = 1.72)


# read in gab1-prr line scan data
linescan_prr <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-PRR.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "Gab1.PRR.0" = "Gab1.PRR") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  filter(is.na(microns) == F) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("Gab1.PRR", "GAB1-PRR", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "GAB1-PRR"))

## supplemental figure 2a: second row, third column
# plotting gab1-prr line scan data
plot_linescan_prr <- ggplot(linescan_prr, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-PRR"); plot_linescan_prr
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-PRR.pdf",
       plot_linescan_prr, units = "in",
       width = 1.49, height = 1.72)


# read in gab1-idr line scan data
linescan_idr <- read.csv("maria_data/maria_figures/slb_gab1LineScans/slb_gab1_lineScan_GAB1-IDR.csv") %>%
  rename("microns" = "Distance..μm.",
         "EGFR.0" = "EGFR",
         "N.WASP.0" = "N.WASP",
         "Gab1.IDR.0" = "Gab1.IDR") %>%
  select(microns, contains("EGFR"), contains("WASP"), contains("Gab")) %>%
  filter(is.na(microns) == F) %>%
  pivot_longer(cols = contains("."), names_to = "condition", values_to = "intensity") %>%
  mutate(condition = gsub("N.WASP", "N-WASP", condition),
         condition = gsub("Gab1.IDR", "GAB1-IDR", condition),
         microns = microns - (max(microns) - min(microns)) / 2) %>%
  separate(condition, into = c("marker", "replicate"), sep = "\\.", remove = F) %>%
  group_by(microns, marker) %>%
  summarize(mean_intensity = mean_cl_boot(intensity)) %>%
  unnest_wider(mean_intensity, names_sep = "_") %>%
  ungroup() %>%
  mutate(marker = fct_relevel(as.factor(marker), "EGFR", "N-WASP", "GAB1-IDR"))

## supplemental figure 2a: first row, fourth column
# plotting gab1-idr line scan data
plot_linescan_idr <- ggplot(linescan_idr, 
                            aes(x = microns, y = mean_intensity_y, color = marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_intensity_ymin, ymax = mean_intensity_ymax),
                linewidth = 0.1) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none") +
  labs(x = "Distance (um)", y = "Intensity", color = "",
       title = "+ 1 uM GAB1-IDR"); plot_linescan_idr
ggsave("maria_data/maria_figures/lineScan_slb_GAB1-IDR.pdf",
       plot_linescan_idr, units = "in",
       width = 1.49, height = 1.72)


## ----sessionInfo------------------------------------------------------------------------------------------------------------
sessionInfo()

