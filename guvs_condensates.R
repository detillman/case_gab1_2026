## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(ggpubr)
library(tidyverse)


## ----loadData---------------------------------------------------------------------------------------------------------------
# making metadata dataframe
guv_metadata <- data.frame("well" = c("d4", "d5", "d6", "e4", "e5", "e6",
                                           "f4", "f5", "f6", "f7", "f8", "f9",
                                           "g4", "g5", "g6", "g7", "g8", "g9"),
                                "gab1" = c("FL", "FL", "FL", "dPH", "dPH", "dPH",
                                           "FL", "FL", "FL", "FL", "FL", "FL",
                                           "dPH", "dPH", "dPH", "dPH", "dPH", "dPH"),
                                "slb" = c("POPC", "POPC", "POPC",
                                          "POPC", "POPC", "POPC",
                                          "PIP2", "PIP2", "PIP2",
                                          "PIP3", "PIP3", "PIP3",
                                          "PIP2", "PIP2", "PIP2",
                                          "PIP3", "PIP3", "PIP3"))

# loading guv roi frap data
guv_condensate_data <- read.csv("guvs/024_frap/024_guvCircle_data.csv") %>%
  separate(file, into = c("well", "condition", "frap", "channel", "fov", "replicate"), 
           remove = F) %>%
  select(-frap) %>%
  left_join(guv_metadata, by = "well") %>%
  mutate(boxIntensityPixel = boxMean * boxAreaPixel,
         boxIntensityPerSquareMicron = boxIntensityPixel / boxAreaMicron) %>%
  mutate(marker = case_when(channel == "488" ~ "pEGFR",
                            channel == "561" ~ "GRB2",
                            channel == "639" ~ "GAB1"),
         phCondition = case_when(condition == "guv" ~ "guv",
                                  condition == "postGrb2" ~ "postGrb2",
                                  condition == "postGab1" & gab1 == "FL" ~ "postGab1-FL",
                                  condition == "postGab1" & gab1 == "dPH" ~ "postGab1-dPH")) %>% 
  mutate(gab1 = fct_relevel(as.factor(gab1), "FL", "dPH"),
         slb = fct_relevel(as.factor(slb), "POPC", "PIP2", "PIP3"),
         condition = fct_relevel(as.factor(condition), "guv", "postGrb2", "postGab1"),
         marker = fct_relevel(as.factor(marker), "pEGFR", "GRB2", "GAB1"),
         phCondition = fct_relevel(as.factor(phCondition), "guv", "postGrb2", 
                                   "postGab1-FL", "postGab1-dPH")) %>%
  unite("comp", "gab1", "slb", "condition", "marker", remove = F)


## ----analyzeData------------------------------------------------------------------------------------------------------------
# summarizing condensate metrics per guv
guv_condensate_summary <- guv_condensate_data %>%
  group_by(file, well, comp, condition, channel, fov, replicate, gab1, slb, marker, phCondition) %>%
  summarize(mean_IntensityPerSquareMicron = mean(boxIntensityPerSquareMicron),
            sd_IntensityPerSquareMicron = sd(boxIntensityPerSquareMicron)) %>%
  mutate(cv_IntensityPerSquareMicron = sd_IntensityPerSquareMicron / mean_IntensityPerSquareMicron,
         cv2_IntensityPerSquareMicron = (cv_IntensityPerSquareMicron)^2) %>%
  ungroup()

# making plotting dataframe for final figures
guv_condensate_summary_plot <- guv_condensate_summary %>%
  mutate(phCondition = case_when(phCondition == "guv" ~ "+EGFR",
                                 phCondition == "postGrb2" ~ "+GRB2/N-WASP",
                                 phCondition == "postGab1-FL" ~ "+GAB1-WT",
                                 phCondition == "postGab1-dPH" ~ "+GAB1-dPH")) %>%
  mutate(phCondition = fct_relevel(as.factor(phCondition), "+EGFR", "+GRB2/N-WASP",
                                   "+GAB1-WT", "+GAB1-dPH"))

## figure 5b: right
gab1_cv_plot <- ggplot(guv_condensate_summary_plot %>% 
         filter(marker == "GAB1") %>%
           filter(phCondition != "+GAB1-dPH"), 
       aes(x = slb, y = cv2_IntensityPerSquareMicron)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")),
                     label = "p.signif", method = "t.test",
                     size = 2) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GUV Type", y = "GAB1-WT CV2"); gab1_cv_plot
ggsave("manuscript/non_maria_figures/guvs/gab1_cv2.pdf", gab1_cv_plot, 
       units = "in", height = 2.29, width = 2.36)

## figure 5a
egfr_popc_cv_plot <- ggplot(guv_condensate_summary_plot %>% 
                              filter(marker == "pEGFR" & slb == "POPC"), 
                            aes(x = phCondition, y = cv2_IntensityPerSquareMicron)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("+EGFR", "+GRB2/N-WASP"),
                                        c("+EGFR", "+GAB1-WT"),
                                        c("+EGFR", "+GAB1-dPH"),
                                        c("+GRB2/N-WASP", "+GAB1-WT"),
                                        c("+GRB2/N-WASP", "+GAB1-dPH"),
                                        c("+GAB1-WT", "+GAB1-dPH")),
                     label = "p.signif", method = "t.test",
                     size = 2) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Condition", y = "EGFR CV2"); egfr_popc_cv_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_popc_cv2.pdf", egfr_popc_cv_plot,
       units = "in", height = 2.29, width = 2.36)

## figure 5b: left
egfr_pips_cv_plot <- ggplot(guv_condensate_summary_plot %>% 
                              filter(marker == "pEGFR" & phCondition == "+GAB1-WT"), 
                            aes(x = slb, y = cv2_IntensityPerSquareMicron)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")),
                     label = "p.signif", method = "t.test",
                     size = 2) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GUV Type", y = "EGFR CV2 (+GAB1-WT)"); egfr_pips_cv_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_pips_cv_plot.pdf", egfr_pips_cv_plot,
       units = "in", height = 2.29, width = 2.36)

## supplemental figure 5c
egfr_pip2_cv_plot <- ggplot(guv_condensate_summary_plot %>% 
                              filter(marker == "pEGFR" & slb == "PIP2"), 
                            aes(x = phCondition, y = cv2_IntensityPerSquareMicron)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("+EGFR", "+GRB2/N-WASP"),
                                        c("+EGFR", "+GAB1-WT"),
                                        c("+EGFR", "+GAB1-dPH"),
                                        c("+GRB2/N-WASP", "+GAB1-WT"),
                                        c("+GRB2/N-WASP", "+GAB1-dPH"),
                                        c("+GAB1-WT", "+GAB1-dPH")),
                     label = "p.signif", method = "t.test",
                     size = 2) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Condition", y = "EGFR CV2"); egfr_pip2_cv_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_pip2_cv2.pdf", egfr_pip2_cv_plot,
       units = "in", height = 2.29, width = 1.75)

## supplemental figure 5d
egfr_pip3_cv_plot <- ggplot(guv_condensate_summary_plot %>% 
                              filter(marker == "pEGFR" & slb == "PIP3"), 
                            aes(x = phCondition, y = cv2_IntensityPerSquareMicron)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("+EGFR", "+GRB2/N-WASP"),
                                        c("+EGFR", "+GAB1-WT"),
                                        c("+EGFR", "+GAB1-dPH"),
                                        c("+GRB2/N-WASP", "+GAB1-WT"),
                                        c("+GRB2/N-WASP", "+GAB1-dPH"),
                                        c("+GAB1-WT", "+GAB1-dPH")),
                     label = "p.signif", method = "t.test",
                     size = 2) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Condition", y = "EGFR CV2"); egfr_pip3_cv_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_pip3_cv2.pdf", egfr_pip3_cv_plot,
       units = "in", height = 2.29, width = 1.75)


## ----sessionInfo------------------------------------------------------------------------------------------------------------
sessionInfo()

