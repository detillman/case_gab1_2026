## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(ggpubr)
library(minpack.lm)
library(tidyverse)


## ----loadData---------------------------------------------------------------------------------------------------------------
# loading 018 
preFRAP_018 <- read_excel("slbFRAP/018_frapAnalysis.xlsx", sheet = "preFRAP") %>%
  unite("Sample", "Condition", "Well", "Channel", "Replicate", "ROI", remove = F) %>%
  mutate(Condition = as.factor(Condition),
         Marker = case_when(Channel == 488 ~ "pEGFR",
                            Channel == 561 ~ "GRB2",
                            Channel == 639 ~ "GAB1"),
         Channel = as.factor(Channel)) %>%
  mutate(Condition = fct_relevel(Condition, "SLB", "postGrb2", "postGab1"),
         Marker = fct_relevel(as.factor(Marker), "pEGFR", "GRB2", "GAB1"))
frap_018 <- read_excel("slbFRAP/018_frapAnalysis.xlsx", sheet = "FRAP") %>%
  unite("Sample", "Condition", "Well", "Channel", "Replicate", "ROI", remove = F) %>%
  mutate(Condition = as.factor(Condition),
         FRAP_Seconds = as.numeric(FRAP_Seconds),
         Ref_Mean = as.numeric(Ref_Mean),
         Marker = case_when(Channel == 488 ~ "pEGFR",
                            Channel == 561 ~ "GRB2",
                            Channel == 639 ~ "GAB1"),
         Channel = as.factor(Channel)) %>%
  mutate(Condition = fct_relevel(Condition, "SLB", "postGrb2", "postGab1"),
         Marker = fct_relevel(as.factor(Marker), "pEGFR", "GRB2", "GAB1"),
         roi_correct = ROI_Mean / Total_Mean,
         ref_correct = Ref_Mean / Total_Mean)


## ----frapQuantIntial--------------------------------------------------------------------------------------------------------
# loading frap data
frap_018 <- read_excel("slbFRAP/018_frapAnalysis.xlsx", sheet = "FRAP") %>%
  unite("Sample", "Condition", "Well", "Channel", "Replicate", "ROI", remove = F) %>%
  mutate(Condition = as.factor(Condition),
         FRAP_Seconds = as.numeric(FRAP_Seconds),
         Ref_Mean = as.numeric(Ref_Mean),
         Back_Mean = as.numeric(Back_Mean),
         Marker = case_when(Channel == 488 ~ "pEGFR",
                            Channel == 561 ~ "GRB2",
                            Channel == 639 ~ "GAB1"),
         Channel = as.factor(Channel)) %>%
  mutate(Condition = fct_relevel(Condition, "SLB", "postGrb2", "postGab1"),
         Marker = fct_relevel(as.factor(Marker), "pEGFR", "GRB2", "GAB1"),
         Ref_Mean = case_when(is.na(Ref_Mean) == T ~ 1, T ~ Ref_Mean)) %>%
  mutate(roi_over_ref = ROI_Mean / Ref_Mean)

# acquiring pre-bleach roi normalization factors
frap_preRoiFactors_018 <- frap_018 %>%
  filter(Frame == 5) %>%
  select(Sample, roi_over_ref) %>%
  rename(norm_factor_roi_pre = roi_over_ref)

# acquiring post-bleach roi normalization factors
frap_postRoiFactors_018 <- frap_018 %>%
  filter(Frame == 6) %>%
  select(Sample, roi_over_ref) %>%
  rename(norm_factor_roi_post = roi_over_ref)

# applying normalization factors
frap_018_norm <- frap_018 %>%
  left_join(frap_preRoiFactors_018, by = "Sample") %>%
  left_join(frap_postRoiFactors_018, by = "Sample") %>%
  mutate(roi_norm = (roi_over_ref - norm_factor_roi_post) / (norm_factor_roi_pre - norm_factor_roi_post)) %>%
  unite(image_name, Condition, Well, Channel, Replicate, remove = F)

# calculating mean post-frap trace  **FOR VISUALIZATION ONLY**
frap_018_mean_post <- frap_018_norm %>%
  filter(is.na(FRAP_Seconds) == F) %>%
  group_by(Condition, Marker, Frame) %>%
  summarize(mean_time = mean(FRAP_Seconds),
            mean_norm = mean_cl_boot(roi_norm)) %>%
  unnest_wider(mean_norm, names_sep = "_")

# visualzing mean frap trace **FOR VISUALIZATION ONLY**
ggplot(frap_018_mean_post, 
       aes(x = mean_time, y = mean_norm_y,
                                color = Marker)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_norm_ymin, ymax = mean_norm_ymax)) +
  geom_line() +
  facet_wrap(~ Condition) +
  theme_classic() +
  scale_color_manual(values = c("#008000", "#800000", "#800080")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

## mean tau guesses
# slb pEGFR: 2.00 seconds -> k = 0.3465736
# postGrb2 pEGFR: 42.0 seconds -> k = 0.0165035
# postGrb2 GRB2: 8.64 seconds -> k = 0.08022537
# postGab1 GRB2: 4.32 seconds -> k = 0.1604507
# postGab1 GAB1: 20.0 seconds -> k = 0.03465736
# mean k = 0.1276821

## mean plateau guesses
# slb pEGFR: 15.975758 seconds -> tau = 1.996970 -> k = 0.3470994
# postGrb2 pEGFR: 110 seconds -> tau = 42 seconds -> k = 0.0165035
# postGrb2 GRB2: 18.36 seconds -> tau = 3.24 seconds -> k = 0.2139343
# postGab1 GRB2: 14.04 seconds -> tau = 2.16 seconds -> k = 0.3209015
# postGab1 GAB1: 18.0 seconds -> tau = 4.0 seconds -> k = 0.1732868
# mean k = 0.2143451


## ----idOutliers-------------------------------------------------------------------------------------------------------------
## plotting individual frap traces to identify outliers
# slb 488 frap for all rois
ggplot(frap_018_norm %>% filter(Condition == "SLB"), 
       aes(x = Total_Seconds, y = roi_norm, color = Sample)) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
# slb 488 frap for all images
ggplot(frap_018_norm %>% filter(Condition == "SLB"), 
       aes(x = Total_Seconds, y = roi_norm, color = image_name)) +
  geom_smooth(span = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
### no outliers for slb 488 frap

# postGrb2 488 frap for all rois
ggplot(frap_018_norm %>% filter(Condition == "postGrb2") %>% filter(Channel == "488"), 
       aes(x = Total_Seconds, y = roi_norm, color = Sample)) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
# postGrb2 488 frap for all images
ggplot(frap_018_norm %>% filter(Condition == "postGrb2") %>% filter(Channel == "488"), 
       aes(x = Total_Seconds, y = roi_norm, color = image_name)) +
  geom_smooth(span = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
### no outliers for postGrb2 488 frap

# postGrb2 561 frap for all rois
ggplot(frap_018_norm %>% filter(Condition == "postGrb2") %>% filter(Channel == "561"), 
       aes(x = Total_Seconds, y = roi_norm, color = Sample)) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
# postGrb2 561 frap for all images
ggplot(frap_018_norm %>% filter(Condition == "postGrb2") %>% filter(Channel == "561"), 
       aes(x = Total_Seconds, y = roi_norm, color = image_name)) +
  geom_smooth(span = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
### no outliers for postGrb2 561 frap

# postGab1 561 frap for all rois
ggplot(frap_018_norm %>% filter(Condition == "postGab1") %>% filter(Channel == "561"), 
       aes(x = Total_Seconds, y = roi_norm, color = Sample)) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

# postGab1 561 frap for all images
ggplot(frap_018_norm %>% filter(Condition == "postGab1") %>% filter(Channel == "561"), 
       aes(x = Total_Seconds, y = roi_norm, color = image_name)) +
  geom_smooth(span = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
### no outliers for postGab1 561 frap

# postGab1 639 frap for all rois
ggplot(frap_018_norm %>% filter(Condition == "postGab1") %>% filter(Channel == "639"), 
       aes(x = Total_Seconds, y = roi_norm, color = Sample)) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

# postGab1 639 frap for all images
ggplot(frap_018_norm %>% filter(Condition == "postGab1") %>% filter(Channel == "639"), 
       aes(x = Total_Seconds, y = roi_norm, color = image_name)) +
  geom_smooth(span = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")
### no outliers for postGab1 561 frap


## ----frapQuantIndivTrim-----------------------------------------------------------------------------------------------------
## mean plateau guesses
# slb pEGFR: 15.975758 seconds -> tau = 1.996970 -> k = 0.3470994
# postGrb2 pEGFR: 110 seconds -> tau = 42 seconds -> k = 0.0165035
# postGrb2 GRB2: 18.36 seconds -> tau = 3.24 seconds -> k = 0.2139343
# postGab1 GRB2: 14.04 seconds -> tau = 2.16 seconds -> k = 0.3209015
# postGab1 GAB1: 18.0 seconds -> tau = 4.0 seconds -> k = 0.1732868
# mean k = 0.2143451

# trimming frap curves based on time when plateau is reached
frap_018_norm_trim <- frap_018_norm %>%
  mutate(trim_keep = case_when(Condition == "SLB" & Marker == "pEGFR" & FRAP_Seconds <= 30 ~ T,
                               Condition == "SLB" & Marker == "pEGFR" & FRAP_Seconds > 30 ~ F,
                               
                               Condition == "postGrb2" & Marker == "pEGFR" & FRAP_Seconds <= 120 ~ T,
                               Condition == "postGrb2" & Marker == "pEGFR" & FRAP_Seconds > 120 ~ F,
                               
                               Condition == "postGrb2" & Marker == "GRB2" & FRAP_Seconds <= 30 ~ T,
                               Condition == "postGrb2" & Marker == "GRB2" & FRAP_Seconds > 30 ~ F,
                               
                               Condition == "postGab1" & Marker == "GRB2" & FRAP_Seconds <= 30 ~ T,
                               Condition == "postGab1" & Marker == "GRB2" & FRAP_Seconds > 30 ~ F,
                               
                               Condition == "postGab1" & Marker == "GAB1" & FRAP_Seconds <= 30 ~ T,
                               Condition == "postGab1" & Marker == "GAB1" & FRAP_Seconds > 30 ~ F,
                               
                               is.na(FRAP_Seconds) == T ~ T)) %>%
  filter(trim_keep == T)

# function to fit frap curves with nls
fit_frap_nls <- function(frap_df) {

  # guessing initial parameters
  bleach_guess <- mean(head(frap_df$roi_norm, 3))
  end_guess <- mean(tail(frap_df$roi_norm, 3))
  if (end_guess > 1) {end_guess <- 1}
  k_guess <- 0.2

  # fitting curve
  frap_fit <- try(nlsLM(roi_norm ~ y_bleach + (y_end - y_bleach) * (1 - exp(-k * FRAP_Seconds)),
                      data = frap_df,
                      start = list(y_bleach = bleach_guess,
                                   y_end = end_guess,
                                   k = k_guess),
                      control = nls.lm.control(maxiter = 100)),
                  silent = TRUE)

  # return null values if fit fails
  if(inherits(frap_fit, "try-error")) {
    return(NULL)
  }
  
  # extract fit coefficients
  frap_coef <- coef(frap_fit)
  frap_rss <- sum(resid(frap_fit)^2)
  frap_rmse <- sqrt(mean(resid(frap_fit)^2))
  
  # save final results
  frap_results <- tibble(y_bleach = unname(frap_coef["y_bleach"]),
                         y_end = unname(frap_coef["y_end"]),
                         k = unname(frap_coef["k"])) %>%
    mutate(t_half = log(2) / k,
           rss = frap_rss,
           rmse = frap_rmse) 
  frap_results
}

# running frap fitting
frap_018_quant_trim <- frap_018_norm_trim %>%
  filter(is.na(FRAP_Seconds) == F) %>%
  group_by(Sample) %>%
  group_modify(~ {
    frap_results <- fit_frap_nls(.x)
    if (is.null(frap_results)) {
      tibble(y_bleach = NaN, y_end = NaN, k = NaN, t_half = NaN)
    }
    else {frap_results}
  }) %>%
  ungroup() %>%
  mutate(mobile_fraction = (y_end - y_bleach) / (1-y_bleach),
         immobile_fraction = 1 - mobile_fraction) %>%
  separate(Sample, into = c("Condition", "Well", "Channel", "Replicate", "ROI"), 
           remove = F)


## ----validFRAP--------------------------------------------------------------------------------------------------------------
# adding validation metrics and removing outliers
frap_018_validate_trim <- frap_018_norm_trim %>%
  left_join(frap_018_quant_trim %>%
              select(-Condition, -Well, -Channel, -Replicate, -ROI), 
            by = "Sample") %>%
  mutate(roi_fit = y_bleach + (y_end - y_bleach) * (1 - exp(-k * FRAP_Seconds))) %>%
  unite("plot_name", "Well", "Replicate", "ROI", remove = F) %>%
  mutate(plot_name = gsub("10", "", plot_name)) %>%
  filter(Sample != "postGrb2_E10_488_2_3",
         Sample != "postGrb2_D10_488_2_5",
         Sample != "postGrb2_E10_488_1_2",
         Sample != "postGrb2_E10_488_2_2",
         
         Sample != "postGrb2_D10_561_1_2",
         
         Sample != "postGab1_D10_561_1_3",
         Sample != "postGab1_E10_561_2_4",
         
         Sample != "postGab1_D10_639_2_1",
         Sample != "postGab1_D10_639_1_3",
         Sample != "postGab1_D10_639_1_5",
         Sample != "postGab1_D10_639_2_3",
         Sample != "postGab1_D10_639_2_4",
         Sample != "postGab1_C10_639_1_1",
         Sample != "postGab1_C10_639_1_4",
         Sample != "postGab1_D10_639_2_2",
         Sample != "postGab1_E10_639_2_2",
         Sample != "postGab1_E10_639_2_4",
         Sample != "postGab1_D10_639_1_3",
         Sample != "postGab1_D10_639_1_2",
         Sample != "postGab1_E10_639_2_1")

# validating slb pEGFR frap for all rois
ggplot(frap_018_validate_trim %>% 
         filter(is.na(FRAP_Seconds) == F) %>% 
         filter(Condition == "SLB"), 
       aes(x = FRAP_Seconds, y = roi_norm, color = plot_name)) +
  geom_line() +
  geom_line(data = frap_018_validate_trim %>% 
              filter(is.na(FRAP_Seconds) == F) %>% 
              filter(Condition == "SLB"), aes(x = FRAP_Seconds, y = roi_fit),
            color = "black") +
  facet_wrap(~ plot_name, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

# validating postGrb2 pEGFR frap for all rois
ggplot(frap_018_validate_trim %>% 
         filter(is.na(FRAP_Seconds) == F) %>% 
         filter(Condition == "postGrb2") %>%
         filter(Marker == "pEGFR"), 
       aes(x = FRAP_Seconds, y = roi_norm, color = plot_name)) +
  geom_line() +
  geom_line(data = frap_018_validate_trim %>% 
              filter(is.na(FRAP_Seconds) == F) %>% 
              filter(Condition == "postGrb2") %>%
              filter(Marker == "pEGFR"), 
            aes(x = FRAP_Seconds, y = roi_fit),
            color = "black") +
  facet_wrap(~ plot_name, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

# validating postGrb2 GRB2 frap for all rois
ggplot(frap_018_validate_trim %>% 
         filter(is.na(FRAP_Seconds) == F) %>% 
         filter(Condition == "postGrb2") %>%
         filter(Marker == "GRB2"), 
       aes(x = FRAP_Seconds, y = roi_norm, color = plot_name)) +
  geom_line() +
  geom_line(data = frap_018_validate_trim %>% 
              filter(is.na(FRAP_Seconds) == F) %>% 
              filter(Condition == "postGrb2") %>%
              filter(Marker == "GRB2"), 
            aes(x = FRAP_Seconds, y = roi_fit),
            color = "black") +
  facet_wrap(~ plot_name, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

# validating postGab1 GRB2 frap for all rois
ggplot(frap_018_validate_trim %>% 
         filter(is.na(FRAP_Seconds) == F) %>% 
         filter(Condition == "postGab1") %>%
         filter(Marker == "GRB2"), 
       aes(x = FRAP_Seconds, y = roi_norm, color = plot_name)) +
  geom_line() +
  geom_line(data = frap_018_validate_trim %>% 
              filter(is.na(FRAP_Seconds) == F) %>% 
              filter(Condition == "postGab1") %>%
              filter(Marker == "GRB2"), 
            aes(x = FRAP_Seconds, y = roi_fit),
            color = "black") +
  facet_wrap(~ plot_name, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")

# validating postGab1 GAB1 frap for all rois
ggplot(frap_018_validate_trim %>% 
         filter(is.na(FRAP_Seconds) == F) %>% 
         filter(Condition == "postGab1") %>%
         filter(Marker == "GAB1"), 
       aes(x = FRAP_Seconds, y = roi_norm, color = plot_name)) +
  geom_line() +
  geom_line(data = frap_018_validate_trim %>% 
              filter(is.na(FRAP_Seconds) == F) %>% 
              filter(Condition == "postGab1") %>%
              filter(Marker == "GAB1"), 
            aes(x = FRAP_Seconds, y = roi_fit),
            color = "black") +
  facet_wrap(~ plot_name, scales = "free_y") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity", color = "")


## ----frapAnalysis-----------------------------------------------------------------------------------------------------------
# calculating mean trimmed validation frap trace  **FOR VISUALIZATION ONLY**
frap_018_mean_validate_trim <- frap_018_validate_trim %>%
  group_by(Condition, Marker, Frame) %>%
  summarize(mean_time = mean(Total_Seconds),
            mean_frapTime = mean(FRAP_Seconds),
            mean_norm = mean_cl_boot(roi_norm)) %>%
  unnest_wider(mean_norm, names_sep = "_")

# organizing quantification for plotting
frap_018_quantPlot <- frap_018_quant_trim %>%
  filter(Sample %in% frap_018_validate_trim$Sample) %>%
  separate(Sample, into = c("condition", "well", "channel", "replicate", "roi"), 
           remove = F) %>%
  mutate(marker = case_when(channel == 488 ~ "pEGFR",
                            channel == 561 ~ "GRB2",
                            channel == 639 ~ "GAB1")) %>%
  mutate(marker = fct_relevel(as.factor(marker), "pEGFR", "GRB2", "GAB1"),
         condition = fct_relevel(as.factor(condition), "SLB", "postGrb2", "postGab1")) %>%
  filter(y_bleach < 1 & is.na(y_bleach) == F) 

## finding representative imgages
# getting per-condition means
frap_mean_summary <- frap_018_quantPlot %>%
  group_by(condition, marker) %>%
  summarize(half_mean = mean(t_half),
            mobile_mean = mean(mobile_fraction))

# joining per-condition means
frap_repImages <- frap_018_quantPlot %>%
  left_join(frap_mean_summary, by = c("condition", "marker"))

## figure 4b
# plotting grb2 and gab1 frap traces
slbFrap_plot <- ggplot(frap_018_mean_validate_trim %>%
                               ungroup() %>%
                               mutate(Marker = case_when(Marker == "pEGFR" ~ "EGFR", 
                                                         (Marker == "GRB2" & Condition == "postGrb2") ~ 
                                                           "GRB2 (Pre-GAB1)",
                                                         (Marker == "GRB2" & Condition == "postGab1") ~
                                                           "GRB2 (Post-GAB1)",
                                                         T ~ Marker)) %>%
                               mutate(Marker = fct_relevel(as.factor(Marker), "EGFR", "GRB2 (Pre-GAB1)",
                                                           "GRB2 (Post-GAB1)", "GAB1")) %>%
                               mutate(Condition = case_when(Condition == "postGrb2" ~ "postGRB2",
                                                            Condition == "postGab1" ~ "postGAB1",
                                                            T ~ Condition)) %>%
                               filter(is.na(mean_frapTime) == F) %>%
                               filter(Marker != "EGFR") %>%
                               filter(mean_frapTime <= 30) %>%
                               filter(Condition == "postGAB1" | Condition == "postGRB2"), 
                             aes(x = mean_frapTime, y = mean_norm_y,
                                 color = Marker)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_norm_ymin, ymax = mean_norm_ymax),
                linewidth = 0.1, width = 0.5) +
  geom_line(linewidth = 0.25) +
  theme_classic(base_size = 7) +
  scale_color_manual(values = c("#BF8A8B", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = c(0.1, 0.9)) +
  labs(x = "Time (seconds)", y = "Normalized FRAP Intensity",
       color = ""); slbFrap_plot
ggsave("manuscript/non_maria_figures/slbFrap_plot.pdf", slbFrap_plot,
       units = "in", height = 1.73, width = 2.97)

## figure 4e
# plotting grb2 and gab1 mobile fractions
  slbMF_plot <- ggplot(frap_018_quantPlot %>%
                       ungroup() %>%
                               mutate(marker = case_when(marker == "pEGFR" ~ "EGFR", 
                                                         (marker == "GRB2" & Condition == "postGrb2") ~ 
                                                           "GRB2 (Pre-GAB1)",
                                                         (marker == "GRB2" & Condition == "postGab1") ~
                                                           "GRB2 (Post-GAB1)",
                                                         T ~ marker)) %>%
                               mutate(marker = fct_relevel(as.factor(marker), "EGFR", "GRB2 (Pre-GAB1)",
                                                           "GRB2 (Post-GAB1)", "GAB1")) %>%
                               mutate(Condition = case_when(Condition == "postGrb2" ~ "postGRB2",
                                                            Condition == "postGab1" ~ "postGAB1",
                                                            T ~ Condition)) %>%
                       filter(marker != "EGFR"),
       aes(x = marker, y = mobile_fraction * 100, fill = marker)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("GRB2 (Pre-GAB1)", "GRB2 (Post-GAB1)"),
                                        c("GRB2 (Pre-GAB1)", "GAB1"),
                                        c("GRB2 (Post-GAB1)", "GAB1")),
                     label = "p.signif", method = "wilcox", # wilcox has no change
                     size = 2) +
  scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  labs(x = "Condition", y = "Mobile Fraction (%)", fill = "") +
  scale_fill_manual(values = c("#BF8A8B", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); slbMF_plot
ggsave("manuscript/non_maria_figures/slbMF_plot.pdf", slbMF_plot, 
       units = "in", height = 1.73, width = 2.97)

## supplemental figure 4a
# plotting grb2 and gab1 half-times
slbHT_plot <- ggplot(frap_018_quantPlot %>%
                       ungroup() %>%
                               mutate(marker = case_when(marker == "pEGFR" ~ "EGFR", 
                                                         (marker == "GRB2" & Condition == "postGrb2") ~ 
                                                           "GRB2 (Pre-GAB1)",
                                                         (marker == "GRB2" & Condition == "postGab1") ~
                                                           "GRB2 (Post-GAB1)",
                                                         T ~ marker)) %>%
                               mutate(marker = fct_relevel(as.factor(marker), "EGFR", "GRB2 (Pre-GAB1)",
                                                           "GRB2 (Post-GAB1)", "GAB1")) %>%
                               mutate(Condition = case_when(Condition == "postGrb2" ~ "postGRB2",
                                                            Condition == "postGab1" ~ "postGAB1",
                                                            T ~ Condition)) %>%
                       filter(marker != "EGFR"),
       aes(x = marker, y = t_half, fill = marker)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("GRB2 (Pre-GAB1)", "GRB2 (Post-GAB1)"),
                                        c("GRB2 (Pre-GAB1)", "GAB1"),
                                        c("GRB2 (Post-GAB1)", "GAB1")),
                     label = "p.signif", method = "wilcox", # wilcox has changes Pre/GAB1 from * to **
                     size = 2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  labs(x = "Condition", y = "t1/2 (seconds)", fill = "") +
  scale_fill_manual(values = c("#BF8A8B", "#7f1416", "#7c287d")) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none"); slbHT_plot
ggsave("manuscript/non_maria_figures/slbHT_plot.pdf", slbHT_plot, 
       units = "in", height = 2.33, width = 1.33)


## ----sessionInfo------------------------------------------------------------------------------------------------------------
sessionInfo()

