## ----setup, include=FALSE---------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readxl)
library(ggpubr)
library(minpack.lm)
library(tidyverse)


## ----loadData---------------------------------------------------------------------------------------------------------------
# making metadata dataframe
guv_frap_metadata <- data.frame("well" = c("d4", "d5", "d6", "e4", "e5", "e6",
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

# loading guv fov frap data
guv_frap_fov <- read.csv("guvs/024_frap/024_guvFOV_data.csv") %>%
  separate(file, into = c("well", "condition", "frap", "channel", "fov")) %>%
  select(-frap)

# loading guv roi frap data
guv_frap_data <- read.csv("guvs/024_frap/024_guvFRAP_data.csv") %>%
  filter(roiMean > 0) %>%
  separate(file, into = c("well", "condition", "frap", "channel", "fov", "replicate"), 
           remove = F) %>%
  select(-frap) %>%
  left_join(guv_frap_fov, 
            by = c("well", "condition", "channel", "fov", "slice", "time")) %>%
  mutate(corMean = roiMean - fovMean,
         frapTime = case_when(channel == "488" ~ 2*slice-12, # post-bleach is slice6; 2 second interval
                              channel == "561" ~ 0.5*slice-3, # post-bleach is slice6; 0.5 second interval
                              channel == "639" ~ 5*slice-30), # post-bleach is slice 6, 5 second interval
         marker = case_when(channel == "488" ~ "pEGFR",
                            channel == "561" ~ "GRB2",
                            channel == "639" ~ "GAB1")) %>% 
  left_join(guv_frap_metadata, by = "well") %>%
  mutate(gab1 = fct_relevel(as.factor(gab1), "FL", "dPH"),
         slb = fct_relevel(as.factor(slb), "POPC", "PIP2", "PIP3"),
         condition = fct_relevel(as.factor(condition), "guv", "postGrb2", "postGab1"),
         marker = fct_relevel(as.factor(marker), "pEGFR", "GRB2", "GAB1")) %>%
  unite("comp", "gab1", "slb", "condition", "marker", remove = F)


## ----frapQuantIntial--------------------------------------------------------------------------------------------------------
# quantifying pre-bleach normalization factors for original and corrected roi values
guv_frap_normPre <- guv_frap_data %>%
  filter(slice == 5) %>%
  select(file, roiMean, corMean) %>%
  rename(normPre_roi = roiMean,
         normPre_cor = corMean)

# quantifying post-bleach normalization factors for original and corrected roi values
guv_frap_normPost <- guv_frap_data %>%
  filter(slice == 6) %>%
  select(file, roiMean, corMean) %>%
  rename(normPost_roi = roiMean,
         normPost_cor = corMean)

# applying pre-bleach normalization factors for original and corrected roi values
guv_frap_norm <- guv_frap_data %>%
  left_join(guv_frap_normPre, by = "file") %>%
  left_join(guv_frap_normPost, by = "file") %>%
  mutate(roi_norm = (roiMean - normPost_roi) / (normPre_roi - normPost_roi),
         cor_norm = (corMean - normPost_cor) / (normPre_cor - normPost_cor))

# calculating mean frap traces  **FOR VISUALIZATION ONLY**
guv_frap_mean <- guv_frap_norm %>%
  group_by(comp, time) %>%
  summarize(mean_time = mean(time),
            mean_frapTime = mean(frapTime),
            mean_roi = mean_cl_boot(roi_norm),
            mean_cor = mean_cl_boot(cor_norm)) %>%
  unnest_wider(mean_roi, names_sep = "_") %>%
  unnest_wider(mean_cor, names_sep = "_") %>%
  separate(comp, into = c("gab1", "slb", "condition", "marker"), remove = F)

# plotting each mean post-frap trace to estimate fitting parameters
guv_comps <- guv_frap_mean %>%
  distinct(comp)

for (i in 1:nrow(guv_comps)) {
  
  # selecting comparison for plotting
  plot_comp <- guv_comps [i, ]
  
  # making dataframe for plotting
  plot_df <- guv_frap_mean %>%
    filter(mean_frapTime >= 0) %>%
    filter(comp == plot_comp$comp)
  
  # plotting corrected frap trace
  print(ggplot(plot_df, aes(x = mean_frapTime, y = mean_cor_y)) +
          geom_point(color = "purple") +
          geom_errorbar(color = "purple",
                        aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax)) +
          geom_line(color = "purple") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black")) +
          labs(x = "Time (seconds)", y = "Corrected ROI FRAP Intensity",
               title = plot_comp))
  
}

### estimated fitting parameters are below - refer to spreadsheet for details
## dividing time at plateau value by 2
# POPC GUVs pEGFR: k = 0.082
# PIP2 GUVs pEGFR: k = .056
# PIP3 GUVs pEGFR: k = 0.081
# POPC postGrb2 pEGFR: k = 0.015
# PIP2 postGrb2 pEGFR: k = 0.020
# PIP3 postGrb2 pEGFR: k = 0.013
# POPC postGrb2 GRB2: k = 0.013
# PIP2 postGrb2 GRB2: k = 0.012
# PIP3 postGrb2 GRB2: k = 0.012
# FL POPC postGab1 pEGFR: k = 0.030
# FL PIP2 postGab1 pEGFR: k = 0.14
# FL PIP3 postGab1 pEGFR: k = 0.030
# FL POPC postGab1 GRB2: k = 0.012
# FL PIP2 postGab1 GRB2: k = 0.012
# FL PIP3 postGab1 GRB2: k = 0.014
# FL POPC postGab1 GAB1: k = 0.020
# FL PIP2 postGab1 GAB1: k = 0.0060
# FL PIP3 postGab1 GAB1: k = 0.0051
# dPH POPC postGab1 pEGFR: k = 0.014
# dPH PIP2 postGab1 pEGFR: k = 0.045
# dPH PIP3 postGab1 pEGFR: k = 0.019
# dPH POPC postGab1 GRB2: k = 0.015
# dPH PIP2 postGab1 GRB2: k = 0.012 
# dPH PIP3 postGab1 GRB2: k = 0.012
# dPH POPC postGab1 GAB1: k = 0.0079
# dPH PIP2 postGab1 GAB1: k = 0.0052
# dPH PIP3 postGab1 GAB1: k = NA

## finding tau based on plateau and baseline values
# POPC GUVs pEGFR: k = 0.30
# PIP2 GUVs pEGFR: k = 0.30
# PIP3 GUVs pEGFR: k = 0.35
# POPC postGrb2 pEGFR: k = 0.033
# PIP2 postGrb2 pEGFR: k = 0.044
# PIP3 postGrb2 pEGFR: k = 0.051
# POPC postGrb2 GRB2: k = 0.053
# PIP2 postGrb2 GRB2: k = 0.027
# PIP3 postGrb2 GRB2: k = 0.030
# FL POPC postGab1 pEGFR: k = 0.13
# FL PIP2 postGab1 pEGFR: k = 0.26
# FL PIP3 postGab1 pEGFR: k = 0.19
# FL POPC postGab1 GRB2: k = 0.058
# FL PIP2 postGab1 GRB2: k = 0.030
# FL PIP3 postGab1 GRB2: k = 0.050
# FL POPC postGab1 GAB1: k = 0.046
# FL PIP2 postGab1 GAB1: k = 0.012
# FL PIP3 postGab1 GAB1: k = 0.026
# dPH POPC postGab1 pEGFR: k = 0.056
# dPH PIP2 postGab1 pEGFR: k = 0.14
# dPH PIP3 postGab1 pEGFR: k = 0.09
# dPH POPC postGab1 GRB2: k = 0.029
# dPH PIP2 postGab1 GRB2: k = 0.024
# dPH PIP3 postGab1 GRB2: k = 0.058
# dPH POPC postGab1 GAB1: k = 0.012
# dPH PIP2 postGab1 GAB1: k = 0.010
# dPH PIP3 postGab1 GAB1: k = NA

### final k guesses
# pEGFR: k = 0.10 -> tau = 6.9 seconds
# GRB2: k = 0.041 -> tau = 17 seconds
# GAB1: k = 0.015 -> tau = 46 seconds


## ----idOutliers-------------------------------------------------------------------------------------------------------------
# plotting individual frap traces to identify outliers
guv_comps <- guv_frap_norm %>%
  distinct(comp)

for (i in 1:nrow(guv_comps)) {
  
  # selecting comparison for plotting
  plot_comp <- guv_comps [i, ]
  
  # making dataframe for plotting
  plot_df <- guv_frap_norm %>%
    filter(comp == plot_comp)
  
  # plotting corrected frap trace
  print(ggplot(plot_df, aes(x = time, y = cor_norm, color = file)) +
          geom_line() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black")) +
          labs(x = "Time (seconds)", y = "Corrected ROI FRAP Intensity",
               title = plot_comp))
  
}

# removing outliers
guv_frap_clean <- guv_frap_norm %>%
  filter(file != "d5_postGab1_frap_488_1_frap3") %>%
  filter(file != "e5_postGab1_frap_561_1_frap4") %>%
  filter(file != "e4_postGab1_frap_561_1_frap2") %>%
  filter(file != "e6_postGrb2_frap_561_1_frap3") %>%
  filter(file != "f6_postGab1_frap_488_1_frap4") %>%
  filter(file != "f6_postGab1_frap_488_2_frap1") %>%
  filter(file != "f6_postGab1_frap_488_2_frap2") %>%
  filter(file != "f6_postGab1_frap_561_3_frap2") %>%
  filter(file != "f4_postGab1_frap_561_3_frap2") %>%
  filter(file != "f6_postGab1_frap_561_2_frap1") %>%
  filter(file != "f4_postGab1_frap_639_2_frap1") %>%
  filter(file != "f4_postGab1_frap_639_2_frap2") %>%
  filter(file != "f4_postGab1_frap_639_2_frap3") %>%
  filter(file != "f5_postGab1_frap_639_1_frap2") %>%
  filter(file != "f5_postGab1_frap_639_2_frap1") %>%
  filter(file != "f5_postGab1_frap_639_2_frap3") %>%
  filter(file != "f5_postGab1_frap_639_2_frap4") %>%
  filter(file != "f5_postGab1_frap_639_2_frap5") %>%
  filter(file != "f6_postGab1_frap_639_1_frap1") %>%
  filter(file != "f6_postGab1_frap_639_1_frap2") %>%
  filter(file != "f6_postGab1_frap_639_1_frap3") %>%
  filter(file != "f4_postGab1_frap_639_1_frap1") %>%
  filter(file != "f5_postGrb2_frap_488_1_frap3") %>%
  filter(file != "g5_postGab1_frap_639_1_frap1") %>%
  filter(file != "g5_postGab1_frap_639_1_frap2") %>%
  filter(file != "g5_postGab1_frap_639_1_frap3") %>%
  filter(file != "g6_postGrb2_frap_488_2_frap5")


## ----frapQuantIndiv---------------------------------------------------------------------------------------------------------
### final k guesses
# pEGFR: k = 0.10 -> tau = 6.9 seconds
# GRB2: k = 0.041 -> tau = 17 seconds
# GAB1: k = 0.015 -> tau = 46 seconds
# mean k = 0.052

# function to fit frap curves with nls
fit_frap_nls <- function(frap_df) {

  # guessing initial parameters
  bleach_guess <- min(frap_df$frap_value)
  end_guess <- max(frap_df$frap_value)
  if (end_guess > 1) {end_guess <- 1}
  k_guess <- 0.05

  # fitting curve
  frap_fit <- try(nlsLM(frap_value ~ y_bleach + (y_end - y_bleach) * (1 - exp(-k * frapTime)),
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

# running frap fitting on corrected rois
guv_frap_quant_cor <- guv_frap_clean %>%
  filter(frapTime >= 0) %>%
  rename(frap_value = cor_norm) %>%
  group_by(file) %>%
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
  separate(file, into = c("well", "condition", "frap", "channel", "fov", "replicate"), 
           remove = F) %>%
  select(-frap)


## ----validFRAP--------------------------------------------------------------------------------------------------------------
# adding validation metrics and removing outliers
guv_frap_validate <- guv_frap_clean %>%
  left_join(guv_frap_quant_roi %>% 
              select(-well, -condition, -channel, -fov, -replicate),
            by = "file") %>%
  left_join(guv_frap_quant_cor %>%
              select(-well, -condition, -channel, -fov, -replicate),
            by = "file", suffix = c("_roi", "_cor")) %>%
  mutate(roi_fit = y_bleach_roi + (y_end_roi - y_bleach_roi) * (1 - exp(-k_roi * frapTime)),
         cor_fit = y_bleach_cor + (y_end_cor - y_bleach_cor) * (1 - exp(-k_cor * frapTime)),
         replicate = gsub("frap", "", replicate)) %>%
  unite("image", "well", "fov", "replicate", remove = F)

# plotting individual frap traces for each condition
guv_comps <- guv_frap_validate %>%
  distinct(comp)

for (i in 1:nrow(guv_comps)) {
  
  # selecting comparison for plotting
  plot_comp <- guv_comps [i, ]
  
  # making dataframe for plotting
  plot_df <- guv_frap_validate %>%
    filter(frapTime >= 0) %>%
    filter(comp == plot_comp)
  
  # plotting corrected frap trace
  print(ggplot(plot_df, aes(x = frapTime, y = cor_norm, color = image)) +
          geom_line() +
          geom_line(data = guv_frap_validate %>% 
                      filter(frapTime >= 0) %>%
                      filter(comp == plot_comp), 
                    aes(x = frapTime, y = cor_fit),
            color = "black") +
          theme_classic() +
          facet_wrap(~ image, scales = "free") +
          theme(axis.text.x = element_text(color = "black"),
                axis.text.y = element_text(color = "black"),
                axis.ticks = element_line(color = "black")) +
          labs(x = "Time (seconds)", y = "Corrected ROI FRAP Intensity",
               title = paste0(plot_comp, ": Corrected")))
  
}

# removing poor frap fits
guv_frap_final <- guv_frap_validate %>%
  filter(file != "d4_guv_frap_488_1_frap2") %>%
  filter(file != "d4_guv_frap_488_2_frap4") %>%
  filter(file != "d5_guv_frap_488_2_frap2") %>%
  filter(file != "d6_guv_frap_488_1_frap1") %>%
  filter(file != "d4_postGab1_frap_488_1_frap3") %>%
  filter(file != "d6_postGab1_frap_488_1_frap5") %>%
  filter(file != "d6_postGab1_frap_488_2_frap3") %>%
  filter(file != "d4_postGab1_frap_639_1_frap5") %>%
  filter(file != "d4_postGab1_frap_639_2_frap3") %>%
  filter(file != "d5_postGab1_frap_639_1_frap3") %>%
  filter(file != "d5_postGab1_frap_639_2_frap1") %>% 
  filter(file != "d5_postGab1_frap_639_2_frap2") %>% 
  filter(file != "d5_postGab1_frap_639_2_frap3") %>%
  filter(file != "d6_postGab1_frap_639_1_frap2") %>%
  filter(file != "d6_postGab1_frap_639_1_frap3") %>%
  filter(file != "d6_postGab1_frap_639_1_frap5") %>%
  filter(file != "d4_postGrb2_frap_488_1_frap1") %>%
  filter(file != "d4_postGrb2_frap_488_1_frap2") %>%
  filter(file != "d4_postGrb2_frap_488_2_frap2") %>%
  filter(file != "d4_postGrb2_frap_488_2_frap3") %>%
  filter(file != "d5_postGrb2_frap_488_1_frap3") %>%
  filter(file != "d6_postGrb2_frap_488_1_frap4") %>%
  filter(file != "d6_postGrb2_frap_488_2_frap2") %>%
  filter(file != "d4_postGrb2_frap_561_1_frap1") %>%
  filter(file != "d6_postGrb2_frap_561_2_frap1") %>%
  filter(file != "d6_postGrb2_frap_561_2_frap3") %>%
  filter(file != "e4_guv_frap_488_1_frap1") %>%
  filter(file != "e4_postGab1_frap_488_1_frap3") %>%
  filter(file != "e4_postGab1_frap_488_2_frap2") %>%
  filter(file != "e4_postGab1_frap_488_2_frap3") %>%
  filter(file != "e5_postGab1_frap_488_1_frap3") %>%
  filter(file != "e5_postGab1_frap_488_2_frap1") %>%
  filter(file != "e6_postGab1_frap_488_2_frap1") %>%
  filter(file != "e4_postGab1_frap_561_1_frap1") %>%
  filter(file != "e5_postGab1_frap_561_1_frap3") %>%
  filter(file != "e5_postGab1_frap_561_2_frap3") %>%
  filter(file != "e4_postGab1_frap_639_1_frap1") %>%
  filter(file != "e4_postGrb2_frap_488_1_frap1") %>%
  filter(file != "e4_postGrb2_frap_488_1_frap3") %>%
  filter(file != "e4_postGrb2_frap_488_1_frap5") %>%
  filter(file != "e5_postGrb2_frap_488_1_frap2") %>%
  filter(file != "e5_postGrb2_frap_488_2_frap2") %>%
  filter(file != "e5_postGrb2_frap_488_2_frap5") %>%
  filter(file != "e6_postGrb2_frap_488_1_frap4") %>%
  filter(file != "e5_postGrb2_frap_561_2_frap3") %>%
  filter(file != "f6_guv_frap_488_1_frap1") %>%
  filter(file != "f4_postGab1_frap_488_1_frap1") %>%
  filter(file != "f5_postGab1_frap_488_1_frap3") %>%
  filter(file != "f5_postGab1_frap_488_1_frap4") %>%
  filter(file != "f5_postGab1_frap_488_2_frap4") %>%
  filter(file != "f6_postGab1_frap_488_1_frap1") %>%
  filter(file != "f6_postGab1_frap_488_1_frap3") %>%
  filter(file != "f6_postGab1_frap_488_2_frap3") %>%
  filter(file != "f5_postGab1_frap_561_2_frap3") %>%
  filter(file != "f4_postGab1_frap_639_1_frap2") %>%
  filter(file != "f5_postGab1_frap_639_1_frap1") %>%
  filter(file != "f6_postGab1_frap_639_2_frap4") %>%
  filter(file != "f6_postGab1_frap_639_2_frap5") %>%
  filter(file != "f5_postGrb2_frap_488_2_frap4") %>%
  filter(file != "f6_postGrb2_frap_488_1_frap1") %>%
  filter(file != "f7_guv_frap_488_2_frap1") %>%
  filter(file != "f7_guv_frap_488_2_frap2") %>%
  filter(file != "f7_guv_frap_488_2_frap3") %>%
  filter(file != "f9_guv_frap_488_2_frap2") %>%
  filter(file != "f9_guv_frap_488_2_frap3") %>%
  filter(file != "f7_postGab1_frap_488_2_frap3") %>%
  filter(file != "f7_postGab1_frap_488_2_frap4") %>%
  filter(file != "f8_postGab1_frap_488_2_frap1") %>%
  filter(file != "f8_postGab1_frap_488_2_frap3") %>%
  filter(file != "f9_postGab1_frap_488_1_frap2") %>%
  filter(file != "f7_postGab1_frap_639_1_frap1") %>%
  filter(file != "f7_postGab1_frap_639_2_frap2") %>%
  filter(file != "f8_postGab1_frap_639_1_frap1") %>%
  filter(file != "f9_postGab1_frap_639_1_frap1") %>%
  filter(file != "f9_postGab1_frap_639_2_frap1") %>%
  filter(file != "f9_postGab1_frap_639_3_frap1") %>%
  filter(file != "f9_postGab1_frap_639_3_frap3") %>%
  filter(file != "f7_postGrb2_frap_488_1_frap3") %>%
  filter(file != "f7_postGrb2_frap_488_1_frap4") %>%
  filter(file != "f7_postGrb2_frap_488_2_frap4") %>%
  filter(file != "f9_postGrb2_frap_488_2_frap1") %>%
  filter(file != "f9_postGrb2_frap_488_2_frap3") %>%
  filter(file != "f7_postGrb2_frap_561_2_frap2") %>%
  filter(file != "f9_postGrb2_frap_561_1_frap3") %>%
  filter(file != "g5_guv_frap_488_2_frap3") %>%
  filter(file != "g4_postGab1_frap_488_1_frap1") %>%
  filter(file != "g4_postGab1_frap_488_1_frap2") %>%
  filter(file != "g4_postGab1_frap_488_1_frap5") %>%
  filter(file != "g5_postGab1_frap_488_1_frap1") %>%
  filter(file != "g5_postGab1_frap_488_1_frap2") %>%
  filter(file != "g5_postGab1_frap_488_1_frap3") %>%
  filter(file != "g5_postGab1_frap_488_1_frap4") %>%
  filter(file != "g5_postGab1_frap_488_1_frap5") %>%
  filter(file != "g5_postGab1_frap_488_2_frap1") %>%
  filter(file != "g5_postGab1_frap_488_2_frap5") %>%
  filter(file != "g6_postGab1_frap_488_1_frap1") %>%
  filter(file != "g6_postGab1_frap_488_1_frap4") %>%
  filter(file != "g6_postGab1_frap_488_2_frap3") %>%
  filter(file != "g4_postGab1_frap_561_2_frap2") %>%
  filter(file != "g5_postGab1_frap_561_1_frap1") %>%
  filter(file != "g6_postGab1_frap_561_1_frap1") %>%
  filter(file != "g6_postGab1_frap_561_1_frap2") %>%
  filter(file != "g6_postGab1_frap_561_1_frap3") %>%
  filter(file != "g6_postGab1_frap_561_2_frap1") %>%
  filter(file != "g4_postGab1_frap_639_2_frap2") %>%
  filter(file != "g4_postGrb2_frap_488_1_frap1") %>%
  filter(file != "g4_postGrb2_frap_488_1_frap2") %>%
  filter(file != "g5_postGrb2_frap_488_1_frap2") %>%
  filter(file != "g4_postGrb2_frap_561_3_frap2") %>%
  filter(file != "g5_postGrb2_frap_561_1_frap1") %>%
  filter(file != "g7_guv_frap_488_2_frap2") %>%
  filter(file != "g8_guv_frap_488_2_frap1") %>%
  filter(file != "g9_guv_frap_488_2_frap1") %>%
  filter(file != "g9_postGab1_frap_488_2_frap1") %>%
  filter(file != "g9_postGab1_frap_488_2_frap2") %>%
  filter(file != "g7_postGab1_frap_561_2_frap1") %>%
  filter(file != "g7_postGab1_frap_561_2_frap2") %>% 
  filter(file != "g7_postGab1_frap_561_2_frap2") %>%
  filter(file != "g9_postGab1_frap_561_1_frap1") %>%
  filter(file != "g9_postGab1_frap_561_1_frap2") %>%
  filter(file != "g7_postGrb2_frap_488_1_frap3") %>%
  filter(file != "g7_postGrb2_frap_488_2_frap2") %>%
  filter(file != "g9_postGrb2_frap_561_1_frap3") %>%
  filter(file != "g9_postGrb2_frap_561_2_frap1")


## ----frapAnalysis-----------------------------------------------------------------------------------------------------------
# calculating mean validated frap traces  **FOR VISUALIZATION ONLY**
guv_frap_final_mean <- guv_frap_final %>%
  mutate(phCondition = case_when(condition == "guv" ~ "guv",
                                  condition == "postGrb2" ~ "postGrb2",
                                  condition == "postGab1" & gab1 == "FL" ~ "postGab1-FL",
                                  condition == "postGab1" & gab1 == "dPH" ~ "postGab1-dPH")) %>%
  group_by(slb, marker, phCondition, time) %>%
  summarize(mean_time = mean(time),
            mean_frapTime = mean(frapTime),
            mean_cor = mean_cl_boot(cor_norm)) %>%
  unnest_wider(mean_cor, names_sep = "_") %>%
  ungroup() %>%
  unite("comp", "marker", "slb", "phCondition", remove = F) %>%
  unite("noSLB_comp", "marker", "phCondition", remove = F) %>%
  mutate(phCondition = fct_relevel(phCondition, "guv", "postGrb2", 
                                   "postGab1-FL", "postGab1-dPH"))

# intial visualization of corrected mean frap traces **FOR VISUALIZATION ONLY**
ggplot(guv_frap_final_mean, aes(x = mean_time, y = mean_cor_y,
                          color = phCondition)) +
  #geom_point() +
  #geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax)) +
  geom_line() +
  facet_wrap(~ marker + slb, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Corrected ROI FRAP Intensity")


# final visualization of corrected mean frap traces by comp **FOR VISUALIZATION ONLY**
ggplot(guv_frap_final_mean %>%
         filter(marker == "pEGFR") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = phCondition)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ slb, scales = "free") +
  scale_color_manual(values = c("#008000", "#800000", "#800080", "#C080C0")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized pEGFR FRAP Intensity",
       color = "Condition")

ggplot(guv_frap_final_mean %>%
         filter(marker == "GRB2") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = phCondition)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ slb, scales = "free") +
  scale_color_manual(values = c("#800000", "#800080", "#C080C0")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized GRB2 FRAP Intensity",
       color = "Condition")



# final visualization of corrected mean frap traces by slb **FOR VISUALIZATION ONLY**
ggplot(guv_frap_final_mean %>%
         filter(marker == "pEGFR") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = slb)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ phCondition, scales = "free") +
  #scale_color_manual(values = c("#008000", "#800000", "#800080", "#C080C0")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized pEGFR FRAP Intensity",
       color = "GUV")

ggplot(guv_frap_final_mean %>%
         filter(marker == "GRB2") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = slb)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), alpha = 0.25) +
  geom_line() +
  facet_wrap(~ phCondition, scales = "free") +
  #scale_color_manual(values = c("#800000", "#800080", "#C080C0")) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized GRB2 FRAP Intensity",
       color = "GUV")

ggplot(guv_frap_final_mean %>%
         filter(marker == "GAB1" & phCondition == "postGab1-FL") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = slb)) +
  geom_point(size = 0.5) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), alpha = 0.25) +
  geom_line() +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Time (seconds)", y = "Normalized GAB1 FRAP Intensity",
       color = "GUV")

# summarizing frap quantification metrics for plots
guv_frap_quant_final <- guv_frap_quant_cor %>%
  filter(file %in% guv_frap_final$file) %>%
  left_join(guv_frap_metadata, by = "well") %>%
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
         phCondition = fct_relevel(as.factor(phCondition), "guv", "postGrb2", 
                                   "postGab1-FL", "postGab1-dPH"),
         marker = fct_relevel(as.factor(marker), "pEGFR", "GRB2", "GAB1")) %>%
  unite("comp", "gab1", "slb", "condition", "marker", remove = F)
  
# visualizing how egfr half-time changes across conditions and slbs
ggplot(guv_frap_quant_final %>% filter(marker == "pEGFR"),
       aes(x = phCondition, y = t_half)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("guv", "postGrb2"),
                                        c("guv", "postGab1-FL"),
                                        c("guv", "postGab1-dPH"),
                                        c("postGrb2", "postGab1-FL"),
                                        c("postGrb2", "postGab1-dPH"),
                                        c("postGab1-FL", "postGab1-dPH")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ slb, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "pEGFR Half-Time (seconds)")

ggplot(guv_frap_quant_final %>% filter(marker == "pEGFR"),
       aes(x = slb, y = t_half)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ phCondition, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "pEGFR Half-Time (seconds)")

# visualizing how grb2 half-time changes across conditions and slbs
ggplot(guv_frap_quant_final %>% filter(marker == "GRB2"),
       aes(x = phCondition, y = t_half)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("postGrb2", "postGab1-FL"),
                                        c("postGrb2", "postGab1-dPH"),
                                        c("postGab1-FL", "postGab1-dPH")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ slb, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "GRB2 Half-Time (seconds)")

ggplot(guv_frap_quant_final %>% filter(marker == "GRB2"),
       aes(x = slb, y = t_half)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ phCondition, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "GRB2 Half-Time (seconds)")

# visualizing how gab1 half-time across slbs
ggplot(guv_frap_quant_final %>% filter(marker == "GAB1") %>% filter(gab1 == "FL"),
       aes(x = slb, y = t_half)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")), 
                     label = "p.format", method = "t.test") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "GAB1 Half-Time (seconds)")


# visualizing how egfr mobile fraction changes across conditions and slbs
ggplot(guv_frap_quant_final %>% filter(marker == "pEGFR"),
       aes(x = phCondition, y = mobile_fraction)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("guv", "postGrb2"),
                                        c("guv", "postGab1-FL"),
                                        c("guv", "postGab1-dPH"),
                                        c("postGrb2", "postGab1-FL"),
                                        c("postGrb2", "postGab1-dPH"),
                                        c("postGab1-FL", "postGab1-dPH")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ slb, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "pEGFR Mobile Fraction")


ggplot(guv_frap_quant_final %>% filter(marker == "pEGFR"),
       aes(x = slb, y = mobile_fraction)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ phCondition, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "pEGFR Mobile Fraction")

# visualizing how grb2 half-time changes across conditions and slbs
ggplot(guv_frap_quant_final %>% filter(marker == "GRB2"),
       aes(x = phCondition, y = mobile_fraction)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("postGrb2", "postGab1-FL"),
                                        c("postGrb2", "postGab1-dPH"),
                                        c("postGab1-FL", "postGab1-dPH")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ slb, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "GRB2 Mobile Fraction")

ggplot(guv_frap_quant_final %>% filter(marker == "GRB2"),
       aes(x = slb, y = mobile_fraction)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")), 
                     label = "p.format", method = "t.test") +
  facet_wrap(~ phCondition, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "GRB2 Mobile Fraction")

# visualizing how gab1 half-time across slbs
ggplot(guv_frap_quant_final %>% filter(marker == "GAB1") %>% filter(gab1 == "FL"),
       aes(x = slb, y = mobile_fraction)) +
  geom_bar(stat = "summary", position = "dodge") +
  geom_errorbar(stat = 'summary', width = 0.25, 
                fun.data = mean_cl_boot,
                position = position_dodge(width = 0.9)) +
  geom_point(position = position_dodge(width = 0.9), size = 0.5) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")), 
                     label = "p.format", method = "t.test") +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "", y = "GAB1 Mobile Fraction")

# GAB1-FL makes pEGFR less mobile across all GUV types
# GAB1-dPH makes pEGFR less mobile on PIP2 GUVs
# GAB1-dPH maeks GRB2 less mobile across all GUV types
# GAB1-FL makes GRB2 less mobile on PIP2 GUVs

# pEGFR is less mobile on PIP2 GUVs after adding GAB1-FL or GAB1-dPH
# GRB2 is less mobile on PIP2 GUVs after adding GAB1-FL or GAB1-dPH


# GAB1-FL makes pEGFR recover more slowly on POPC GUVs
# GAB1-dPH makes pEGFR recover more slowly on POPC GUVs
# GAB1-FL makes GRB2 recover more slowly on PIP3 GUVs
# GAB1-dPH makes GRB2 recover more slowly on PIP3 GUVs
# GRB2 recovers more slowly on PIP2 GUVs


## finding representative images
guv_frap_quant_summary <- guv_frap_quant_final %>%
  group_by(phCondition, marker, slb) %>%
  summarize(mean_mobile = mean(mobile_fraction),
            mean_half = mean(t_half))

## making dataframe for plotting final frap traces
guv_frap_final_mean_plot <- guv_frap_final_mean %>%
  mutate(phCondition = case_when(phCondition == "guv" ~ "+EGFR",
                                 phCondition == "postGrb2" ~ "+GRB2/N-WASP",
                                 phCondition == "postGab1-FL" ~ "+GAB1-WT",
                                 phCondition == "postGab1-dPH" ~ "+GAB1-dPH")) %>%
  mutate(phCondition = fct_relevel(as.factor(phCondition), "+EGFR", "+GRB2/N-WASP",
                                   "+GAB1-WT", "+GAB1-dPH"))

## final figures
# egfr corrected mean frap traces across conditions **FOR VISUALIZATION ONLY**
frap_egfr_popc_plot <- ggplot(guv_frap_final_mean_plot %>%
         filter(marker == "pEGFR" & slb == "POPC") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = phCondition)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), 
                #alpha = 0.25,
                linewidth = 0.1, width = 0.5) +
  geom_line(linewidth = 0.25) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d", "#BE94BE")) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "Time (seconds)", y = "Normalized EGFR FRAP Intensity",
       color = ""); frap_egfr_popc_plot
ggsave("manuscript/non_maria_figures/guvs/frap_egfr_popc.pdf", frap_egfr_popc_plot,
       units = "in", height = 1.64, width = 2.36)

frap_egfr_pip2_plot <- ggplot(guv_frap_final_mean_plot %>%
         filter(marker == "pEGFR" & slb == "PIP2") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = phCondition)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), 
                #alpha = 0.25,
                linewidth = 0.1, width = 0.5) +
  geom_line(linewidth = 0.25) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d", "#BE94BE")) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "Time (seconds)", y = "Normalized EGFR FRAP Intensity",
       color = ""); frap_egfr_pip2_plot
ggsave("manuscript/non_maria_figures/guvs/frap_egfr_pip2.pdf", frap_egfr_pip2_plot,
       units = "in", height = 1.64, width = 2.36)

frap_egfr_pip3_plot <- ggplot(guv_frap_final_mean_plot %>%
         filter(marker == "pEGFR" & slb == "PIP3") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = phCondition)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), 
                #alpha = 0.25,
                linewidth = 0.1, width = 0.5) +
  geom_line(linewidth = 0.25) +
  scale_color_manual(values = c("#0c8140", "#7f1416", "#7c287d", "#BE94BE")) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "none") +
  labs(x = "Time (seconds)", y = "Normalized EGFR FRAP Intensity",
       color = ""); frap_egfr_pip3_plot
ggsave("manuscript/non_maria_figures/guvs/frap_egfr_pip3.pdf", frap_egfr_pip3_plot,
       units = "in", height = 1.64, width = 2.36)


## making final dataframe for frap quantification plots
guv_frap_quant_final_plot <- guv_frap_quant_final %>%
  mutate(phCondition = case_when(phCondition == "guv" ~ "+EGFR",
                                 phCondition == "postGrb2" ~ "+GRB2/N-WASP",
                                 phCondition == "postGab1-FL" ~ "+GAB1-WT",
                                 phCondition == "postGab1-dPH" ~ "+GAB1-dPH")) %>%
  mutate(phCondition = fct_relevel(as.factor(phCondition), "+EGFR", "+GRB2/N-WASP",
                                   "+GAB1-WT", "+GAB1-dPH"))

## visualizing how egfr mobile fraction changes across conditions 
egfr_popc_mf_plot <- ggplot(guv_frap_quant_final_plot %>% 
                              filter(marker == "pEGFR" & slb == "POPC"),
       aes(x = phCondition, y = mobile_fraction*100)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("+EGFR", "+GRB2/N-WASP"),
                                        c("+EGFR", "+GAB1-WT"),
                                        c("+EGFR", "+GAB1-dPH"),
                                        c("+GRB2/N-WASP", "+GAB1-WT"),
                                        c("+GRB2/N-WASP", "+GAB1-dPH"),
                                        c("+GAB1-WT", "+GAB1-dPH")),
                     label = "p.signif", method = "t.test", # wilcox changes?
                     size = 2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Condition", y = "EGFR Mobile Fraction (%)"); egfr_popc_mf_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_popc_mf.pdf", egfr_popc_mf_plot,
       units = "in", height = 2.05, width = 2.36)

egfr_pip2_mf_plot <- ggplot(guv_frap_quant_final_plot %>% 
                              filter(marker == "pEGFR" & slb == "PIP2"),
       aes(x = phCondition, y = mobile_fraction*100)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("+EGFR", "+GRB2/N-WASP"),
                                        c("+EGFR", "+GAB1-WT"),
                                        c("+EGFR", "+GAB1-dPH"),
                                        c("+GRB2/N-WASP", "+GAB1-WT"),
                                        c("+GRB2/N-WASP", "+GAB1-dPH"),
                                        c("+GAB1-WT", "+GAB1-dPH")),
                     label = "p.signif", method = "t.test", # wilcox changes?
                     size = 2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Condition", y = "EGFR Mobile Fraction (%)"); egfr_pip2_mf_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_pip2_mf.pdf", egfr_pip2_mf_plot,
       units = "in", height = 2.05, width = 2.36)

egfr_pip3_mf_plot <- ggplot(guv_frap_quant_final_plot %>% 
                              filter(marker == "pEGFR" & slb == "PIP3"),
       aes(x = phCondition, y = mobile_fraction*100)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("+EGFR", "+GRB2/N-WASP"),
                                        c("+EGFR", "+GAB1-WT"),
                                        c("+EGFR", "+GAB1-dPH"),
                                        c("+GRB2/N-WASP", "+GAB1-WT"),
                                        c("+GRB2/N-WASP", "+GAB1-dPH"),
                                        c("+GAB1-WT", "+GAB1-dPH")),
                     label = "p.signif", method = "t.test", # wilcox changes?
                     size = 2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "Condition", y = "EGFR Mobile Fraction (%)"); egfr_pip3_mf_plot
ggsave("manuscript/non_maria_figures/guvs/egfr_pip3_mf.pdf", egfr_pip3_mf_plot,
       units = "in", height = 2.05, width = 2.36)

## visualizing how gab1 mobile fraction changes across conditions 
gab1_mf_plot <- ggplot(guv_frap_quant_final_plot %>% 
                              filter(marker == "GAB1" & phCondition == "+GAB1-WT"),
       aes(x = slb, y = mobile_fraction*100)) +
  geom_boxplot(outlier.size = 0.25, color = "black", linewidth = 0.25,
               outlier.stroke = 0.25) +
  stat_compare_means(comparisons = list(c("POPC", "PIP2"),
                                        c("POPC", "PIP3"),
                                        c("PIP2", "PIP3")),
                     label = "p.signif", method = "t.test", # wilcox changes?
                     size = 2) +
  #scale_y_continuous(trans = "log10") +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black")) +
  labs(x = "GUV Type", y = "GAB1-WT Mobile Fraction (%)"); gab1_mf_plot
ggsave("manuscript/non_maria_figures/guvs/gab1_mf.pdf", gab1_mf_plot,
       units = "in", height = 2.29, width = 2.36)

# gab1 corrected mean frap traces across conditions **FOR VISUALIZATION ONLY**
frap_gab1_plot <- ggplot(guv_frap_final_mean_plot %>%
         filter(marker == "GAB1" & phCondition == "+GAB1-WT") %>%
         filter(mean_frapTime >= 0), 
       aes(x = mean_frapTime, y = mean_cor_y,
                          color = slb)) +
  geom_point(size = 0.25) +
  geom_errorbar(aes(ymin = mean_cor_ymin, ymax = mean_cor_ymax), 
                #alpha = 0.25,
                linewidth = 0.1, width = 0.5) +
  geom_line(linewidth = 0.25) +
  scale_color_manual(values = c("#A9C1E1", "#5382C3", "#3F6293")) +
  theme_classic(base_size = 7) +
  theme(axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        legend.position = "inside",
        legend.position.inside = c(0.1, 0.9)) +
  labs(x = "Time (seconds)", y = "Normalized GAB1-WT FRAP Intensity",
       color = ""); frap_gab1_plot
ggsave("manuscript/non_maria_figures/guvs/frap_gab1.pdf", frap_gab1_plot,
       units = "in", height = 2.29, width = 2.36)


## ----sessionInfo------------------------------------------------------------------------------------------------------------
sessionInfo()

