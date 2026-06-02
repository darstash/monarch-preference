# PROJECT: Monarch caterpillar preference
# AUTHORS: Cooper Pryor & Ashley Darst
# DATE: September 24, 2025
# PURPOSE: This script analyzes monarch caterpillar preference for pesticides after initial exposure, correcting for turning bias.

# R version 4.4.1 (2024-06-14)

# Package versions:
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggpubr_0.6.2    car_3.1-3       carData_3.0-5   betareg_3.2-4   emmeans_2.0.1   lubridate_1.9.4
# [7] forcats_1.0.1   stringr_1.6.0   dplyr_1.1.4     purrr_1.2.1     readr_2.1.6     tidyr_1.3.2    
# [13] tibble_3.3.1    ggplot2_4.0.1   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] sandwich_3.1-1     generics_0.1.4     rstatix_0.7.3      stringi_1.8.7      lattice_0.22-7    
# [6] hms_1.1.4          magrittr_2.0.4     grid_4.4.1         estimability_1.5.1 timechange_0.3.0  
# [11] RColorBrewer_1.1-3 mvtnorm_1.3-3      Matrix_1.7-4       backports_1.5.0    nnet_7.3-20       
# [16] Formula_1.2-5      survival_3.8-6     multcomp_1.4-29    scales_1.4.0       TH.data_1.1-5     
# [21] modeltools_0.2-24  codetools_0.2-20   abind_1.4-8        cli_3.6.5          rlang_1.1.7       
# [26] splines_4.4.1      withr_3.0.2        flexmix_2.3-20     tools_4.4.1        tzdb_0.5.0        
# [31] ggsignif_0.6.4     coda_0.19-4.1      broom_1.0.12       vctrs_0.7.1        R6_2.6.1          
# [36] stats4_4.4.1       zoo_1.8-15         lifecycle_1.0.5    MASS_7.3-65        pkgconfig_2.0.3   
# [41] pillar_1.11.1      gtable_0.3.6       glue_1.8.0         lmtest_0.9-40      tidyselect_1.2.1  
# [46] rstudioapi_0.18.0  farver_2.1.2       xtable_1.8-4       labeling_0.4.3     compiler_4.4.1    
# [51] S7_0.2.1

# Read in data
# Main caterpillar datasheet
cat_m_data <- read.csv("monarch_datasheet.csv")
# Surfactant test datasheet
cat_surf_data <- read.csv("silwet_caterpillar_datasheet.csv")

# Load libraries
library(tidyverse)
library(emmeans)
library(betareg)
library(car)
library(ggpubr)

# Data cleaning ----

# Check dataset
View(cat_m_data)
str(cat_m_data) # Check variable types

# Check surfactant dataset
View(cat_surf_data)
str(cat_surf_data) # Check variable types

# Convert dates to date format
cat_m_data$setup_date <- as.Date(cat_m_data$setup_date)
cat_m_data$date_hatched <- as.Date(cat_m_data$date_hatched)
cat_m_data$date_assay <- as.Date(cat_m_data$date_assay)
cat_surf_data$date_assay <- as.Date(cat_surf_data$date_assay)
cat_m_data$food_treatment <- as.factor(cat_m_data$food_treatment)

# Exclude assay in which larva died from surfactant dataset
cat_surf_data <- cat_surf_data[cat_surf_data$caterpillar_id != "m_sw_1", ]

# Calculations and analyses ----

# Adding columns representing the area of each leaf disc eaten and the preference indices.
cat_m_data <- cat_m_data %>%
  mutate(
    # Area eaten for each leaf disc.
    cont_area_eaten = ifelse(control_area_remaining > 1.76714586764, 0, 1.76714586764 - control_area_remaining),
    pest_area_eaten = ifelse(pesticide_area_remaining > 1.76714586764, 0, 1.76714586764 - pesticide_area_remaining),
    # Pesticide preference index using proportion.
    pest_pref = ((pest_area_eaten) / (cont_area_eaten + pest_area_eaten)),
    # Directional preference index using proportion.
    lr_pref = ((ifelse(cat_p_lr == "left", cont_area_eaten, pest_area_eaten)) / (cont_area_eaten + pest_area_eaten)),
    # Area of left disc eaten.
    left_eaten = (ifelse(cat_p_lr == "left", pest_area_eaten, cont_area_eaten)),
    # Area of right disc eaten.
    right_eaten = (ifelse(cat_p_lr == "right", pest_area_eaten, cont_area_eaten))
  )

# Beta regression is bounded at (0, 1) so need to transform data to prevent zeros and ones using the formula in the betareg documentation
cat_m_data$pest_pref_trans <- (cat_m_data$pest_pref * (21 - 1) + 0.5) / 21

# Re-run model with transformed response value
data_model_trans <- betareg(pest_pref_trans ~ cat_p_lr + food_treatment, data = cat_m_data, link = "logit")

summary(data_model_trans)

# Check model diagnostics
par(mfrow = c(3, 2))
set.seed(123)
plot(data_model_trans, which = 1:4, type = "pearson")
plot(data_model_trans, which = 5, type = "deviance", sub.caption = "")
plot(data_model_trans, which = 1, type = "deviance", sub.caption = "")

# Try a log-log family to help with extremes
data_model_log <- betareg(pest_pref_trans ~ cat_p_lr + food_treatment, data = cat_m_data, link = "loglog")
summary(data_model_log)
emm <- emmeans(data_model_log, pairwise ~ food_treatment, infer = T)
summary(emm)
Anova(data_model_log, type = "II") # can safely ignore these warnings
em <- emmeans(data_model_log, ~ food_treatment)
contrast(em, method = "identity", null = 0.5)

# Check model diagnostics
par(mfrow = c(3, 2))
set.seed(123)
plot(data_model_log, which = 1:4, type = "pearson")
plot(data_model_log, which = 5, type = "deviance", sub.caption = "")
plot(data_model_log, which = 1, type = "deviance", sub.caption = "")
plot(data_model_log, which = 1:6)

# Is the model improved with logit or log-log link functions?
summary(data_model_trans)$pseudo.r.squared
summary(data_model_log)$pseudo.r.squared # Better R-squared

AIC(data_model_trans, data_model_log) # AIC not different (need at least a difference of 2)

# Calculating the difference between left and right leaf discs eaten for Shapiro-Wilk test.
lr_diff <- cat_m_data$left_eaten - cat_m_data$right_eaten

# Shapiro-Wilk test to determine normality.
shapiro.test(lr_diff) # go with normal distribution

# t-test to determine if there is a significant difference between left and right discs eaten.
t.test(cat_m_data$left_eaten, cat_m_data$right_eaten, paired = TRUE, alternative = "two.sided")

# Adding columns to the surfactant datasheet representing area eaten and preference index.
cat_surf_data <- cat_surf_data %>%
  mutate(
    # Area eaten for each leaf disc.
    water_area_eaten = ifelse(end_water_area > 1.76714586764, 0, 1.76714586764 - end_water_area),
    surf_area_eaten = ifelse(end_silwet_area > 1.76714586764, 0, 1.76714586764 - end_silwet_area),
    # Surfactant preference index using proportion.
    surf_pref = ((surf_area_eaten) / (water_area_eaten + surf_area_eaten)))

# Calculating the difference between surfactant treated and water treated discs for a Shapiro-Wilk test.
surf_diff <- cat_surf_data$surf_area_eaten - cat_surf_data$water_area_eaten

# Shapiro-Wilk test to determine normality.
shapiro.test(surf_diff) # go with normal distribution

# t-test to determine if there is a significant difference between surfactant treated and water treated discs eaten.
t.test(cat_surf_data$surf_area_eaten, cat_surf_data$water_area_eaten, paired = TRUE, alternative = "two.sided")

# Plots ----

## Figure 2 ----
# pdf("plot_fig_2.pdf", width = 5, height = 5)
ggsave("plot_fig_2.tiff", width = 5, height = 5, units = "in", dpi = 600)
cat_m_data %>%
  ggplot(aes(x = food_treatment, y = pest_pref_trans)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.2) +
  stat_summary(fun.data = "mean_se") +
  theme_classic() +
  labs(x = "Initial Food Treatment", y = "Preference Index") +
  scale_x_discrete(labels = c("control" = "Control", 
                              "pesticide" = "Pesticide")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 14))
dev.off()

## Figure 3 ----
# Plot of disc placement and area eaten.
d <- data.frame(Left = cat_m_data$left_eaten, Right = cat_m_data$right_eaten)
g <- ggpaired(d, cond1 = "Left", cond2 = "Right",
         fill = "condition",
         xlab = "Disc Placement",
         palette = c("gray90","gray63"))
# pdf("plot_fig_3.pdf", width = 5, height = 5)
ggsave("plot_fig_3.tiff", width = 5, height = 5, units = "in", dpi = 600)
g + labs(y = bquote('Area Eaten'~(cm^2))) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 14))
dev.off()

##Figure S1 ----
# Plot of silwet and area eaten.
s <- data.frame(control = cat_surf_data$water_area_eaten, surfactant = cat_surf_data$surf_area_eaten)
silwet.plot <- ggpaired(s, cond1 = "control", cond2 = "surfactant",
              fill = "condition",
              xlab = "Treatment",
              palette = c("gray90","gray63")) +
  scale_x_discrete(labels = c("control" = "Control", 
                              "surfactant" = "Surfactant"))
jpeg("plot_fig_s1.jpg", width = 5, height = 5, units = "in", res = 300)
silwet.plot + labs(y = bquote('Area Eaten'~(cm^2))) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 14))
dev.off()
