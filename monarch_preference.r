# PROJECT: Monarch caterpillar preference
# AUTHORS: Cooper Pryor & Ashley Darst
# DATE: September 24, 2025
# REVISED: September 24, 2026
# PURPOSE: This script analyzes monarch caterpillar preference for pesticides after initial exposure, correcting for turning bias.

# R version 4.6.1 (2026-06-24) -- "Happy Hop"
# Platform: x86_64-apple-darwin20
# Running under: macOS Sequoia 15.7.9

# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] DHARMa_0.5.0    glmmTMB_1.1.14  ggpubr_1.0.0    car_3.1-5       carData_3.0-6   emmeans_2.0.3  
# [7] lubridate_1.9.5 forcats_1.0.1   stringr_1.6.0   dplyr_1.2.1     purrr_1.2.2     readr_2.2.0    
# [13] tidyr_1.3.2     tibble_3.3.1    ggplot2_4.0.3   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.2.1    farver_2.1.2        S7_0.2.2            fastmap_1.2.0       TH.data_1.1-5      
# [6] promises_1.5.0      digest_0.6.39       mime_0.13           estimability_2.0.0  timechange_0.4.0   
# [11] lifecycle_1.0.5     survival_3.8-6      magrittr_2.0.5      compiler_4.6.1      rlang_1.3.0        
# [16] tools_4.6.1         ggsignif_0.6.4      labeling_0.4.3      plyr_1.8.9          RColorBrewer_1.1-3 
# [21] gap.datasets_0.0.6  multcomp_1.4-31     abind_1.4-8         withr_3.0.3         numDeriv_2016.8-1.1
# [26] grid_4.6.1          xtable_1.8-8        scales_1.4.0        iterators_1.0.14    MASS_7.3-65        
# [31] cli_3.6.6           mvtnorm_1.4-2       ragg_1.5.2          reformulas_0.4.4    generics_0.1.4     
# [36] otel_0.2.0          rstudioapi_0.19.0   tzdb_0.5.0          minqa_1.2.8         splines_4.6.1      
# [41] parallel_4.6.1      vctrs_0.7.3         boot_1.3-32         Matrix_1.7-5        sandwich_3.1-2     
# [46] hms_1.1.4           rstatix_1.0.0       qgam_2.0.0          Formula_1.2-5       systemfonts_1.3.2  
# [51] foreach_1.5.2       gap_1.15.2          glue_1.8.1          nloptr_2.2.1        codetools_0.2-20   
# [56] stringi_1.8.7       gtable_0.3.6        later_1.4.8         lme4_2.0-1          pillar_1.11.1      
# [61] htmltools_0.5.9     R6_2.6.1            TMB_1.9.21          textshaping_1.0.5   Rdpack_2.6.6       
# [66] doParallel_1.0.17   shiny_1.14.0        lattice_0.22-9      rbibutils_2.4.1     backports_1.5.1    
# [71] broom_1.0.13        httpuv_1.6.17       Rcpp_1.1.2          nlme_3.1-169        mgcv_1.9-4         
# [76] zoo_1.8-15          pkgconfig_2.0.3   

# Read in data
# Main caterpillar datasheet
cat_m_data <- read.csv("monarch_datasheet.csv")
# Surfactant test datasheet
cat_surf_data <- read.csv("silwet_caterpillar_datasheet.csv")

# Load libraries
library(tidyverse)
library(emmeans)
library(car)
library(ggpubr)
library(glmmTMB)
library(DHARMa)

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

# Use ordbeta to capture true 0s and 1s
data_model_ord <- glmmTMB(pest_pref ~ cat_p_lr + food_treatment, data = cat_m_data, family = ordbeta)
summary(data_model_ord)
Anova(data_model_ord, type = "II")
# Estimated marginal means on the response scale
emmeans(data_model_ord, ~ food_treatment, type = "response", infer = T)
# On the link scale, 0.5 = 0
emmeans(data_model_ord, ~ food_treatment, infer = T, null = 0)

# Check model
res <- simulateResiduals(data_model_ord)
plot(res) # not a major concern the residual test is signif because the individual variables are NS
testUniformity(res)
testDispersion(res)
testOutliers(res)
plotResiduals(res, cat_m_data$cat_p_lr)
plotResiduals(res, cat_m_data$food_treatment)

# Calculating the difference between control and pesticide leaf discs eaten
mean(cat_m_data$cont_area_eaten - cat_m_data$pest_area_eaten)
0.3225719 / 1.767145

# Wilcoxon signed-rank test due to small sample size (double check, but stick with ordbeta results)
wilcox.test(cat_m_data$pest_pref, mu = 0.5, alternative = "two.sided")

# Calculating the difference between left and right leaf discs eaten
mean(cat_m_data$right_eaten - cat_m_data$left_eaten)
0.6918156 / 1.767145

# Wilcoxon signed-rank test due to small sample size
wilcox.test(cat_m_data$lr_pref, mu = 0.5, alternative = "two.sided")
median(cat_m_data$lr_pref)
quantile(cat_m_data$lr_pref, c(0.25, 0.75))

# Adding columns to the surfactant datasheet representing area eaten and preference index.
cat_surf_data <- cat_surf_data %>%
  mutate(
    # Area eaten for each leaf disc.
    water_area_eaten = ifelse(end_water_area > 1.76714586764, 0, 1.76714586764 - end_water_area),
    surf_area_eaten = ifelse(end_silwet_area > 1.76714586764, 0, 1.76714586764 - end_silwet_area),
    # Surfactant preference index using proportion.
    surf_pref = ((surf_area_eaten) / (water_area_eaten + surf_area_eaten)))

# Wilcoxon signed-rank test due to small sample size
wilcox.test(cat_surf_data$surf_pref, mu = 0.5, alternative = "two.sided")
median(cat_surf_data$surf_pref)
quantile(cat_surf_data$surf_pref, c(0.25, 0.75))

# Adding columns representing the area of each leaf disc eaten and the preference indices.
cat_surf_data <- cat_surf_data %>%
  mutate(
    # Directional preference index using proportion.
    lr_pref = ((ifelse(cat_sw_lr == "left", water_area_eaten, surf_area_eaten)) / (water_area_eaten + surf_area_eaten)),
    # Area of left disc eaten.
    left_eaten = (ifelse(cat_sw_lr == "left", surf_area_eaten, water_area_eaten)),
    # Area of right disc eaten.
    right_eaten = (ifelse(cat_sw_lr == "right", surf_area_eaten, water_area_eaten))
  )

wilcox.test(cat_surf_data$lr_pref, mu = 0.5, alternative = "two.sided")
median(cat_surf_data$surf_pref)
quantile(cat_surf_data$surf_pref, c(0.25, 0.75))

# Ordbeta for surfactant preference left/right
data_surf_ord <- glmmTMB(surf_pref ~ cat_sw_lr, data = cat_surf_data, family = ordbeta)
summary(data_surf_ord)
Anova(data_surf_ord, type = "II")
emmeans(data_surf_ord, specs = ~ 1, type = "response", infer = T)
emmeans(data_surf_ord, specs = ~ 1, infer = T)
# Estimated marginal means on the response scale
emmeans(data_surf_ord, ~ cat_sw_lr, type = "response", infer = T)
emmeans(data_surf_ord, ~ cat_sw_lr, infer = T)

## Leaf disc weight

# Initial feeding treatment
cat_m_data %>%
  summarize(mean = mean(disc_weight),
            min = min(disc_weight),
            max = max(disc_weight))

# Preference assay
cat_m_data %>%
  summarise(mean_value = mean(c(start_control_weight, start_pesticide_weight), na.rm = TRUE),
            min_value = min(c(start_control_weight, start_pesticide_weight), na.rm = TRUE),
            max_value = max(c(start_control_weight, start_pesticide_weight), na.rm = TRUE))

# Silwet preference assay
cat_surf_data %>%
  summarise(mean_value = mean(c(start_water_weight, start_silwet_weight), na.rm = TRUE),
            min_value = min(c(start_water_weight, start_silwet_weight), na.rm = TRUE),
            max_value = max(c(start_water_weight, start_silwet_weight), na.rm = TRUE))

# Plots ----

## Figure 2a ----
# pdf("plot_fig_2.pdf", width = 5, height = 5)
ggsave("plot_fig_2.tiff", width = 5, height = 5, units = "in", dpi = 600)
cat_m_data %>%
  ggplot(aes(x = food_treatment, y = pest_pref)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.2) +
  stat_summary(fun.data = "mean_se") +
  theme_classic() +
  labs(x = "Initial Food Treatment", y = "Preference Index") +
  scale_x_discrete(labels = c("control" = "Control", 
                              "pesticide" = "Pesticide")) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  theme(text = element_text(size = 14))
dev.off()

## Figure 2b ----
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

##Figure S2 ----
# Plot of preference assay model assumptions.
jpeg("plot_fig_s2.jpg", width = 5, height = 5, units = "in", res = 300)
testUniformity(res)
dev.off()

