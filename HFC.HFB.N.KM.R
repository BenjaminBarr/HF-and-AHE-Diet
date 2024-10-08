install.packages("ggdark")
install.packages("ggplot2")
install.packages("readxl")
install.packages("reshape2")
install.packages("matrixStats")
install.packages("dplyr")
install.packages("AICcmodavg")
install.packages("tidyverse")
install.packages("ggpubr")
install.packages("rstatix")
install.packages("emmeans")
install.packages("lme4")


#install.packages(c("lubridate", "ggsurvfit", "gtsummary", "tidycmprsk"))
library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

#devtools::install_github("zabore/condsurv")
#library(condsurv)

# Data Read In

library(ggplot2)
library(readxl)
library(reshape2)
library(matrixStats)
library(dplyr)
library(AICcmodavg)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(ggdark)
library(emmeans)
library(lme4)
library(survival)
library(knitr)
library(tibble)
library(cowplot)
library(ggtext)

#READ ME#
## Shorthand for groups in this assessment are as follows:
  # 3 = HFBN
  # 4 = HFB
  # 5 = HFCN
  # 6 = HFC

#Event data (Censored[0] and Deaths[1])
getwd()
hflive <- read.csv("./Data/hf.live.csv")

names(hflive)

# Generate Kaplan-Meier Model
hf.km.model <- survfit(data = hflive, Surv(week, death) ~ group)

### ALL ###
hflive$group <- factor(hflive$group, levels = c("group 6", "group 5", "group 4", "group 3"))
#hflive$group <- relevel(hflive, ref = "group 3")
female.plot <- survfit2(Surv(week, death) ~ group, data = hflive[which(hflive$sex == "female"),]) %>%
  ggsurvfit(linewidth = 2, linetype_aes = T) +
  scale_color_manual(values = c("brown1", "deepskyblue3", "chartreuse4", "purple"), 
                     labels = c("HFC", "HFCN", "HFB", "HFBN"), name = "Diets: ") +
  labs(
    x = "Weeks",
    y = "Overall Survival Probability",
    title = "Females") +
  guides(linetype = "none") +
  theme_ggsurvfit_KMunicate() +  
  geom_richtext(fill = NA, label.color = NA, y = .50, x = 20, 
                label = "Log-rank *p* = 0.5", 
                size = 5, family = "sans") +
  scale_ggsurvfit(y_scales = list(limits = c(.4, 1), breaks = seq(.2, 1, by = .2)))
female.plot

## Males ##
male.plot <- survfit2(Surv(week, death) ~ group, data = hflive[which(hflive$sex == "male"),]) %>%
  ggsurvfit(linewidth = 2, linetype_aes = T) +
  scale_color_manual(values = c("brown1", "deepskyblue3", "chartreuse4", "purple"), 
                     labels = c("HFC", "HFCN", "HFB", "HFBN"), name = "Diets: ") +
  labs(
    x = "Weeks",
    y = "Overall Survival Probability",
    title = "Males") +
  guides(linetype = "none") +
  theme_ggsurvfit_KMunicate() +  
  geom_richtext(fill = NA, label.color = NA, y = .50, x = 20, 
                                               label = "Log-rank *p* = 0.066", 
                                               size = 5, family = "sans") +
  scale_ggsurvfit(y_scales = list(limits = c(.4, 1), breaks = seq(.2, 1, by = .2)))

mfkm <- female.plot + male.plot + 
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

title <- ggdraw() + 
  draw_label(
    "C3H/HeJ Survivability Over 18-Months",
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )
# Add title to female and male KM
plot_grid(title, mfkm,
          ncol = 1,
          rel_heights = c(0.1, 1))

# Assess significance using survdiff
survdiff(formula = Surv(week, death) ~ group, data = hflive[which(hflive$sex == "male"),])
survdiff(formula = Surv(week, death) ~ group, data = hflive[which(hflive$sex == "female"),])

coxph(formula = Surv(week, death) ~ group, data = hflive[which(hflive$sex == "male"),]) %>%
  tbl_regression(exp = T)

summary(hf.km.model)
