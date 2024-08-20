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

## Read in Data
hf.mice <- read.csv("./hf.mice.csv")

hf.mice$series <- paste(hf.mice$group," ",hf.mice$sex,"", sep = "")
# Calculate weekly mean
hf.avg <- hf.mice %>%
  group_by(sex, group, week) %>%
  summarise("mean" = mean(weight), "sd" = sd(weight))
# Generate color label for weekly means
hf.avg$series <- paste(hf.avg$group," ",hf.avg$sex,"", sep = "")

# Plot the data using LOESS
ggplot(data = hf.mice, aes(week, weight, group = group, col = group)) +
  geom_smooth() +
  geom_point(data = hf.avg, aes(week, mean, group = group, col = group)) +
  facet_wrap(~sex,
             labeller = as_labeller(c(female = "Females",
                                      male = "Males"))) +
  scale_x_continuous(limits = c(0, 73),
                     breaks = c(4, 24, 48, 72)) +
  scale_color_manual(label = c("HFBN", "HFB", "HFCN", "HFC"),
                     values = c("purple", "chartreuse4", "deepskyblue3", "brown1")) +
  labs(title = "Weekly Change in Total Mass by LOESS",
       y = "Mean Mass (g)", x = "Week", color = "Diets") +
  theme(plot.title = element_text(family = "Fira Sans Condensed", hjust = 0.5, size = 20),
        plot.background = element_rect(fill = "white"),
        panel.border = element_rect(linetype = "solid", fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 18),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.text = element_text(size = 15, color = "black"),
        axis.ticks = element_blank(),
        legend.key = element_blank(),
        legend.position = "right")
