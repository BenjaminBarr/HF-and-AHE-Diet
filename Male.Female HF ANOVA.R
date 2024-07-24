install.packages('ggpattern')
install.packages("patchwork")

#Three Way interactions
install.packages(c("sjPlot", "sjmisc"))
install.packages("ggeffects")


#Three-Way Interactions
library(sjPlot)
library(sjmisc)
library(effects)
library(ggeffects)

#
library(here) # here makes a project transportable
library(janitor) # clean_names
library(readxl) # read excel, duh!
library(data.table) # magical data frames
library(magrittr) # pipes
library(stringr) # string functions
library(forcats) # factor functions

# analysis packages
library(emmeans) # the workhorse for inference
library(nlme) # gls and some lmm
library(lme4) # linear mixed models
library(lmerTest) # linear mixed model inference
library(afex) # ANOVA linear models
library(glmmTMB) # generalized linear models
library(MASS) # negative binomial and some other functions
library(car) # model checking and ANOVA
library(DHARMa) # model checking

# graphing packages
library(ggsci) # color palettes
library(ggpubr) # publication quality plots
library(ggforce) # better jitter

library(knitr) # kable tables
library(kableExtra) # kable_styling tables

library(insight)
library(lazyWeave)

library(cowplot) # combine plots
library(patchwork)# combine plots
# Data Read In

library(reshape2)
library(matrixStats)
library(dplyr)
library(AICcmodavg)
library(rstatix)
library(ggdark)
library(multcompView)
library(tidyverse)
#STRIPED BARS
library(ggpattern)
##Functions for normalization
source_path <- here("raw-R work", "ggplot_the_model.R")
source(source_path)

## Data Start ##
getwd()

hfdt <- read.csv("./Data/HFandN MRI.csv")

hf.summary <- hfdt %>% group_by(group, sex, m.start,) %>%
  summarise(mean_mass = mean(weight),
            se_mass = sd(weight)/sqrt(length(weight)),
            mean_fat = mean(fat),
            se_fat = sd(fat)/sqrt(length(fat)),
            mean_lean = mean(lean),
            se_lean = sd(lean)/sqrt(length(lean)),
            n = n())
hf.summary

write.csv(hf.summary, "./Data/HFSummary.csv")

hf6.lm <- lm(weight ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 6),])
hf12.lm <- lm(weight ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 12),])
hf18.lm <- lm(weight ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 18),])

anova(hf6.lm)
anova(hf12.lm)
anova(hf18.lm)

anova(hf6.lm)

hf6.em <- emmeans(hf6.lm, specs = c("sex", "ph", "protein"))
hf6.em

hf6.em_dt <- data.table(summary(hf6.em))
hf6.em_dt

hf6.pairs <- contrast(hf6.em,
                     method = "revpairwise",
                     simple = "each",
                     combine = T,
                     adjust = "sidak") %>%
  summary(infer = T, )

hf6.pairs <- na.omit(hf6.pairs)
hf6.pairs_dt <- data.table(hf6.pairs)

group1 <- hf6.pairs_dt$contrast

group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))
group1 <- as.factor(group1)

group2 <- hf6.pairs_dt$contrast
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))
group2 <- as.factor(group2)
group2

hf6.pairs_dt$group1 <- group1
hf6.pairs_dt$group2 <- group2
hf6.pairs_dt


hf6.pairs_dt[, p_rounded := p_round(p.value,
                                   digits = 2)]
hf6.pairs_dt[, p_pretty := p_format(p_rounded,
                                   digits = 2,
                                   accuracy = 1e-04,
                                   add.p = TRUE)]
hf6.pairs_dt

levels(hf6.pairs)

nrow(hf6.pairs_dt)

hf6.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                      "group 4", "group 4", "group 3", "group 3")

hf6.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                          "group 5", "group 5", "group 3", "group 3",
                          "group 4", "group 4", "group 3", "group 3")


hf6.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                          "group 6", "group 6", "group 4", "group 4",
                          "group 6", "group 6", "group 5", "group 5")


ggbarplot(data = hf6.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap(facets = "sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf6.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf6.pairs_dt[-c(1:4)],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 50, to = 64, by = 2)),
                     tip.length = 0.01) +
  labs(title = "HF Month 6 Mass", y ="Weight")

### HF 12 Mass
hf12.em <- emmeans(hf12.lm, specs = c("sex", "ph", "protein"))
hf12.em

hf12.em_dt <- data.table(summary(hf12.em))
hf12.em_dt

hf12.pairs <- contrast(hf12.em,
                       method = "revpairwise",
                       simple = "each",
                       combine = T,
                       adjust = "sidak") %>%
  summary(infer = T, )

hf12.pairs <- na.omit(hf12.pairs)
hf12.pairs_dt <- data.table(hf12.pairs)

group1 <- hf12.pairs_dt$contrast

group1 <- gsub(" - fatc0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))
group1 <- as.factor(group1)

group2 <- hf12.pairs_dt$contrast
group2 <- gsub("fatc1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))
group2 <- as.factor(group2)
group2

hf12.pairs_dt$group1 <- group1
hf12.pairs_dt$group2 <- group2
hf12.pairs_dt


hf12.pairs_dt[, p_rounded := p_round(p.value,
                                     digits = 2)]
hf12.pairs_dt[, p_pretty := p_format(p_rounded,
                                     digits = 2,
                                     accuracy = 1e-04,
                                     add.p = TRUE)]
hf12.pairs_dt

nrow(hf12.pairs_dt)

hf12.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                      "group 4", "group 4", "group 3", "group 3")

hf12.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                          "group 5", "group 5", "group 3", "group 3",
                          "group 4", "group 4", "group 3", "group 3")


hf12.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                          "group 6", "group 6", "group 4", "group 4",
                          "group 6", "group 6", "group 5", "group 5")


ggbarplot(data = hf12.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap(facets = "sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf12.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf12.pairs_dt[-c(1:4)],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 50, to = 64, by = 2)),
                     tip.length = 0.01) +
  labs(title = "HF Month 12 Mass", y ="Weight")

#### 18 HF Mass
hf18.em <- emmeans(hf18.lm, specs = c("sex", "ph", "protein"))
hf18.em

hf18.em_dt <- data.table(summary(hf18.em))
hf18.em_dt

hf18.pairs <- contrast(hf18.em,
                       method = "revpairwise",
                       simple = "each",
                       combine = T,
                       adjust = "sidak") %>%
  summary(infer = T, )

hf18.pairs <- na.omit(hf18.pairs)
hf18.pairs_dt <- data.table(hf18.pairs)

group1 <- hf18.pairs_dt$contrast

group1 <- gsub(" - fatc0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))
group1 <- as.factor(group1)

group2 <- hf18.pairs_dt$contrast
group2 <- gsub("fatc1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))
group2 <- as.factor(group2)
group2

hf18.pairs_dt$group1 <- group1
hf18.pairs_dt$group2 <- group2
hf18.pairs_dt


hf18.pairs_dt[, p_rounded := p_round(p.value,
                                     digits = 2)]
hf18.pairs_dt[, p_pretty := p_format(p_rounded,
                                     digits = 2,
                                     accuracy = 1e-04,
                                     add.p = TRUE)]

nrow(hf18.pairs_dt)

hf18.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                      "group 4", "group 4", "group 3", "group 3")

hf18.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                          "group 5", "group 5", "group 3", "group 3",
                          "group 4", "group 4", "group 3", "group 3")


hf18.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                          "group 6", "group 6", "group 4", "group 4",
                          "group 6", "group 6", "group 5", "group 5")


ggbarplot(data = hf18.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap(facets = "sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf18.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf18.pairs_dt[-c(1:4)],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 50, to = 64, by = 2)),
                     tip.length = 0.01) +
  labs(title = "HF Month 18 Mass", y ="Weight")

## ALL MASS ##
hf.em_dt <- rbind(hf6.em_dt, hf12.em_dt, hf18.em_dt)
hf.em_dt$age <- rep(c(6, 12, 18), each = 8)

hf6.pairs_dt$age <- 6
hf12.pairs_dt$age <- 12
hf18.pairs_dt$age <- 18
hf.pairs <- rbind(hf6.pairs_dt, hf12.pairs_dt, hf18.pairs_dt)
hf.pairs

hf.em_dt$xlab <- paste(hf.em_dt$group," ",hf.em_dt$age,"", sep = "")
hf.pairs$xlab1 <- paste(hf.pairs$group1, hf.pairs$age)
hf.pairs$xlab2 <- paste(hf.pairs$group2, hf.pairs$age)

hf.pairs.sig <- hf.pairs[which(hf.pairs$p.value <= 0.05),]
hf.em_dt <- hf.em_dt[order(hf.em_dt$group),]


## Anova significance ##


hf6 <- hfdt[which(hfdt$m.start == 6),]
hf12 <- hfdt[which(hfdt$m.start == 12),]
hf18 <- hfdt[which(hfdt$m.start == 18),]

anova(hf18.lm) ##No significant interactions...


#Facet Labels
sex_label <- c("Female", "Male")
names(sex_label) <- c(0, 1)

## Order the Graphs ##

hf.em_dt$group <- as.factor(hf.em_dt$group)
hf.em_dt <- hf.em_dt[order(hf.em_dt$group, decreasing = TRUE),]

#All mass
ggbarplot(data = hf.em_dt, 
          x = "xlab", 
          y = "emmean", 
          fill = "group") +
  facet_wrap(facets = "sex", labeller = labeller(sex = sex_label)) +
  theme_pubr() +
  geom_errorbar(data = hf.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE), width = 0.1, color = "black") +
  scale_fill_manual(values = c("brown1", "deepskyblue3", "chartreuse4", "purple"),
                    breaks = c("group 6", "group 5", "group 4", "group 3"),
                    labels = c("HFC", "HFCN", "HFB", "HFBN"),
                    name = "Diets") +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 12)) +
  labs(title = "Mass Over Time", y ="Total Mass (g)", x = "Collection Timepoint")


##HF Fat

hf6f.lm <- lm(fat ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 6),])
hf12f.lm <- lm(fat ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 12),])
hf18f.lm <- lm(fat ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 18),])

anova(hf6f.lm)
anova(hf12f.lm)
anova(hf18f.lm)

hf6f.em <- emmeans(hf6f.lm, specs = c("sex", "ph", "protein"))
hf6f.em

hf6f.em_dt <- data.table(summary(hf6f.em))
hf6f.em_dt

hf6f.pairs <- contrast(hf6f.em,
                      method = "revpairwise",
                      simple = "each",
                      combine = T,
                      adjust = "sidak") %>%
  summary(infer = T, )

hf6f.pairs <- na.omit(hf6f.pairs)
hf6f.pairs_dt <- data.table(hf6f.pairs)

group1 <- hf6f.pairs_dt$contrast

group1 <- gsub(" - protein0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - ph0", "", as.character(group1))

group1 <- as.factor(group1)

group2 <- hf6f.pairs_dt$contrast
group2 <- gsub("protein1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("ph1 - ", "", as.character(group2))

group2 <- as.factor(group2)
group2

hf6f.pairs_dt$group1 <- group1
hf6f.pairs_dt$group2 <- group2
hf6f.pairs_dt


hf6f.pairs_dt[, p_rounded := p_round(p.value,
                                    digits = 2)]
hf6f.pairs_dt[, p_pretty := p_format(p_rounded,
                                    digits = 2,
                                    accuracy = 1e-04,
                                    add.p = TRUE)]
hf6f.pairs_dt

hf6f.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                       "group 4", "group 4", "group 3", "group 3")

hf6f.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                           "group 5", "group 5", "group 3", "group 3",
                           "group 4", "group 4", "group 3", "group 3")


hf6f.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                           "group 6", "group 6", "group 4", "group 4",
                           "group 6", "group 6", "group 5", "group 5")

ggbarplot(data = hf6f.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf6f.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf6f.pairs_dt[-c(1:4),],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 30, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 6 Fat Mass", y ="Weight")

### C 12 Fat
anova(hf12f.lm)

hf12f.em <- emmeans(hf12f.lm, specs = c("sex", "ph", "protein"))
hf12f.em

hf12f.em_dt <- data.table(summary(hf12f.em))
hf12f.em_dt

hf12f.pairs <- contrast(hf12f.em,
                       method = "revpairwise",
                       simple = "each",
                       combine = T,
                       adjust = "sidak") %>%
  summary(infer = T, )

hf12f.pairs <- na.omit(hf12f.pairs)
hf12f.pairs_dt <- data.table(hf12f.pairs)

group1 <- hf12f.pairs_dt$contrast

group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))

group1 <- as.factor(group1)

group2 <- hf12f.pairs_dt$contrast
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))

group2 <- as.factor(group2)
group2

hf12f.pairs_dt$group1 <- group1
hf12f.pairs_dt$group2 <- group2
hf12f.pairs_dt


hf12f.pairs_dt[, p_rounded := p_round(p.value,
                                     digits = 2)]
hf12f.pairs_dt[, p_pretty := p_format(p_rounded,
                                     digits = 2,
                                     accuracy = 1e-04,
                                     add.p = TRUE)]
hf12f.pairs_dt

hf12f.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                       "group 4", "group 4", "group 3", "group 3")

hf12f.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                           "group 5", "group 5", "group 3", "group 3",
                           "group 4", "group 4", "group 3", "group 3")


hf12f.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                           "group 6", "group 6", "group 4", "group 4",
                           "group 6", "group 6", "group 5", "group 5")

ggbarplot(data = hf12f.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf12f.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf12f.pairs_dt[-c(1:4),],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 30, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 12 Mass", y ="Weight")

#### 18 C Fat
anova(hf18f.lm)

hf18f.em <- emmeans(hf18f.lm, specs = c("sex", "ph", "protein"))
hf18f.em

hf18f.em_dt <- data.table(summary(hf18f.em))
hf18f.em_dt

hf18f.pairs <- contrast(hf18f.em,
                       method = "revpairwise",
                       simple = "each",
                       combine = T,
                       adjust = "sidak") %>%
  summary(infer = T, )

hf18f.pairs <- na.omit(hf18f.pairs)
hf18f.pairs_dt <- data.table(hf18f.pairs)

group1 <- hf18f.pairs_dt$contrast

group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))

group1 <- as.factor(group1)

group2 <- hf18f.pairs_dt$contrast
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))

group2 <- as.factor(group2)
group2

hf18f.pairs_dt$group1 <- group1
hf18f.pairs_dt$group2 <- group2
hf18f.pairs_dt


hf18f.pairs_dt[, p_rounded := p_round(p.value,
                                     digits = 2)]
hf18f.pairs_dt[, p_pretty := p_format(p_rounded,
                                     digits = 2,
                                     accuracy = 1e-04,
                                     add.p = TRUE)]
hf18f.pairs_dt

hf18f.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                        "group 4", "group 4", "group 3", "group 3")

hf18f.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                            "group 5", "group 5", "group 3", "group 3",
                            "group 4", "group 4", "group 3", "group 3")


hf18f.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                            "group 6", "group 6", "group 4", "group 4",
                            "group 6", "group 6", "group 5", "group 5")

ggbarplot(data = hf18f.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf18f.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf18f.pairs_dt[-c(1:4),],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 30, to = 44, by = 2)),
                     tip.length = 0.01) +
  labs(title = "Month 18 Fat Mass", y ="Weight")

## ALL FAT ##
f.em_dt <- rbind(hf6f.em_dt, hf12f.em_dt, hf18f.em_dt)
f.em_dt$age <- rep(c(6, 12, 18), each = 8)

hf6f.pairs_dt$age <- 6
hf12f.pairs_dt$age <- 12
hf18f.pairs_dt$age <- 18
f.pairs <- rbind(hf6f.pairs_dt, hf12f.pairs_dt, hf18f.pairs_dt)
f.pairs

f.em_dt$xlab <- paste(f.em_dt$group," ",f.em_dt$age,"", sep = "")
f.pairs$xlab1 <- paste(f.pairs$group1, f.pairs$age)
f.pairs$xlab2 <- paste(f.pairs$group2, f.pairs$age)

f.pairs.sig <- f.pairs[which(f.pairs$p.value <= 0.05),]
f.em_dt <- f.em_dt[order(f.em_dt$group),]

nrow(f.pairs)

## Order the graphs ##

f.em_dt$group <- as.factor(f.em_dt$group)
f.em_dt <- f.em_dt[order(f.em_dt$group, decreasing = T),]

#All Fat Graph
ggbarplot(data = f.em_dt, 
          x = "xlab", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex", labeller = labeller(sex = sex_label))+
  theme_pubr() +
  geom_errorbar(data = f.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE), width = 0.1, color = "black") +
  scale_fill_manual(values = c("brown1", "deepskyblue3", "chartreuse4", "purple"),
                    breaks = c("group 6", "group 5", "group 4", "group 3"),
                    labels = c("HFC", "HFCN", "HFB", "HFBN"),
                    name = "Diets") +
  stat_pvalue_manual(data = f.pairs.sig[-c(1:2),],
                     label = "p_pretty",
                     xmin = "xlab1",
                     xmax = "xlab2",
                     y.position = c(seq(from = 16, to = 20, by = 2)),
                     tip.length = 0.01) +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 4)) +
  labs(title = "Fat Mass Over Time", y ="Fat Mass (g)", x = "Collection Timepoint (Months)")

##1st Lean

hf6l.lm <- lm(lean ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 6),])
hf12l.lm <- lm(lean ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 12),])
hf18l.lm <- lm(lean ~ sex * ph * protein, data = hfdt[which(hfdt$m.start == 18),])

anova(hf6l.lm)
anova(hf12l.lm)
anova(hf18l.lm)

hf6l.em <- emmeans(hf6l.lm, specs = c("sex", "ph", "protein"))
hf6l.em

hf6l.em_dt <- data.table(summary(hf6l.em))
hf6l.em_dt

hf6l.pairs <- contrast(hf6l.em,
                        method = "revpairwise",
                        simple = "each",
                        combine = T,
                        adjust = "sidak") %>%
  summary(infer = T, )

hf6l.pairs <- na.omit(hf6l.pairs)
hf6l.pairs_dt <- data.table(hf6l.pairs)

group1 <- hf6l.pairs_dt$contrast

group1 <- gsub(" - protein0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - ph0", "", as.character(group1))

group1 <- as.factor(group1)

group2 <- hf6l.pairs_dt$contrast
group2 <- gsub("protein1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("ph1 - ", "", as.character(group2))

group2 <- as.factor(group2)
group2

hf6l.pairs_dt$group1 <- group1
hf6l.pairs_dt$group2 <- group2
hf6l.pairs_dt


hf6l.pairs_dt[, p_rounded := p_round(p.value,
                                      digits = 2)]
hf6l.pairs_dt[, p_pretty := p_format(p_rounded,
                                      digits = 2,
                                      accuracy = 1e-04,
                                      add.p = TRUE)]
hf6l.pairs_dt

hf6l.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                       "group 4", "group 4", "group 3", "group 3")

hf6l.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                           "group 5", "group 5", "group 3", "group 3",
                           "group 4", "group 4", "group 3", "group 3")


hf6l.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                           "group 6", "group 6", "group 4", "group 4",
                           "group 6", "group 6", "group 5", "group 5")

ggbarplot(data = hf6l.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf6l.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf6l.pairs_dt[-c(1:4),],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 30, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 6 Lean Mass", y ="Weight")

### C 12 Lean
anova(hf12l.lm)

hf12l.em <- emmeans(hf12l.lm, specs = c("sex", "ph", "protein"))
hf12l.em

hf12l.em_dt <- data.table(summary(hf12l.em))
hf12l.em_dt

hf12l.pairs <- contrast(hf12l.em,
                         method = "revpairwise",
                         simple = "each",
                         combine = T,
                         adjust = "sidak") %>%
  summary(infer = T, )

hf12l.pairs <- na.omit(hf12l.pairs)
hf12l.pairs_dt <- data.table(hf12l.pairs)

group1 <- hf12l.pairs_dt$contrast

group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))

group1 <- as.factor(group1)

group2 <- hf12l.pairs_dt$contrast
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))

group2 <- as.factor(group2)
group2

hf12l.pairs_dt$group1 <- group1
hf12l.pairs_dt$group2 <- group2
hf12l.pairs_dt


hf12l.pairs_dt[, p_rounded := p_round(p.value,
                                       digits = 2)]
hf12l.pairs_dt[, p_pretty := p_format(p_rounded,
                                       digits = 2,
                                       accuracy = 1e-04,
                                       add.p = TRUE)]
hf12l.pairs_dt

hf12l.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                        "group 4", "group 4", "group 3", "group 3")

hf12l.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                            "group 5", "group 5", "group 3", "group 3",
                            "group 4", "group 4", "group 3", "group 3")


hf12l.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                            "group 6", "group 6", "group 4", "group 4",
                            "group 6", "group 6", "group 5", "group 5")

ggbarplot(data = hf12l.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf12l.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf12l.pairs_dt[-c(1:4),],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 30, to = 65, by = 5)),
                     tip.length = 0.01) +
  labs(title = "Month 12 Lean Mass", y ="Weight")

#### 18 Lean Mass
anova(hf18l.lm)

hf18l.em <- emmeans(hf18l.lm, specs = c("sex", "ph", "protein"))
hf18l.em

hf18l.em_dt <- data.table(summary(hf18l.em))
hf18l.em_dt

hf18l.pairs <- contrast(hf18l.em,
                         method = "revpairwise",
                         simple = "each",
                         combine = T,
                         adjust = "sidak") %>%
  summary(infer = T, )

hf18l.pairs <- na.omit(hf18l.pairs)
hf18l.pairs_dt <- data.table(hf18l.pairs)

group1 <- hf18l.pairs_dt$contrast

group1 <- gsub(" - ph0", "", as.character(group1))
group1 <- gsub(" - sex0", "", as.character(group1))
group1 <- gsub(" - protein0", "", as.character(group1))

group1 <- as.factor(group1)

group2 <- hf18l.pairs_dt$contrast
group2 <- gsub("ph1 - ", "", as.character(group2))
group2 <- gsub("sex1 - ", "", as.character(group2))
group2 <- gsub("protein1 - ", "", as.character(group2))

group2 <- as.factor(group2)
group2

hf18l.pairs_dt$group1 <- group1
hf18l.pairs_dt$group2 <- group2
hf18l.pairs_dt


hf18l.pairs_dt[, p_rounded := p_round(p.value,
                                       digits = 2)]
hf18l.pairs_dt[, p_pretty := p_format(p_rounded,
                                       digits = 2,
                                       accuracy = 1e-04,
                                       add.p = TRUE)]
hf18l.pairs_dt


hf18l.em_dt$group <- c("group 6", "group 6", "group 5", "group 5", 
                        "group 4", "group 4", "group 3", "group 3")

hf18l.pairs_dt$group1 <- c("group 6", "group 5", "group 4", "group 3",
                            "group 5", "group 5", "group 3", "group 3",
                            "group 4", "group 4", "group 3", "group 3")


hf18l.pairs_dt$group2 <- c("group 6", "group 5", "group 4", "group 3",
                            "group 6", "group 6", "group 4", "group 4",
                            "group 6", "group 6", "group 5", "group 5")

ggbarplot(data = hf18l.em_dt, 
          x = "group", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex") +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title.x = element_blank()) +
  geom_errorbar(data = hf18l.em_dt, aes(y = emmean, ymin = lower.CL, ymax = upper.CL), width = 0.1, color = "black") +
  stat_pvalue_manual(data = hf18l.pairs_dt[-c(1:4),],
                     label = "p_pretty",
                     xmin = "group1",
                     xmax = "group2",
                     y.position = c(seq(from = 30, to = 44, by = 2)),
                     tip.length = 0.01) +
  labs(title = "Month 18 Lean Mass", y ="Weight")

## ALL Lean ##
l.em_dt <- rbind(hf6l.em_dt, hf12l.em_dt, hf18l.em_dt)
l.em_dt$age <- rep(c(6, 12, 18), each = 8)

hf6l.pairs_dt$age <- 6
hf12l.pairs_dt$age <- 12
hf18l.pairs_dt$age <- 18
l.pairs <- rbind(hf6l.pairs_dt, hf12l.pairs_dt, hf18l.pairs_dt)
l.pairs

l.em_dt$xlab <- paste(l.em_dt$group," ",l.em_dt$age,"", sep = "")
l.pairs$xlab1 <- paste(l.pairs$group1, l.pairs$age)
l.pairs$xlab2 <- paste(l.pairs$group2, l.pairs$age)

l.pairs.sig <- l.pairs[which(l.pairs$p.value <= 0.05),]
l.em_dt <- l.em_dt[order(l.em_dt$group),]

nrow(l.pairs)

#All Lean Graph

l.em_dt$group <- as.factor(l.em_dt$group)
l.em_dt <- l.em_dt[order(l.em_dt$group, decreasing = T),]

ggbarplot(data = l.em_dt, 
          x = "xlab", 
          y = "emmean", 
          fill = "group") +
  facet_wrap("sex", labeller = labeller(sex = sex_label))+
  theme_pubr() +
  geom_errorbar(data = l.em_dt, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE), width = 0.1, color = "black") +
  scale_fill_manual(values = c("brown1", "deepskyblue3", "chartreuse4", "purple"),
                    breaks = c("group 6", "group 5", "group 4", "group 3"),
                    labels = c("HFC", "HFCN", "HFB", "HFBN"),
                    name = "Diets") +
  stat_pvalue_manual(data = l.pairs.sig[-c(1:4, 9:12,15:18),],
                     label = "p_pretty",
                     xmin = "xlab1",
                     xmax = "xlab2",
                     y.position = c(seq(from = 33, to = 54, by = 3)),
                     tip.length = 0.01) +
  scale_x_discrete(labels = rep(c(6, 12, 18), times = 4)) +
  labs(title = "Lean Mass Over Time", y ="Lean Mass (g)", x = "Collection Timepoint (Months)")



#Facet Wrap#
hf.em_dt$cat <- "Total Mass (g)"
l.em_dt$cat <- "Lean Mass (g)"
f.em_dt$cat <- "Fat Mass (g)"
all.em_dt <- rbind(hf.em_dt, l.em_dt, f.em_dt)
all.em_dt$xlab <- paste(all.em_dt$group,"",all.em_dt$age,"", sep = "")
all.em_dt$group <- as.factor(all.em_dt$group)

all.em_dt <- all.em_dt[order(all.em_dt$age),]

all.em_dt$xlab <- paste(all.em_dt$age, rep(c("T", "L", "F"), each = 8, times = 3))

all.em_dt$xlab <- as.factor(all.em_dt$xlab)

levels(all.em_dt$xlab)

is.numeric(all.em_dt$emmean)

glab <- c(`group 3` = "HFBN", `group 4` = "HFB", 
          `group 5` = "HFCN", `group 6` = "HFC")
#names(glab) <-levels(all.em_dt$group)
names(glab) <- as.factor(names(glab))

is.character(names(glab))
is.factor(all.em_dt$group)

all.em_dt$group <- as.factor(all.em_dt$group)

gcolor <- c("brown1", "deepskyblue3", "chartreuse4", "purple")

#### Plot Combining ####
## Females First ##

levels(all.em_dt$cat)
all.em_dt$cat <- as.factor(all.em_dt$cat)
all.em_dt$cat <- factor(all.em_dt$cat, levels = rev(levels(all.em_dt$cat)))
all.em_dt$cat

levels(all.em_dt$group)
females <- all.em_dt[which(all.em_dt$sex == 0)]
males <- all.em_dt[which(all.em_dt$sex == 1)]

sex_label <- c("Female", "Male")
names(sex_label) <- c(0, 1)

group_label <- c("HFC", "HFCN",
                 "HFB", "HFBN")
names(group_label) <- c("group 6", "group 5", "group 4", "group 3")

#FEMALES#
## Plot Figure ##
fplt <- ggplot(data = females, aes(x = factor(xlab, levels = c("6 T", "6 L", "6 F",
                                                         "12 T", "12 L", "12 F",
                                                         "18 T", "18 L", "18 F")), 
                             y = emmean)) +
  #facet_wrap(vars(group), labeller = labeller(group = group_label))+
  geom_col_pattern(aes(pattern = cat, pattern_angle = cat),
                   pattern_fill = "white",
                   colour = "black",
                   fill = rep(gcolor, each = 9),
                   pattern_color = "black", 
                   pattern_angle = 45,
                   pattern_spacing = 0.075, pattern_density = rep(c(.35, .35, .75), times = 12),
                   pattern_key_scale_factor = 0.3) +
  scale_pattern_discrete(choices = c("none", "stripe", "circle")) +
  scale_pattern_spacing_discrete(range = c(0.01, 0.1)) +
  facet_wrap(~factor(group, levels = c("group 6", "group 5", 
                                       "group 4", "group 3")),
             labeller = as_labeller(glab)) +
  scale_fill_manual(values = gcolor, name = "Diets", labels = c("HFC", "HFCN", 
                                                                "HFB", "HFBN")) +
  geom_errorbar(data = females, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = xlab), width = 0.1, color = "black") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 52)) +
  scale_x_discrete(labels = rep(c(6, 12, 18), each = 3)) +
  guides(pattern_angle = "none") +
  theme_pubr() +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  labs(title = "Female C3H/HeJ Echo MRI Analysis", y ="Mass (g)", x = "Age (Months)")

fplt

#Type legend

ftype <- ggplot(data = females, aes(x = factor(xlab, levels = c("6 T", "6 L", "6 F",
                                                         "12 T", "12 L", "12 F",
                                                         "18 T", "18 L", "18 F")), 
                             y = emmean)) +
  geom_col_pattern(aes(pattern = cat, pattern_angle = cat), 
                   colour = 'black',
                   fill = 'purple', #switch fill and pattern fill for stripes
                   pattern_fill = 'white',
                   pattern_color = "black", pattern_angle = 45,
                   pattern_spacing = 0.05, pattern_density = .75,
                   pattern_key_scale_factor = 0.3) +
  scale_pattern_discrete(choices = c("none", "stripe", "circle"), name = "Mass Type") +
  scale_pattern_spacing_discrete(range = c(0.01, 0.1)) +
    facet_wrap(~factor(group, levels = c("group 3", "group 4", 
                                         "group 5", "group 6")),
               labeller = as_labeller(glab)) +
  geom_errorbar(data = females, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = xlab), width = 0.1, color = "black") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 52)) +
  guides(pattern_angle = "none") +
  theme_pubr() +
  theme(legend.box = "vertical",
        legend.title = element_blank()) +
  labs(title = "HF Female C3H/HeJ Echo MRI Analysis", y ="Mass (g)", x = "Age (Months)")

ftype

types = get_plot_component(last_plot(), "guide-box-top")
#Diet can be removed from legs

(legs <- plot_grid(types, plot_spacer() + theme_void(),
                   nrow = 2, ncol = 1, rel_heights = c(1, 2)))

plot_grid(fplt, types, nrow = 2, rel_widths = c(.8, .2), rel_heights = c(.8, .1))

#MALES#
## Plot Figure ##
mplt <- ggplot(data = males, aes(x = factor(xlab, levels = c("6 T", "6 L", "6 F",
                                                               "12 T", "12 L", "12 F",
                                                               "18 T", "18 L", "18 F")), 
                                   y = emmean)) +
  #facet_wrap(vars(group), labeller = labeller(group = group_label))+
  geom_col_pattern(aes(pattern = cat, pattern_angle = cat),
                   pattern_fill = "white",
                   colour = "black",
                   fill = rep(gcolor, each = 9),
                   pattern_color = "black", 
                   pattern_angle = 45,
                   pattern_spacing = 0.075, pattern_density = rep(c(.35, .35, .75), times = 12),
                   pattern_key_scale_factor = 0.3) +
  scale_pattern_discrete(choices = c("none", "stripe", "circle")) +
  scale_pattern_spacing_discrete(range = c(0.01, 0.1)) +
  facet_wrap(~factor(group, levels = c("group 6", "group 5", 
                                      "group 4", "group 3")),
             labeller = as_labeller(glab)) +
  scale_fill_manual(values = gcolor, name = "Diets", labels = c("HFC", "HFCN", 
                                                                "HFB", "HFBN")) +
  geom_errorbar(data = males, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = xlab), width = 0.1, color = "black") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 52)) +
  scale_x_discrete(labels = rep(c(6, 12, 18), each = 3)) +
  guides(pattern_angle = "none") +
  theme_pubr() +
  theme(legend.position = "none",
        legend.title = element_blank()) +
  labs(title = "Male C3H/HeJ Echo MRI Analysis", y ="Mass (g)", x = "Age (Months)")


mplt

#Type legend

mtype <- ggplot(data = males, aes(x = factor(xlab, levels = c("6 T", "6 L", "6 F",
                                                              "12 T", "12 L", "12 F",
                                                              "18 T", "18 L", "18 F")), 
                                  y = emmean)) +
  geom_col_pattern(aes(pattern = cat, pattern_angle = cat), 
                   colour = 'black',
                   fill = 'purple', #switch fill and pattern fill for stripes
                   pattern_fill = 'white',
                   pattern_color = "black", pattern_angle = 45,
                   pattern_spacing = 0.05, pattern_density = .75,
                   pattern_key_scale_factor = 0.3) +
  scale_pattern_discrete(choices = c("none", "stripe", "circle"), name = "Mass Type") +
  scale_pattern_spacing_discrete(range = c(0.01, 0.1)) +
  facet_wrap(~factor(group, levels = c("group 3", "group 4", 
                                       "group 5", "group 6")),
             labeller = as_labeller(glab)) +
  geom_errorbar(data = males, aes(y = emmean, ymin = emmean - SE, ymax = emmean + SE, x = xlab), width = 0.1, color = "black") +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 52)) +
  scale_x_discrete(labels = rep(c(6, 12, 18), each = 3)) +
  guides(pattern_angle = "none") +
  theme_pubr() +
  theme(legend.box = "vertical",
        legend.title = element_blank()) +
  labs(title = "HF Female C3H/HeJ Echo MRI Analysis", y ="Mass (g)", x = "Age (Months)")

mtype

type = get_plot_component(last_plot(), "guide-box-top")
#Diet can be removed from legs

plot_grid(mplt, type, nrow = 2, rel_widths = c(.8, .2), rel_heights = c(.8, .1))


#### Edit Tables ####

# Tables #
write.csv(anova(hf6.lm), "./Data/hf.anova.6.csv")
write.csv(anova(hf12.lm), "./Data/hf.anova.12.csv")
write.csv(anova(hf18.lm), "./Data/hf.anova.18.csv")
write.csv(anova(hf6f.lm), "./Data/hf.anova.f6.csv")
write.csv(anova(hf12f.lm), "./Data/hf.anova.f12.csv")
write.csv(anova(hf18f.lm), "./Data/hf.anova.f18.csv")
write.csv(anova(hf6l.lm),"./Data/hf.anova.l6.csv")
write.csv(anova(hf12l.lm),"./Data/hf.anova.l12.csv")
write.csv(anova(hf18l.lm),"./Data/hf.anova.l18.csv")
write.csv(hf.pairs.sig, "./Data/hfn.total.sig.diff.csv")
write.csv(f.pairs.sig, "./Data/hfn.fat.sig.diff.csv")
write.csv(l.pairs.sig, "./Data/hfn.lean.sig.diff.csv")

