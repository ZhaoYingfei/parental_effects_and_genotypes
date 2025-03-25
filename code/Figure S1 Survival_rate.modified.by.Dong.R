# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2024-12-11 22:17:10
# @Last Modified by: dbc
# @Last Modified time: 2024-12-13 15:30:59
# @Description: survival analysis of clonal offspring

# Cleaning memory
cat("\014")
rm(list = ls())
gc()

# Load necessary packages
library(broom)
library(car)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)
library(readr)
library(readxl)
library(tidyverse)
library(skimr)

# set work directory
setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()


# Read data with improved column specification
data_02gen <- read_excel("2022herbivory_data241212.xlsx", sheet = "2nd.gen")

data_02gen <- data_02gen %>% mutate(across(fw1st_mg:sla, ~ as.numeric(.)))

skim(data_02gen)
str(data_02gen)

# Data preprocessing
data_02gen$treatment1st <- as.factor(data_02gen$treatment1st)
data_02gen$gene <- as.factor(data_02gen$gene)
data_02gen$rep <- as.factor(data_02gen$rep)
data_02gen$nest.id <- as.factor(data_02gen$nest.id)

levels(data_02gen$treatment1st)
levels(data_02gen$gene)

data_02gen$fw1st_g <- data_02gen$fw1st_mg / 1000

# ==========================================================================
# analysis of survival rate
# ==========================================================================
# Data preprocessing
table(data_02gen$survival)

# Fit full model
model_full <- glm(survival ~ fw1st_g + gene + treatment1st + treatment1st:gene,
                  family = binomial("cloglog"),
                  data = data_02gen)

# Remove terms sequentially using update
models <- list(
               full = model_full,
               m3   = update(model_full, . ~ . - treatment1st:gene),
               m2   = update(model_full, . ~ . - treatment1st:gene - treatment1st),
               m1   = update(model_full, . ~ . - treatment1st:gene - treatment1st - gene),
               null = update(model_full, . ~ . - treatment1st:gene - treatment1st - gene - fw1st_g)
               )

# Perform chi-square tests
chisq_tests <- list(
                    gxh = anova(models$m3, models$full, test = "Chisq"),
                    treatment = anova(models$m2, models$m3, test = "Chisq"),
                    gene = anova(models$m1, models$m2, test = "Chisq"),
                    fw1st_mg = anova(models$null, models$m1, test = "Chisq")
                    )


# Create results data frame
chisq_results <- data.frame(
                            effect = c("G x H", "Herbivory (H)", "Gene (G)", "initial fresh weight"),
                            df = sapply(chisq_tests, function(x) tidy(x)$df[2]),
                            chi_square = sapply(chisq_tests, function(x) tidy(x)$deviance[2]),
                            p_value = sapply(chisq_tests, function(x) tidy(x)$p.value[2])
                            )


# Format numeric values: fixed decimal places for chi-square and p-values
chisq_results <- chisq_results %>%
mutate(chi_square = format(chi_square, nsmall = 2, scientific = FALSE),
       p_value = format( p_value, nsmall = 4, scientific = FALSE))

chisq_results

# Export chi-square test results to CSV file
write_csv(chisq_results, "./results/chisq_results.survival.csv")

# ==========================================================================
# plotting
# ==========================================================================
# Calculate survival rates
survival_rate <- data_02gen %>%
group_by(treatment1st) %>%
summarise(survival_count = sum(survival == 1, na.rm = TRUE),
          total_count = n(),
          survival_rate = survival_count / total_count * 100,
          .groups = "drop")

print(survival_rate)

survival_line <- data_02gen %>%
group_by(treatment1st, gene) %>%
summarise(survival_count = sum(survival == 1, na.rm = TRUE),
          total_count = n(),
          survival_rate = survival_count / total_count * 100,
          .groups = "drop")

print(survival_line)

# Create plot with improved aesthetics
plot_sr <- ggplot(survival_rate, aes(x = treatment1st, y = survival_rate, fill = treatment1st)) +
geom_bar(stat = "identity", width = 0.5) +
geom_line(data = survival_line, aes(x = treatment1st, y = survival_rate, group = gene),
          color = "grey50", size = 1.2) +
geom_point(data = survival_line, aes(x = treatment1st, y = survival_rate, group = gene),
           color = "grey50", size = 2) +
geom_text_repel(data = survival_line,
                aes(x = treatment1st,
                    y = survival_rate,
                    label = as.character(gene)),
                color = "black",
                size = 4,
                segment.size = NA) +
labs(x = "Parental", y = "Survival rate (%)") +
scale_fill_manual(values = c("Control" = "#599CB4", "Herbivory" = "#C25759")) +
scale_y_continuous(expand = c(0, 0), limits = c(0, 110), breaks = seq(0, 100, by = 20)) +
theme_minimal(base_family = "serif") +
theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(size = 1, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      axis.ticks = element_line(size = 1, color = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      legend.position = "none"
      )
#
print(plot_sr)

# Save plot
ggexport(plot_sr, filename = "./results/figure3.6survival.png",
         width = 1000,
         height = 1000,
         pointsize = 12,
         res = 300)
