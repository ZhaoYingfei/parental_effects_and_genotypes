# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2024-12-11 22:17:10
# @Last Modified by: dbc
# @Last Modified time: 2024-12-13 15:33:44
# @Description: leaf loss analysis of clonal offspring

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

#Calculating losses and loss rates
data_02gen <- data_02gen %>% mutate(loss.fw = 1- fw.leaf.t1/fw.leaf.t0)
data_02gen <- data_02gen %>% mutate(loss.la = 1- leafarea.t1/leafarea.t0)

range(data_02gen$loss.fw, na.rm = TRUE)
range(data_02gen$loss.la, na.rm = TRUE)

# Function to fit models and perform chi-square tests
analyze_loss <- function(response_var, data) {
  # Fit full model
  model_full <- glm(as.formula(paste(response_var, "~ fw1st_g + gene + treatment1st + treatment1st:gene")),
                    family = quasibinomial(),
                    data = data)

  # Create nested models
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
  var_name <- if(response_var == "loss.fw") "loss of fresh mass" else "loss of leaf area"

  results <- data.frame(
    var = rep(var_name, 4),
    effect = c("G x H", "Herbivory (H)", "Gene (G)", "initial fresh weight"),
    df = sapply(chisq_tests, function(x) tidy(x)$df[2]),
    chi_square = sapply(chisq_tests, function(x) tidy(x)$deviance[2]),
    p_value = sapply(chisq_tests, function(x) tidy(x)$p.value[2])
  ) %>%
    mutate(
      chi_square = format(chi_square, nsmall = 2, scientific = FALSE),
      p_value = format(p_value, nsmall = 4, scientific = FALSE)
    )

  return(results)
}

# Analyze both response variables
chisq_results.fw <- analyze_loss("loss.fw", data_02gen)
chisq_results.la <- analyze_loss("loss.la", data_02gen)

# Combine results
chisq_results.loss <- bind_rows(chisq_results.fw, chisq_results.la)

# Export results
write_csv(chisq_results.loss, "./results/chisq_results.loss.csv")

# ==========================================================================
# plotting loss data
# ==========================================================================
# Define variables for plotting with their y-axis labels and limits
plot_info <- tibble(
  variable = c("loss.fw",
               "loss.la"),
  y_label = c("Loss of fresh mass (%)",
              "Loss of leaf area (%)"),
  y_limit = c(70, 
              20)
)

# sapply(plot_info$variable, function(var) {
#   range(data_02gen[[var]], na.rm = TRUE)
# })

# Convert data to long format and calculate means
data_long <- data_02gen %>%
  pivot_longer(
    cols = all_of(plot_info$variable),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  group_by(gene, treatment1st, variable) %>%
  summarize(mean = mean(value) * 100, .groups = "drop") # Convert to Percent

# Define color mapping
color_values <- c("Herbivory" = "#C25759", "Control" = "#599CB4")

# Define base theme to avoid code repetition
base_theme <- theme_minimal(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.ticks.x = element_line(linewidth = 1, color = "black"),
    axis.ticks.y = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none"
  )

# Define universal plotting function
create_plot <- function(var, y_label, y_lim, data, mean_data) {
  mean_var <- mean_data %>% filter(variable == var)

  subset_data <- data %>%
    filter(!is.na(!!sym(var))) %>%
    mutate(!!sym(var) := !!sym(var) * 100) # Convert to percentage

  p <- ggplot(subset_data, aes(x = treatment1st, y = !!sym(var), fill = treatment1st)) +
    geom_boxplot(aes(colour = treatment1st), size = 1.2, width = 0.5, outlier.shape = NA, alpha = 0.5) +
    geom_line(data = mean_var, aes(x = treatment1st, y = mean, group = gene),
              color = "grey50", size = 1.2) +
    geom_point(data = mean_var, aes(x = treatment1st, y = mean),
               color = "grey50", size = 2) +
    geom_text_repel(data = mean_var,
                    aes(x = treatment1st, y = mean, label = as.character(gene)),
                    color = "black", size = 4, max.overlaps = Inf) +
    scale_colour_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    labs(x = "Parental", y = y_label) +
    ylim(y_lim, 100) +
    base_theme

    p <- p +
      theme(axis.text.x = element_text(size = 14, color = "black"),
            axis.title.x = element_text(size = 14, color = "black"))

  return(p)
}

# Generate all plots
plots <- plot_info %>%
  mutate(plot = pmap(list(variable, y_label, y_limit),
                     ~ create_plot(..1, ..2, ..3, data_02gen, data_long))) %>%
  pull(plot)

# Combine plots
figure_loss <- wrap_plots(plots, ncol = 2, byrow = FALSE) +
  plot_annotation(tag_levels = 'A')

figure_loss

# Export figures
ggexport(figure_loss, filename = "./results/figure_loss.png",
         width = 2200,
         height = 1000,
         pointsize = 12,
         res = 300)