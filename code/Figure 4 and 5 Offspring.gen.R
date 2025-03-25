# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2024-12-13 10:20:23
# @Last Modified by: dbc
# @Last Modified time: 2024-12-13 15:31:24
# @Description: analysis of offspring generation

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
library(lmPerm)
library(patchwork)
library(purrr)
library(rcompanion)
library(readr)
library(readxl)
library(tidyverse)
library(skimr)

# set work directory
setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()


# Read data with improved column specification
data_02gen <- read_excel("2022herbivory_data241212.xlsx", sheet = "2nd.gen")

data_02gen <- data_02gen %>% mutate(across(fw1st_mg:sla, ~as.numeric(.)))

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
# Assumption tests
# ==========================================================================
# Homogeneity of variance tests
# Define the variables to perform Levene's Test on
variables <- c(
               "total.mass",
               "total.leaf.mass",
               "stolon.mass",
               "petiole.mass",
               "root.mass",
               "no.ramets",
               "no.leaves",
               "stolon.length",
               "total.leaf.area",
               "mla",
               "sla",
               "ratio")


# Normalize the specified column
# data_02gen <- data_02gen %>%
#  mutate(across(all_of(variables), 
#                  ~ ifelse(. > 0, sqrt(.), NA)))
#                  ~ ifelse(. > 0, log(.), NA)))

# data_02gen <- data_02gen %>%
#   group_by(gene, treatment1st) %>%  # Replace group_variable with your grouping variable
#   mutate(across(all_of(variables), 
#                  ~ (.-mean(., na.rm = TRUE)) / sd(., na.rm = TRUE))) %>%
#   ungroup()


# Function to perform Levene's Test for a given variable
perform_levene_test <- function(var, data) {
  # Construct the formula for Levene's Test
  formula <- as.formula(paste(var, "~ gene * treatment1st"))

  # Subset data to exclude NA values for the current variable
  subset_data <- data %>% filter(!is.na(.data[[var]]))

  # Perform Levene's Test
  test_result <- leveneTest(formula, data = subset_data)

  # Extract the p-value
  test_result.tidy <- tidy(test_result) %>% mutate(variable = var, .before = 1)

  return(test_result.tidy)

}

# Apply Levene's Test to all variables and store p-values
levene_results_df <- map_dfr(variables, ~ perform_levene_test(.x, data = data_02gen))

levene_results_df <- levene_results_df %>%
mutate(statistic = round(statistic, 2),
       p.value = round(p.value, 3))

# > levene_results_df
#                var statistic p.value df df.residual
# 1       total.mass      0.87  0.5724 11         134
# 2  total.leaf.mass      1.72  0.0755 11         134
# 3      stolon.mass      0.98  0.4652 11         134
# 4     petiole.mass      0.69  0.7470 11         134
# 5        root.mass      1.83  0.0550 11         134
# 6        no.ramets      0.44  0.9337 11         134
# 7        no.leaves      1.61  0.1027 11         134
# 8    stolon.length      0.76  0.6780 11         134
# 9  total.leaf.area      1.06  0.3989 11         134
# 10             mla      3.70  0.0001 11         134
# 11             sla      2.92  0.0017 11         134
# 12           ratio      1.18  0.3033 11         134

write_csv(levene_results_df, "./results/levene_results_df.csv")

# Normality Test
# Function to perform normality test for a given variable
perform_normality_test <- function(var, data) {

  # Subset data to exclude NA values for the current variable
  subset_data <- data %>% filter(!is.na(.data[[var]])) %>% pull(var)

    # Perform Shapiro-Wilk test
    test_result <- shapiro.test(subset_data)

    # Extract the p-value
    test_result.tidy <- tidy(test_result) %>% mutate(variable = var, .before = 1)

    return(test_result.tidy)
}

# Apply normality test to all variables
normality_results_df <- map_dfr(variables, ~ perform_normality_test(.x, data = data_02gen))

#
normality_results_df <- normality_results_df %>%
mutate(statistic = round(statistic, 2),
       p.value = round(p.value, 3))

# > normality_results_df
# # A tibble: 12 × 4
#    variable        statistic p.value method
#    <chr>               <dbl>   <dbl> <chr>
#  1 total.mass           0.81   0     Shapiro-Wilk normality test
#  2 total.leaf.mass      0.76   0     Shapiro-Wilk normality test
#  3 stolon.mass          0.74   0     Shapiro-Wilk normality test
#  4 petiole.mass         0.81   0     Shapiro-Wilk normality test
#  5 root.mass            0.55   0     Shapiro-Wilk normality test
#  6 no.ramets            0.96   0.001 Shapiro-Wilk normality test
#  7 no.leaves            0.9    0     Shapiro-Wilk normality test
#  8 stolon.length        0.93   0     Shapiro-Wilk normality test
#  9 total.leaf.area      0.84   0     Shapiro-Wilk normality test
# 10 mla                  0.92   0     Shapiro-Wilk normality test
# 11 sla                  0.72   0     Shapiro-Wilk normality test
# 12 ratio                0.44   0     Shapiro-Wilk normality test

# Save results
write_csv(normality_results_df, "./results/normality_results_df.csv")

## =========================================================================
# permutation ANOVA: models with covariates
# ==========================================================================
# Function
perform_aovp <- function(var, data) {

  # Construct the formula for Levene's Test
  formula <- as.formula(paste(var, "~ fw1st_g + gene + treatment1st + gene:treatment1st"))

  # Subset data to exclude NA values for the current variable
  subset_data <- data %>% filter(!is.na(.data[[var]]))

  # Permutation test
  set.seed(2024)
  perm_anova <- aovp(formula, data = subset_data, np = 999)

  # Extract the p-value
  anova_result <- tidy(perm_anova) %>% mutate(variable = var, .before = 1)

  return(anova_result)
}

# Apply normality test to all variables
aovp_results_df <- map_dfr(variables, ~ perform_aovp(.x, data = data_02gen))

aovp_results_df <- aovp_results_df %>%
rename(
    sum.sq = `R Sum Sq`,
    mean.sq = `R Mean Sq`,
    iterations = Iter,
    p.value = `Pr(Prob)`
  ) %>%
mutate(
       sum.sq = round(sum.sq, 3),
       mean.sq = round(mean.sq, 3),
       p.value = round(p.value, 4))

aovp_results_df

# Save results to CSV file
write_csv(aovp_results_df, "./results/offspring.performance.anova.csv")

# ==========================================================================
# Permutation ANOVA and pairwise comparison
# noting: pairwise comparison did not account for cov effects
# ==========================================================================
# Main analysis function
aovp_posthoc <- function(data, var, cov, factors, np = 999) {

  # Build formula
  formula <- as.formula(paste(var, "~", paste(cov, factors[1], factors[2], paste0(factors[1], ":", factors[2]), sep = "+")))

  # Subset data to exclude NA values for the current variable
  subset_data <- data %>% filter(!is.na(.data[[var]]))

  # Perform permutation ANOVA
  set.seed(2024)
  perm_anova <- aovp(formula, data = subset_data, np = np)
  anova_result <- tidy(perm_anova) %>% mutate(variable = var, .before = 1)

  # Post hoc comparison for each main effect
  posthoc_results <- map(factors, function(factor) {
    posthoc_formula <- as.formula(paste(var, "~", factor))
    pairwisePermutationTest(posthoc_formula, data = subset_data, method = "fdr")
  })

  names(posthoc_results) <- factors

  # If there is an interaction effect, perform simple main effect analysis
  if (length(factors) > 1) {
    interaction_term <- paste(factors, collapse = ":")

    if (interaction_term %in% anova_result$term) {
      # if (anova_result$'Pr(Prob)'[anova_result$term == interaction_term] < 0.05) {
      # Create interaction effect combinations
      subset_data$interaction <- interaction(subset_data[[factors[1]]], subset_data[[factors[2]]])

      # Post hoc comparison for interaction effect
      posthoc_interaction <- pairwisePermutationTest(as.formula(paste(var, "~ interaction")),
                                                       data = subset_data, method = "fdr")

      posthoc_results$interaction <- posthoc_interaction

      posthoc_results_df <- posthoc_results %>%
                            bind_rows(.id = "source") %>%
                            mutate(across(everything(), as.character)) %>%
                            mutate(variable = var, .before = 1)
      # }
    }
  }

  list(anova = anova_result, posthoc = posthoc_results_df)
}

# Perform analysis for each subset
anova_results.posthoc <- map(variables, ~ aovp_posthoc(data = data_02gen, var = .x, cov = "fw1st_g", factors = c("gene", "treatment1st"))) %>% set_names(variables)

# Extract the parameters from results
extract_data <- function(anova_results.posthoc) {
  result_list <- map(names(anova_results.posthoc), function(subset_name) {
    subset <- anova_results.posthoc[[subset_name]]

    # Post-hoc results
    anova_df <- subset$anova %>%
    mutate(result_type = "anova_df", .before = 1)

    posthoc_df <- subset$posthoc %>%
    mutate(result_type = "posthoc_df", .before = 1)

      # Combine all results
    all_results <- bind_rows(anova_df, posthoc_df)
    all_results

  }
  )

  bind_rows(result_list)
}

# Use function to extract data
anova_results.posthoc_df <- extract_data(anova_results.posthoc)

# print(anova_results.posthoc_df)

# Save results to CSV file
write_csv(anova_results.posthoc_df, "./results/offspring.performance.posthoc.csv")

# ==========================================================================
# plotting
# ==========================================================================
# Define variables for plotting with their y-axis labels and limits
plot_info <- tibble(
                    variable = c("total.mass",
                                 "total.leaf.mass",
                                 "petiole.mass",
                                 "stolon.mass",
                                 "root.mass",
                                 "ratio",
                                 "no.ramets",
                                 "stolon.length",
                                 "total.leaf.area",
                                 "no.leaves",
                                 "mla",
                                 "sla"),
                    y_label = c("Total mass (g)",
                                "Leaf mass (g)",
                                "Petiole mass (g)",
                                "Stem mass (g)",
                                "Root mass (g)",
                                "Root-to-shoot ratio",
                                "Number of ramets",
                                "Stem length (cm)",
                                "Total leaf area (cm²)",
                                "Number of leaves",
                                "Mean leaf area (cm²)",
                                "Specific leaf area (cm²/g)"),
                     y_limit = c(0.10,
                                0.02,
                                0.006,
                                0.05,
                                0.02,
                                0.8,
                                10,
                                20,
                                10,
                                8,
                                4,
                                8000)
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
  summarize(mean = mean(value), .groups = "drop")

# Define color mapping
color_values <- c("Herbivory" = "#C25759", "Control" = "#599CB4")

# Define base theme to avoid code repetition
base_theme <- theme_minimal(base_family = "serif") +
  theme(
    panel.grid = element_blank(),            # Remove grid lines
    panel.background = element_blank(),      # Remove background color
    axis.line = element_line(linewidth = 1, color = "black"), # Add axis line
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.ticks.x = element_line(linewidth = 1, color = "black"),
    axis.ticks.y = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none"                  # Remove Legend
  )

# Define universal plotting function
create_plot <- function(var, y_label, y_lim, data, mean_data) {
  mean_var <- mean_data %>% filter(variable == var)

  # Determine x-axis label position
  x_var_label <- if(var %in% c("petiole.mass", "ratio", "sla", "total.leaf.area")) {
    "Parental"
  } else {
    ""
  }

  #
  subset_data <- data %>% filter(!is.na(!!sym(var)))

  # # Extracting max values
  # max_val <- max(subset_data[[var]], na.rm = TRUE)
  # # max values * 1.1
  # adjusted_limit <- max(y_lim, max_val * 1.1)

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
    labs(x = x_var_label, y = y_label) +
    ylim(0, y_lim) +
    base_theme

  # Adjust theme based on row position
  if(var %in% c("petiole.mass", "ratio", "sla", "total.leaf.area")){
    p <- p +
      theme(axis.text.x = element_text(size = 14, color = "black"),
            axis.title.x = element_text(size = 14, color = "black"))
  } else {
    p <- p +
      theme(axis.text.x = element_blank(),
            axis.title.x = element_blank())
  }

  return(p)
}

# Generate all plots
plots <- plot_info %>%
  mutate(plot = pmap(list(variable, y_label, y_limit), ~ create_plot(..1, ..2, ..3, data_02gen, data_long))) %>%
  pull(plot)

# Split plots into figure_3.7 and figure_3.8
figure_3.7_vars <- c("total.mass", "total.leaf.mass", "petiole.mass", "stolon.mass", "root.mass", "ratio")
figure_3.8_vars <- c("no.ramets", "stolon.length", "total.leaf.area",
                     "no.leaves", "mla", "sla")

# Extract respective plots
figure_4_plots <- plots[which(plot_info$variable %in% figure_4_vars)]
figure_5_plots <- plots[which(plot_info$variable %in% figure_5_vars)]

# Combine plots
figure_4 <- wrap_plots(figure_4_plots, ncol = 2, byrow = FALSE) +
  plot_annotation(tag_levels = 'A')

figure_5 <- wrap_plots(figure_5_plots, ncol = 2, byrow = FALSE) +
  plot_annotation(tag_levels = 'A')

# Export figures
ggexport(figure_4, filename = "./results/figure4.png",
     width = 2200,
     height = 3000,
     pointsize = 12,
     res = 300)

ggexport(figure_5, filename = "./results/figure5.png",
     width = 2200,
     height = 3000,
     pointsize = 12,
     res = 300)
