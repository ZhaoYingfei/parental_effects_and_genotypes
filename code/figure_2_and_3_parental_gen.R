# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2024-12-13 10:20:23
# @Last Modified by: dbc
# @Last Modified time: 2024-12-13 15:31:22
# @Description: analysis of parental generation

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
data_01gen <- read_excel("2022herbivory_data250325.xlsx", sheet = "1st.gen")
str(data_01gen)

data_01gen <- data_01gen %>%
# Remove rows where 'no.leaves' contains NA values
filter(!is.na(no.leaves)) %>%
# Convert all columns from 'stem.length' to 'sla' to numeric type
mutate(across(stem.length:sla, ~as.numeric(.)))

skim(data_01gen)
str(data_01gen)

# Data preprocessing
data_01gen$treatment1st <- as.factor(data_01gen$treatment1st)
data_01gen$gene <- as.factor(data_01gen$gene)
data_01gen$rep <- as.factor(data_01gen$rep)
data_01gen$nest.id <- as.factor(data_01gen$nest.id)

levels(data_01gen$treatment1st)
levels(data_01gen$gene)

# ==========================================================================
# Assumption tests
# ==========================================================================
# Homogeneity of variance tests
# Define the variables to perform Levene's Test on
variables <- c(
               "total.mass",
               "total.leaf.mass",
               "stem.mass",
               "petiole.mass",
               "root.mass",
               "no.ramets",
               "no.leaves",
               "stem.length",
               "total.leaf.area",
               "mla",
               "sla",
               "ratio")

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
levene_results_df <- map_dfr(variables, ~ perform_levene_test(.x, data = data_01gen))

levene_results_df <- levene_results_df %>%
mutate(statistic = round(statistic, 2),
       p.value = round(p.value, 3))

# > levene_results_df
# # A tibble: 12 × 5
#    variable        statistic p.value    df df.residual
#    <chr>               <dbl>   <dbl> <int>       <int>
#  1 total.mass           1      0.452    11          83
#  2 total.leaf.mass      0.66   0.775    11          83
#  3 stem.mass          0.69   0.746    11          83
#  4 petiole.mass         1.57   0.123    11          83
#  5 root.mass            0.39   0.954    11          83
#  6 no.ramets            0.74   0.693    11          83
#  7 no.leaves            0.59   0.833    11          83
#  8 stem.length        0.57   0.846    11          83
#  9 total.leaf.area      1.08   0.387    11          83
# 10 mla                  4.58   0        11          83
# 11 sla                  1.12   0.354    11          83
# 12 ratio                1.48   0.155    11          83

write_csv(levene_results_df, "./results/levene_results_df.parental.csv")

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
normality_results_df <- map_dfr(variables, ~ perform_normality_test(.x, data = data_01gen))

normality_results_df <- normality_results_df %>%
mutate(statistic = round(statistic, 2),
       p.value = round(p.value, 3))

# > normality_results_df
# # A tibble: 12 × 4
#    variable        statistic p.value method
#    <chr>               <dbl>   <dbl> <chr>
#  1 total.mass           0.9    0     Shapiro-Wilk normality test
#  2 total.leaf.mass      0.92   0     Shapiro-Wilk normality test
#  3 stem.mass          0.91   0     Shapiro-Wilk normality test
#  4 petiole.mass         0.35   0     Shapiro-Wilk normality test
#  5 root.mass            0.95   0.001 Shapiro-Wilk normality test
#  6 no.ramets            0.97   0.043 Shapiro-Wilk normality test
#  7 no.leaves            0.98   0.111 Shapiro-Wilk normality test
#  8 stem.length        0.98   0.097 Shapiro-Wilk normality test
#  9 total.leaf.area      0.83   0     Shapiro-Wilk normality test
# 10 mla                  0.82   0     Shapiro-Wilk normality test
# 11 sla                  0.79   0     Shapiro-Wilk normality test
# 12 ratio                0.99   0.37  Shapiro-Wilk normality test

# Save results
write_csv(normality_results_df, "./results/normality_results_df.parental.csv")

## =========================================================================
# permutation ANOVA: models with covariates
# ==========================================================================
# Function
perform_aovp <- function(var, data) {

  # Construct the formula for Levene's Test
  formula <- as.formula(paste(var, "~ gene + treatment1st + gene:treatment1st"))

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
aovp_results_df <- map_dfr(variables, ~ perform_aovp(.x, data = data_01gen))

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
write_csv(aovp_results_df, "./results/parental.performance.anova.csv")

# ==========================================================================
# Permutation ANOVA and pairwise comparison
# noting: pairwise comparison did not account for cov effects
# ==========================================================================
# Main analysis function
aovp_posthoc <- function(data, var, factors, np = 999) {

  # Build formula
  formula <- as.formula(paste(var, "~", paste(factors[1], factors[2], paste0(factors[1], ":", factors[2]), sep = "+")))

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
anova_results.posthoc <- map(variables, ~ aovp_posthoc(data = data_01gen, var = .x, factors = c("gene", "treatment1st"))) %>% set_names(variables)

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
write_csv(anova_results.posthoc_df, "./results/parental.performance.posthoc.csv")

# ==========================================================================
# plotting
# ==========================================================================
# Define variables for plotting with their y-axis labels and limits
plot_info <- tibble(
                    variable = c("total.mass",
                                 "total.leaf.mass",
                                 "petiole.mass",
                                 "stem.mass",
                                 "root.mass",
                                 "ratio",
                                 "no.ramets",
                                 "stem.length",
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
                     y_limit = c(1,
                                0.3,
                                0.2,
                                0.5,
                                0.1,
                                0.3,
                                30,
                                100,
                                80,
                                25,
                                10,
                                600)
                                )

# sapply(plot_info$variable, function(var) {
#   range(data_01gen[[var]], na.rm = TRUE)
# })

# Convert data to long format and calculate means
data_long <- data_01gen %>%
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
    labs(x = "", y = y_label) +
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
  mutate(plot = pmap(list(variable, y_label, y_limit), ~ create_plot(..1, ..2, ..3, data_01gen, data_long))) %>%
  pull(plot)

# Split plots into figure_3.2 and figure_3.3
figure_2_vars <- c("total.mass", "total.leaf.mass", "petiole.mass", "stem.mass", "root.mass", "ratio")
figure_3_vars <- c("no.ramets", "stem.length", "total.leaf.area",
                     "no.leaves", "mla", "sla")

# Extract respective plots
figure_2_plots <- plots[which(plot_info$variable %in% figure_2_vars)]
figure_3_plots <- plots[which(plot_info$variable %in% figure_3_vars)]

# Combine plots
figure_2 <- wrap_plots(figure_2_plots, ncol = 2, byrow = FALSE) +
  plot_annotation(tag_levels = 'A')

figure_3 <- wrap_plots(figure_3_plots, ncol = 2, byrow = FALSE) +
  plot_annotation(tag_levels = 'A')

# Export figures

ggexport(figure_2, filename = "./results/figure2.png",
     width = 2200,
     height = 3000,
     pointsize = 12,
     res = 300)

ggexport(figure_3, filename = "./results/figure3.png",
     width = 2200,
     height = 3000,
     pointsize = 12,
     res = 300)

