# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2025-12-08
# @Last Modified by: dbc
# @Last Modified time: 2025-12-08
# @Description: analysis of offspring generation whthin initial weight

# Cleaning memory
cat("\014")
rm(list = ls())
gc()

# Load necessary packages
library(broom)
library(ggpubr)
library(patchwork)
library(readxl)
library(tidyverse)
library(skimr)
library(cowplot)
library(lme4)
library(lmerTest)
library(emmeans)

# set work directory
setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()


# Read data
data_02gen <- read_excel("./data/2022herbivory_data250325.xlsx", sheet = "2nd.gen")

data_02gen <- data_02gen %>% mutate(across(fw1st_mg:sla, ~as.numeric(.)))

skim(data_02gen)
str(data_02gen)

# Data preprocessing
data_02gen$treatment1st <- as.factosr(data_02gen$treatment1st)
data_02gen$gene <- as.factor(data_02gen$gene)
data_02gen$rep <- as.factor(data_02gen$rep)
data_02gen$nest.id <- as.factor(data_02gen$nest.id)

levels(data_02gen$treatment1st)
levels(data_02gen$gene)

data_02gen$fw1st_g <- data_02gen$fw1st_mg / 1000

# Define the variables
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

# Performing Mixed-Effects Model (LMM) analysis using `lmer` with nested random effect
## =========================================================================
# Permutation ANOVA with nested error term
# ==========================================================================
# Function to fit LMM for each dependent variable
fit_lmm <- function(var, data) {
  # Ensure nest.id is a factor type
  if (!is.factor(data$nest.id)) {
    data$nest.id <- as.factor(data$nest.id)
  }
  
  fml <- as.formula(paste(var, "~ fw1st_g + gene * treatment1st + (1 | nest.id)"))
  subset_data <- data[!is.na(data[[var]]), , drop = FALSE]
  
  set.seed(2024)
  fit <- lmer(fml, data = subset_data)
  
  return(fit)
}

lmm_fits <- setNames(
  purrr::map(variables, ~ fit_lmm(.x, data_02gen)),
  variables
)

anova_results <- purrr::map(lmm_fits, ~ anova(.x, type = "III"))
anova_results

# Defining function to convert p-values to significance stars
# ==========================================================================
# Function to convert P-values to stars
# ==========================================================================
# 0.05: *; 0.01: **; 0.001: ***
p_to_stars <- function(p) {
  sapply(p, function(x) {
    if (x <= 0.001) {
      return("***")
    } else if (x <= 0.01) {
      return("**")
    } else if (x <= 0.05) {
      return("*")
    } else {
      return("ns") # non-significant
    }
  })
}

# Defining lookup table for variable names to descriptive trait names
# variable → Traits 
traits_lut <- tibble(
  variable = c("total.mass","total.leaf.mass","petiole.mass","stem.mass",
               "root.mass","ratio","no.ramets","stem.length","total.leaf.area",
               "no.leaves","mla","sla"),
  Traits   = c("Total mass","Leaf mass","Petiole mass","Stem mass","Root mass",
               "Root-to-shoot ratio","Number of ramets","Stem length",
               "Total leaf area","Number of leaves","Mean leaf area",
               "Specific leaf area")
)

format_anova_table <- function(anova_df, var_name) {
  
  df_clean <- anova_df %>%
    tibble::rownames_to_column(var = "Effect") %>%
    mutate(Variable = var_name)
  
  df_clean <- df_clean %>%
    dplyr::select(
      Variable,
      Effect,
      NumDF,
      DenDF,
      `Sum Sq`,
      `Mean Sq`,
      `F value`,
      P_Value = starts_with("Pr") 
    ) %>%
    
    mutate(Stars = p_to_stars(P_Value)) %>%
    
    mutate(
      `Sum Sq` = round(`Sum Sq`, 3),
      `Mean Sq` = round(`Mean Sq`, 3),
      `F value` = round(`F value`, 3),
      P_Value = round(P_Value, 3) 
    ) %>%
    
    mutate(
      Effect = case_when(
        Effect == "gene" ~ "Genotype",
        Effect == "treatment1st" ~ "Treatment (1st Gen)",
        Effect == "gene:treatment1st" ~ "Genotype x Treatment",
        TRUE ~ Effect
      )
    )
  
  return(df_clean)
}

final_anova_table01 <- purrr::map2_dfr(
  .x = anova_results, 
  .y = names(anova_results), 
  .f = format_anova_table
)

final_anova_table <- final_anova_table01 %>%
  left_join(traits_lut, by = "Variable") %>%
  dplyr::select(-Variable) %>%
  dplyr::select(Traits, everything())

# Exporting the final ANOVA table to a CSV file
write_csv(final_anova_table, "./results/offspring.performance.anova.withinfw.csv")

# Extracting summaries of fixed effects and variance components from LMMs
fixed_effects_summary <- purrr::map(lmm_fits, ~ summary(.x)$coefficients)
fixed_effects_summary

random_effects_variance <- purrr::map(lmm_fits, function(fit) {
  s <- summary(fit)
  var_comp <- data.frame(s$varcor)
  result <- var_comp %>%
    dplyr::select(grp, vcov) %>%
    dplyr::rename(Group = grp, Variance = vcov)
  return(result)
})

# Define the residual diagnostics function
perform_diagnostics <- function(model_fit, var_name) {
  # Extract the residuals and fitted values from the model
  # .resid is the raw residuals
  # .fitted is the predicted values from the model
  model_data <- augment(model_fit)
  
  # --- Diagnostic Plot 1: Normality of Residuals (QQ Plot) ---
  qq_plot <- ggplot(model_data, aes(sample = .resid)) +
    stat_qq() +
    stat_qq_line() +
    labs(
      title = paste0("Normality Test: QQ Plot (", var_name, ")"),
      x = "Theoretical Quantiles",
      y = "Standardized Residuals"
    ) +
    theme_bw()

  # --- Diagnostic Plot 2: Homoscedasticity (Residuals vs. Fitted) ---
  resid_fit_plot <- ggplot(model_data, aes(x = .fitted, y = .resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    geom_smooth(method = "loess", se = FALSE, color = "blue") + # Add smoothing curve
    labs(
      title = paste0("Homoscedasticity Test: Residuals vs. Fitted (", var_name, ")"),
      x = "Fitted Values",
      y = "Residuals"
    ) +
    theme_bw()

  return(list(qq_plot = qq_plot, resid_fit_plot = resid_fit_plot))
}

# Executing diagnostics and storing the results (plots)
diagnostic_plots <- purrr::map2(lmm_fits, variables, perform_diagnostics)
library(gridExtra)

# Displaying all QQ plots and Residuals vs. Fitted plots in separate windows/grids
plots_qq <- purrr::map(diagnostic_plots, ~ .x$qq_plot)
plots_resid_fit <- purrr::map(diagnostic_plots, ~ .x$resid_fit_plot)

# Arrange plots_qq in a 2-column grid and save to a PNG file
png("./results/offspring.withinfw.plots_qq.png", width = 1600, height = 1200)
grid.arrange(grobs = plots_qq, ncol = 4)
dev.off()

# Arrange plots_resid_fit in a 2-column grid and save to a PNG file
png("./results/offspring.withinfw.plots_resid_fit.png", width = 1600, height = 1200)
grid.arrange(grobs = plots_resid_fit, ncol = 4)
dev.off()  