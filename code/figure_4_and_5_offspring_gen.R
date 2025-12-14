# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2025-12-08
# @Last Modified by: dbc
# @Last Modified time: 2025-12-08
# @Description: analysis of offspring generation

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
library(gridExtra)
library(broom.mixed)
library(multcompView)
library(multcomp)

# set work directory
setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()


# Read data
data_02gen <- read_excel("./data/2022herbivory_data250325.xlsx", sheet = "2nd.gen")

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

# data_02gen$fw1st_g <- data_02gen$fw1st_mg / 1000

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
  
  fml <- as.formula(paste(var, "~ gene * treatment1st + (1 | nest.id)"))
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
  Variable = c("total.mass","total.leaf.mass","petiole.mass","stem.mass",
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
write_csv(final_anova_table, "./results/offspring.performance.anova.csv")

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
png("./results/offspring.plots_qq.png", width = 1600, height = 1200)
grid.arrange(grobs = plots_qq, ncol = 4)
dev.off()

# Arrange plots_resid_fit in a 2-column grid and save to a PNG file
png("./results/offspring.plots_resid_fit.png", width = 1600, height = 1200)
grid.arrange(grobs = plots_resid_fit, ncol = 4)
dev.off()  

# Performing post-hoc tests (Tukey-adjusted multiple comparisons) for main and simple effects
# ==========================================================================
# Post-hoc Tests (Multiple Comparisons)
# ==========================================================================

# Define a function to calculate EMMs and generate Letter Groupings (CLD) and simple effect pairs
get_posthoc_results <- function(model, var_name) {
  
  # ----------------------------------------------------------------------
  # Request 1: Differences between Genotypes (Main Effect) - CLD
  # ----------------------------------------------------------------------
  emm_gene <- emmeans(model, ~ gene)
  cld_gene <- cld(emm_gene, 
                  alpha = 0.05, 
                  Letters = letters, 
                  adjust = "none", 
                  reversed = TRUE) 
  
  # ----------------------------------------------------------------------
  # Request 2: Differences between Treatments WITHIN Genotypes
  # ----------------------------------------------------------------------
  emm_interaction_pairs <- pairs(emmeans(model, ~ treatment1st | gene), 
                                 adjust = "none")
  
  results_df <- as.data.frame(summary(emm_interaction_pairs))

  results_df$stars <- p_to_stars(results_df$p.value)
  
  results_df <- results_df %>% 
    dplyr::select(gene, contrast, estimate, SE, t.ratio, df, p.value, stars)
  
  # Return a list containing both results
  return(list(
    Variable = var_name,
    Genotype_Main_Effect_CLD = cld_gene,
    Treatment_Within_Genotype_Pairs = results_df
  ))
}

posthoc_outputs <- purrr::map2(lmm_fits, variables, get_posthoc_results)

# Extract CLD results
results_main_effect_cld <- purrr::map_dfr(posthoc_outputs, function(item) {
  df <- as.data.frame(item$Genotype_Main_Effect_CLD)
  df$Variable <- item$Variable 
  df$.group <- trimws(df$.group)
  return(df %>% dplyr::select(Variable, everything()))
})

# Extract P-value and star results
results_simple_effect_stars <- purrr::map_dfr(posthoc_outputs, function(item) {
  df <- as.data.frame(item$Treatment_Within_Genotype_Pairs)
  df$Variable <- item$Variable 
  return(df %>% dplyr::select(Variable, everything()))
})

write_csv(results_main_effect_cld, "./results/offspring.posthoc.gene.csv")
write_csv(results_simple_effect_stars, "./results/offspring.posthoc.interaction.csv")

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
                     y_limit = c(0.05,
                                0.02,
                                0.006,
                                0.05,
                                0.025,
                                1,
                                10,
                                20,
                                10,
                                5,
                                4,
                                10000)
                                )

sig_group <- tibble::tribble(
  ~variable,     ~gene,   ~label,
# "total.mass","G6","a",
# "total.mass","G2","a",
# "total.mass","G1","a",
# "total.mass","G4","a",
# "total.mass","G3","a",
# "total.mass","G5","a",
# "total.leaf.mass","G1","a",
# "total.leaf.mass","G6","ab",
# "total.leaf.mass","G5","ab",
# "total.leaf.mass","G2","ab",
# "total.leaf.mass","G3","b",
# "total.leaf.mass","G4","b",
# "stem.mass","G4","a",
# "stem.mass","G6","a",
# "stem.mass","G3","a",
# "stem.mass","G2","a",
# "stem.mass","G1","a",
# "stem.mass","G5","a",
# "petiole.mass","G4","a",
# "petiole.mass","G1","a",
# "petiole.mass","G6","ab",
# "petiole.mass","G5","ab",
# "petiole.mass","G2","ab",
# "petiole.mass","G3","b",
# "root.mass","G2","a",
# "root.mass","G6","ab",
# "root.mass","G1","ab",
# "root.mass","G3","ab",
# "root.mass","G5","b",
# "root.mass","G4","b",
# "no.ramets","G1","a",
# "no.ramets","G6","a",
# "no.ramets","G2","a",
# "no.ramets","G5","a",
# "no.ramets","G3","a",
# "no.ramets","G4","a",
"no.leaves","G6","a",
"no.leaves","G1","a",
"no.leaves","G2","ab",
"no.leaves","G5","ab",
"no.leaves","G3","bc",
"no.leaves","G4","c",
# "stem.length","G6","a",
# "stem.length","G1","ab",
# "stem.length","G2","ab",
# "stem.length","G3","ab",
# "stem.length","G5","b",
# "stem.length","G4","ab",
# "total.leaf.area","G1","a",
# "total.leaf.area","G6","a",
# "total.leaf.area","G4","ab",
# "total.leaf.area","G2","ab",
# "total.leaf.area","G5","ab",
# "total.leaf.area","G3","b",
"mla","G4","a",
"mla","G1","b",
"mla","G6","bc",
"mla","G5","bc",
"mla","G2","bc",
"mla","G3","c",
"sla","G4","a",
"sla","G6","b",
"sla","G2","b",
"sla","G3","b",
"sla","G1","b",
"sla","G5","b",
"ratio","G2","a",
"ratio","G3","b",
"ratio","G1","b",
"ratio","G6","b",
"ratio","G5","b",
"ratio","G4","b",
)

sig_inter <- tibble::tribble(
  ~variable,     ~gene,   ~label,
"root.mass","G2","**",
"sla","G4","***",
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
  summarize(mean = mean(value),
            se = sd(value) / sqrt(n()),
            .groups = "drop")

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
    legend.position = "none"               
  )

# Define universal plotting function
create_grouped_barplot <- function(var, y_label, y_lim, data) {
    subset_data <- data %>% filter(variable == var)
  
    sig_data1 <- sig_group %>%
    filter(variable == var) %>%
    left_join(subset_data %>% group_by(gene) %>% summarize(y = y_lim), by = "gene")
    sig_data2 <- sig_inter %>%
    filter(variable == var) %>%
    left_join(subset_data %>% group_by(gene) %>% summarize(y = max(mean) + 0.05), by = "gene")

    ggplot(subset_data, aes(x = gene, y = mean, fill = treatment1st)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
      geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                        position = position_dodge(width = 0.8),
                        width = 0.2) +
      scale_fill_manual(values = color_values) +
      labs(x = NULL, y = y_label) +
      coord_cartesian(ylim = c(0, y_lim)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
      geom_text(data = sig_data1, aes(x = gene, y = y , label = label), 
              inherit.aes = FALSE, size = 5, family = "serif") +
      geom_text(data = sig_data2, aes(x = gene, y = y , label = label), 
              inherit.aes = FALSE, size = 5, family = "serif") +
      base_theme +
      theme(axis.text.x = element_text(size = 14, color = "black"),
            axis.title.x = element_blank())
  }

bar_plots <- plot_info %>%
    mutate(plot = pmap(list(variable, y_label, y_limit), ~ create_grouped_barplot(..1, ..2, ..3, data_long))) %>%
    pull(plot)

# Extract Legend
legend_plot <- ggplot(data_long %>% filter(variable == "total.mass"),
                      aes(x = gene, y = mean, fill = treatment1st)) +
  geom_col(position = position_dodge(0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = color_values) +
  theme_minimal(base_family = "serif") +
  labs(fill = "Parental") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, family = "serif"),
    legend.text = element_text(size = 12, family = "serif"),
    legend.key.size = unit(0.5, "cm")
  )+
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5)) 

legend_grob <- cowplot::get_legend(legend_plot)

# Split plots into figure_4 and figure_5
figure_4_vars <- c("total.mass", "total.leaf.mass", "petiole.mass", "stem.mass", "root.mass", "ratio")
figure_5_vars <- c("no.ramets", "stem.length", "total.leaf.area", "no.leaves", "mla", "sla")

bold_tag_theme <- theme(plot.tag = element_text(face = 'bold'))

# Generate main figure
main_plot_4 <- wrap_plots(bar_plots[plot_info$variable %in% figure_4_vars], ncol = 2) +
    plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')& bold_tag_theme
main_plot_5 <- wrap_plots(bar_plots[plot_info$variable %in% figure_5_vars], ncol = 2)+
    plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')& bold_tag_theme

# Combine pictures
figure_4 <- ggdraw() +
  draw_plot(main_plot_4, x = 0, y = 0, width = 1, height = 1) +
  draw_grob(legend_grob, x = 0.80, y = 0.75, width = 0.22, height = 0.3)

figure_5 <- ggdraw() +
  draw_plot(main_plot_5, x = 0, y = 0, width = 1, height = 1) +
  draw_grob(legend_grob, x = 0.80, y = 0.75, width = 0.22, height = 0.3)

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
