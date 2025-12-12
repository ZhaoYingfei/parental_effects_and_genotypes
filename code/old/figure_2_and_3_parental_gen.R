# -*- coding: utf-8 -*-
# @Author: Zhao Yingfei
# @Date: 2025-10-01 
# @Last Modified by: dbc
# @Last Modified time: 2025-10-01
# @Description: analysis of parental generation

# Cleaning memory
cat("\014")
rm(list = ls())
gc()

# Load necessary packages
library(broom)
library(car)
library(ggpubr)
library(lmPerm)
library(patchwork)
library(rcompanion)
library(readxl)
library(tidyverse)
library(skimr)
library(cowplot)
library(multcompView)
library(permuco)

# set work directory
setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()


# # Read data
data_01gen <- read_excel("./data/2022herbivory_data250325.xlsx", sheet = "1st.gen")
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

# write_csv(levene_results_df, "./results/parental.levene.results.csv")

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
# write_csv(normality_results_df, "./results/parental.normality.result.csv")

## =========================================================================
# Permutation ANOVA with nested error term
# ==========================================================================
# Function
perform_aovp <- function(var, data) {
  fml <- as.formula(paste(var, "~ gene + treatment1st + gene:treatment1st + Error(nest.id/(gene))"))
  subset_data <- data[!is.na(data[[var]]), , drop = FALSE]
  set.seed(2024)
  fit <- permuco::aovperm(fml, data = subset_data, np = 999)
  summary(fit)
}

aovp_results <- setNames(
  purrr::map(variables, ~ perform_aovp(.x, data_01gen)),
  variables
)
aovp_results_df0 <- imap_dfr(aovp_results, ~{
  df <- as.data.frame(.x, check.names = FALSE)
  df$effect <- rownames(df)
  df$variable <- .y
  rownames(df) <- NULL
  df
})

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

# Rename and keep the specified column
aovp_results_df <- aovp_results_df0 %>%
  rename(
    parametric.p.value = `parametric P(>F)`,
    resampled.p.value  = `resampled P(>F)`
  )  %>%
  left_join(traits_lut, by = "variable") %>%
  select(Traits, effect, dfn, dfd, `F`, `parametric.p.value`, `resampled.p.value`) %>%
  mutate(
    dfn = as.integer(round(dfn, 0)),
    dfd = as.integer(round(dfd, 0)),
    `F` = round(`F`, 3),
    parametric.p.value = round(parametric.p.value, 3),
    resampled.p.value  = round(resampled.p.value, 3)
  )
#Save results to CSV file
write_csv(aovp_results_df, "./results/parental.performance.anova.csv")

# ==========================================================================
# Permutation ANOVA and pairwise comparison
# Note: pairwise tests ignore the random/nested structure and other factors
# ==========================================================================
# permuco::aovperm currently does not have a built-in emmeans interface, 
# and the authors confirm that it is not very suitable; 
# therefore, post hoc tests are usually performed using pairwise permutation tests for "simple effects"
aovp_posthoc <- function(data, var, factors) {

  # Build formula
  formula <- as.formula(paste(var, "~ gene + treatment1st + gene:treatment1st + Error(nest.id/(gene))"))

  # Perform permutation ANOVA
  set.seed(2024)
  perm_anova <- aovp(formula, data = data, 999)
  anova_result <- tidy(perm_anova) %>% mutate(variable = var, .before = 1)

  # Post hoc comparison for each main effect
  posthoc_results <- map(factors, function(factor) {
    posthoc_formula <- as.formula(paste(var, "~", factor))
    pairwisePermutationTest(posthoc_formula, data = data, method = "fdr")
  })

  names(posthoc_results) <- factors

  # If there is an interaction effect, perform simple main effect analysis
  if (length(factors) > 1) {
    interaction_term <- paste(factors, collapse = ":")

    if (interaction_term %in% anova_result$term) {
      # if (anova_result$'Pr(Prob)'[anova_result$term == interaction_term] < 0.05) {
      # Create interaction effect combinations
      data$interaction <- interaction(data[[factors[1]]], data[[factors[2]]])

      # Post hoc comparison for interaction effect
      posthoc_interaction <- pairwisePermutationTest(as.formula(paste(var, "~ interaction")),
                                                       data = data, method = "fdr")

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
    mutate(result_type = "anova_df", .before = 1)%>%
    mutate(p.value = as.numeric(p.value))

    posthoc_df <- subset$posthoc %>%
    mutate(result_type = "posthoc_df", .before = 1)%>%
    mutate(p.value = as.numeric(p.value))
  
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
# write_csv(anova_results.posthoc_df, "./results/parental.performance.posthoc.csv")

# Filter post hoc results for gene main effect
gene_main_effect <- anova_results.posthoc_df %>%
  filter(result_type == "posthoc_df", source == "gene") %>%
  mutate(groups = str_extract_all(Comparison, "G\\d+")) %>%
  unnest(groups) %>%
  group_by(variable, Comparison) %>%
  mutate(
    group1 = first(groups),
    group2 = last(groups)
  ) %>%
  ungroup() %>%
  distinct(variable, group1, group2, p.value) %>%
  mutate(p.value = as.numeric(p.value))

# Function to generate letter groupings
get_group_letters <- function(df) {
  mat <- matrix(1, nrow = 6, ncol = 6)
  rownames(mat) <- colnames(mat) <- sort(unique(c(df$group1, df$group2)))

  for (i in seq_len(nrow(df))) {
    g1 <- df$group1[i]
    g2 <- df$group2[i]
    pval <- as.numeric(df$p.value[i])
    mat[g1, g2] <- pval
    mat[g2, g1] <- pval
  }

  signif_mat <- mat < 0.1
  diag(signif_mat) <- FALSE

  letters <- multcompLetters(signif_mat)$Letters
  return(tibble(group = names(letters), group_letter = letters))
}

# Apply grouping by variable
group_results <- gene_main_effect %>%
  group_by(variable) %>%
  group_split() %>%
  map_df(~ {
    tibble(
      variable = unique(.x$variable),
      get_group_letters(.x)
    )
  })

# Print results
print(group_results)

# Save results to CSV file
write_csv(group_results, "./results/parental.posthoc.gene.csv")

# Comparison of Control vs Herbivory within each genotype G1–G6 under interaction extraction
genewise_treatment_comp <- anova_results.posthoc_df %>%
  filter(source == "interaction") %>%
  filter(str_detect(Comparison, "^G[1-6]\\.Control - G[1-6]\\.Herbivory = 0$")) %>%
  filter(str_extract(Comparison, "^G[1-6]") == str_extract(Comparison, "(?<=- )G[1-6]")) %>%
  mutate(
    gene = str_extract(Comparison, "^G[1-6]"),
    comparison = "Control vs Herbivory"
  ) %>%
  select(variable, gene, comparison, Stat, p.value, p.adjust) %>%
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE             ~ ""
  ))
# Save results to CSV file
write_csv(genewise_treatment_comp, "./results/parental.posthoc.interaction.csv")

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
                                0.3,
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

sig_group <- tibble::tribble(
  ~variable,     ~gene,   ~label,
 "mla" ,             "G1" , "a",
 "mla" ,             "G2" , "a",
 "mla" ,             "G3" , "a",
 "mla" ,             "G4" , "b",
 "mla" ,             "G5" , "a",
 "mla" ,             "G6" , "a",
 "no.leaves" ,       "G1" , "b",
 "no.leaves" ,       "G2" , "b",
 "no.leaves" ,       "G3" , "c",
 "no.leaves" ,       "G4" , "a",
 "no.leaves" ,       "G5" , "b",
 "no.leaves" ,       "G6" , "b",
 "no.ramets" ,       "G1" , "b",
 "no.ramets" ,       "G2" , "ab",
 "no.ramets" ,       "G3" , "c",
 "no.ramets" ,       "G4" , "a",
 "no.ramets" ,       "G5" , "ab",
 "no.ramets" ,       "G6" , "bc",
 "petiole.mass" ,    "G1" , "ab",
 "petiole.mass" ,    "G2" , "a",
 "petiole.mass" ,    "G3" , "bc",
 "petiole.mass" ,    "G4" , "c",
 "petiole.mass" ,    "G5" , "a",
 "petiole.mass" ,    "G6" , "abc",
 "ratio" ,           "G1" , "b",
 "ratio" ,           "G2" , "b",
 "ratio" ,           "G3" , "b",
 "ratio" ,           "G4" , "a",
 "ratio" ,           "G5" , "b",
 "ratio" ,           "G6" , "b",
 "root.mass" ,       "G1" , "ab",
 "root.mass" ,       "G2" , "a",
 "root.mass" ,       "G3" , "b",
 "root.mass" ,       "G4" , "ab",
 "root.mass" ,       "G5" , "a",
 "root.mass" ,       "G6" , "ab",
 # "sla" ,             "G1" , "",
 # "sla" ,             "G2" , "",
 # "sla" ,             "G3" , "",
 # "sla" ,             "G4" , "",
 # "sla" ,             "G5" , "",
 # "sla" ,             "G6" , "",
 "stem.length" ,     "G1" , "b",
 "stem.length" ,     "G2" , "ab",
 "stem.length" ,     "G3" , "b",
 "stem.length" ,     "G4" , "ab",
 "stem.length" ,     "G5" , "a",
 "stem.length" ,     "G6" , "ab",
 "stem.mass" ,       "G1" , "bcd",
 "stem.mass" ,       "G2" , "ab",
 "stem.mass" ,       "G3" , "cd",
 "stem.mass" ,       "G4" , "d",
 "stem.mass" ,       "G5" , "a",
 "stem.mass" ,       "G6" , "abc",
 "total.leaf.area" , "G1" , "a",
 "total.leaf.area" , "G2" , "a",
 "total.leaf.area" , "G3" , "ab",
 "total.leaf.area" , "G4" , "b",
 "total.leaf.area" , "G5" , "a",
 "total.leaf.area" , "G6" , "a",
 "total.leaf.mass" , "G1" , "abc",
 "total.leaf.mass" , "G2" , "ab",
 "total.leaf.mass" , "G3" , "c",
 "total.leaf.mass" , "G4" , "bc",
 "total.leaf.mass" , "G5" , "a",
 "total.leaf.mass" , "G6" , "abc",
 "total.mass" ,      "G1" , "ab",
 "total.mass" ,      "G2" , "a",
 "total.mass" ,      "G3" , "bc",
 "total.mass" ,      "G4" , "c",
 "total.mass" ,      "G5" , "a",
 "total.mass" ,      "G6" , "ab"
)

sig_inter <- tibble::tribble(
  ~variable,     ~gene,   ~label,
# "total.mass","G3","*",
# "total.leaf.mass","G2","*",
# "total.leaf.mass","G3","**",
# "petiole.mass","G3","*",
# "no.ramets","G2","*",
# "no.ramets","G3","**",
# "no.leaves","G2","*",
# "no.leaves","G3","**",
# "mla","G1","**",
# "sla","G1","***",
# "sla","G2","**",
# "sla","G3","**"
)

#  sapply(plot_info$variable, function(var) {
#  range(data_01gen[[var]], na.rm = TRUE)
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
    left_join(subset_data %>% group_by(gene) %>% summarize(y = max(mean) * 1.2), by = "gene")

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
  labs(fill = NULL) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12, family = "serif"),
    legend.key.size = unit(0.5, "cm")
  )

legend_grob <- cowplot::get_legend(legend_plot)

# Split plots into figure_2 and figure_3
figure_2_vars <- c("total.mass", "total.leaf.mass", "petiole.mass", "stem.mass", "root.mass", "ratio")
figure_3_vars <- c("no.ramets", "stem.length", "total.leaf.area", "no.leaves", "mla", "sla")

bold_tag_theme <- theme(plot.tag = element_text(face = 'bold'))

# Generate main figure
main_plot_2 <- wrap_plots(bar_plots[plot_info$variable %in% figure_2_vars], ncol = 2) +
    plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')& bold_tag_theme
main_plot_3 <- wrap_plots(bar_plots[plot_info$variable %in% figure_3_vars], ncol = 2)+
    plot_annotation(tag_levels = 'A', tag_prefix = '(', tag_suffix = ')')& bold_tag_theme

# Combine pictures
figure_2 <- ggdraw() +
  draw_plot(main_plot_2, x = 0, y = 0, width = 1, height = 1) +
  draw_grob(legend_grob, x = 0.80, y = 0.75, width = 0.22, height = 0.3)

figure_3 <- ggdraw() +
  draw_plot(main_plot_3, x = 0, y = 0, width = 1, height = 1) +
  draw_grob(legend_grob, x = 0.80, y = 0.75, width = 0.22, height = 0.3)

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

