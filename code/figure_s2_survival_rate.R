# Cleaning memory
cat("\014")
rm(list = ls())
gc()

# Load packages
library(broom)
library(ggpubr)
library(tidyverse)
library(skimr)
library(emmeans)
library(cowplot)   
library(grid)  
library(readxl)  

# set work directory
setwd("C:/Users/qianh/Desktop/R/chapter3")
getwd()

# Import data
data_02gen <- read_excel("./data/2022herbivory_data250325.xlsx", sheet = "2nd.gen")

data_02gen <- data_02gen %>% mutate(across(fw1st_mg:sla, ~ as.numeric(.)))

skim(data_02gen)
str(data_02gen)

# Factorize variables
data_02gen$treatment1st <- as.factor(data_02gen$treatment1st)
data_02gen$gene <- as.factor(data_02gen$gene)
data_02gen$rep <- as.factor(data_02gen$rep)
data_02gen$nest.id <- as.factor(data_02gen$nest.id)

levels(data_02gen$treatment1st)
levels(data_02gen$gene)

# Convert fresh weight (mg → g)
data_02gen$fw1st_g <- data_02gen$fw1st_mg / 1000

# ==========================================================================
# Survival analysis
# ==========================================================================
table(data_02gen$survival)

# Fit full model
model_full <- glm(survival ~ fw1st_g + gene + treatment1st + treatment1st:gene + (1 | nest.id/gene),
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
write_csv(chisq_results, "./results/chisq_results.offspring.survival.csv")

# Estimated marginal means
emm_gene <- emmeans(model_full, ~ gene)
emm_inter <- emmeans(model_full, ~ treatment1st | gene)

# Post-hoc tests
pairs(emm_gene, adjust = "holm")
pairs(emm_inter, adjust = "holm")
# Grouping (from Tukey results)
# G1 b | G2 b | G3 b | G4 a | G5 b | G6 b
# G2 shows * under treatment contrast

# ==========================================================================
# plotting
# ==========================================================================
# Summarize survival rate (%) and binomial SE (%)
subset_data <- data_02gen %>%
  dplyr::group_by(gene, treatment1st) %>%
  dplyr::summarise(
    survival_count = sum(survival == 1, na.rm = TRUE),
    total_count    = dplyr::n(),
    mean           = survival_count / total_count * 100,                     
    se             = sqrt((mean/100) * (1 - mean/100) / total_count) * 100,  
    .groups = "drop"
  )

color_values <- c("Control" = "#599CB4", "Herbivory" = "#C25759")

# Legend (extracted separately)
legend_plot <- ggplot(subset_data, aes(x = gene, y = mean, fill = treatment1st)) +
  geom_col(position = position_dodge(0.8), width = 0.8, color = "black") +
  scale_fill_manual(values = color_values) +
  theme_minimal(base_family = "serif") +
  labs(fill = "Parental") +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, family = "serif"),
    legend.text = element_text(size = 12, family = "serif"),
    legend.key.size = unit(0.5, "cm"),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  guides(fill = guide_legend(title.position = "left", title.hjust = 0.5, nrow = 1, byrow = TRUE))


legend_grob <- cowplot::get_legend(legend_plot)

# Group letters
sig_group <- tibble::tribble(
  ~gene, ~label,
  "G1", "b",
  "G2", "b",
  "G3", "b",
  "G4", "a",
  "G5", "b",
  "G6", "b"
)

# # Additional star mark
# star_data <- tibble(gene = "G2", y = 98, label = "*")

# Letter position
y_top_fallback <- subset_data %>%
  group_by(gene) %>%
  summarise(y = 105, .groups = "drop")

sig_data1 <- sig_group %>%
  left_join(y_top_fallback, by = "gene")

# Main plot
plot_sr <- ggplot(subset_data, aes(x = gene, y = mean, fill = treatment1st)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  geom_errorbar(
    aes(ymin = mean - se, ymax = mean + se),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_values) +
  labs(x = NULL, y = "Survival rate (%)", fill = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 100), oob = scales::rescale_none) +
  theme_minimal(base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(linewidth = 1, color = "black"),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black"),
    axis.ticks.x = element_line(linewidth = 1, color = "black"),
    axis.ticks.y = element_line(linewidth = 1, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    legend.position = "none"  
  )+  
  geom_text(data = sig_data1, aes(x = gene, y = y, label = label),
    inherit.aes = FALSE, size = 5, family = "serif") +
  # geom_text(data = star_data, aes(x = gene, y = y, label = label),
  #   inherit.aes = FALSE, size = 6, family = "serif") +
  coord_cartesian(ylim = c(0, 105), clip = "off")

# Combine main plot with legend
final_plot <- cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1.15)) + 
  cowplot::draw_plot(plot_sr + theme(legend.position = "none"),
                     x = 0, y = 0, width = 1, height = 1) +
  cowplot::draw_plot(legend_grob,
                     x = 0.05, y = 1.1, width = 0.96, height = 0.07,
                     hjust = 0, vjust = 1)

# Export figure
ggexport(final_plot, filename = "./results/figure_s2_survival_rate.png",
     width = 1000,
     height = 1000,
     pointsize = 12,
     res = 300)

