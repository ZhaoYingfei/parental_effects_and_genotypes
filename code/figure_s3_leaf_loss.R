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

# Fit models for emmeans
fit_glm_quasibin <- function(response_var, data) {
  glm(
    as.formula(paste(response_var, "~ fw1st_g + gene + treatment1st + treatment1st:gene")),
    family = quasibinomial(),
    data = data
  )
}

run_emmeans_blocks <- function(model, response_var, data) {
  fw_mean <- mean(data$fw1st_g, na.rm = TRUE)

  # (A) Overall gene effect
  emm_gene <- emmeans(
    model,
    ~ gene,
    weights = "proportional",
    at = list(fw1st_g = fw_mean)
  )
  emm_gene_resp <- summary(emm_gene, type = "response")
  pairs_gene    <- summary(pairs(emm_gene, adjust = "holm"))
  cld_gene      <- multcomp::cld(emm_gene, adjust = "holm", Letters = letters)

  # (B) Treatment within gene
  emm_trt_by_gene <- emmeans(
    model,
    ~ treatment1st | gene,
    at = list(fw1st_g = fw_mean)
  )
  emm_trt_by_gene_resp <- summary(emm_trt_by_gene, type = "response")
  pairs_trt_by_gene <- pairs(emm_trt_by_gene, adjust = "holm")


  # Return all results
  list(
    response            = response_var,
    emm_gene_resp       = emm_gene_resp,
    pairs_gene          = pairs_gene,
    cld_gene            = cld_gene,
    emm_trt_by_gene_resp = emm_trt_by_gene_resp,
    pairs_trt_by_gene   = pairs_trt_by_gene
  )
}

model_loss_fw <- fit_glm_quasibin("loss.fw", data_02gen)
model_loss_la <- fit_glm_quasibin("loss.la", data_02gen)

emm_out_fw <- run_emmeans_blocks(model_loss_fw, "loss.fw", data_02gen)
emm_out_la <- run_emmeans_blocks(model_loss_la, "loss.la", data_02gen)

# emm_out_fw$emm_gene_resp       
# emm_out_fw$pairs_gene       
emm_out_fw$cld_gene       
# emm_out_fw$emm_trt_by_gene_resp
emm_out_fw$pairs_trt_by_gene  

# fw Grouping (from Tukey results)
# G1 ab | G2 ab | G3 b | G4 a | G5 b | G6 b

# emm_out_la$emm_gene_resp  
# emm_out_la$pairs_gene        
emm_out_la$cld_gene          
# emm_out_la$emm_trt_by_gene_resp 
emm_out_la$pairs_trt_by_gene  

# la Grouping (from Tukey results)
# G1 ab | G2 ab | G3 ab | G4 a | G5 b | G6 ab

# ==========================================================================
# plotting loss data
# ==========================================================================
# Summarize means and SE (percentage)
subset_data <- data_02gen %>%
  dplyr::group_by(gene, treatment1st) %>%
  dplyr::summarise(
    mean_fw = mean(loss.fw, na.rm = TRUE) * 100,                         # 转百分比
    se_fw   = (sd(loss.fw, na.rm = TRUE) / sqrt(dplyr::n())) * 100,      # 标准误也转百分比
    mean_la = mean(loss.la, na.rm = TRUE) * 100,
    se_la   = (sd(loss.la, na.rm = TRUE) / sqrt(dplyr::n())) * 100,
    .groups = "drop"
  )

color_values <- c("Control" = "#599CB4", "Herbivory" = "#C25759")

# Extract legend
legend_plot <- ggplot(subset_data, aes(x = gene, y = mean_fw, fill = treatment1st)) +
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
  )+
  guides(fill = guide_legend(title.position = "left", title.hjust = 0.5, nrow = 1, byrow = TRUE))

legend_grob <- cowplot::get_legend(legend_plot)

sig_group_fw <- tibble::tribble(
  ~gene, ~label,
  "G1", "ab",
  "G2", "ab",
  "G3", "b",
  "G4", "a",
  "G5", "b",
  "G6", "b"
)
sig_group_la <- tibble::tribble(
  ~gene, ~label,
  "G1", "ab",
  "G2", "ab",
  "G3", "ab",
  "G4", "a",
  "G5", "b",
  "G6", "ab"
)

y_top_fallback <- subset_data %>%
  group_by(gene) %>%
  summarise(y = 105, .groups = "drop")

sig_data_fw <- sig_group_fw %>%
  left_join(y_top_fallback, by = "gene")
sig_data_la <- sig_group_la %>%
  left_join(y_top_fallback, by = "gene")

# Plot loss of fresh weight
plot_fw <- ggplot(subset_data, aes(x = gene, y = mean_fw, fill = treatment1st)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  geom_errorbar(
    aes(ymin = mean_fw - se_fw, ymax = mean_fw + se_fw),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_values) +
  labs(x = NULL, y = "Loss of fresh mass (%)", fill = NULL) +   
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
  ) +
  geom_text(data = sig_data_fw, aes(x = gene, y = y, label = label),
            inherit.aes = FALSE, size = 5, family = "serif") +
  coord_cartesian(clip = "off")

# Plot loss of leaf area
plot_la <- ggplot(subset_data, aes(x = gene, y = mean_la, fill = treatment1st)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.8, color = "black") +
  geom_errorbar(
    aes(ymin = mean_la - se_la, ymax = mean_la + se_la),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  scale_fill_manual(values = color_values) +
  labs(x = NULL, y = "Loss of fresh mass (%)", fill = NULL) +   
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
  ) +
  geom_text(data = sig_data_la, aes(x = gene, y = y, label = label),
            inherit.aes = FALSE, size = 5, family = "serif") +
  coord_cartesian(clip = "off")

# Combine plots
combined_plots <- plot_grid(
  plot_fw,
  plot_la,
  ncol = 2,
  labels = c("A", "B"),  
  align = "v",
  label_fontfamily = "serif" 
)
# Combine main plot with legend
final_plot <- cowplot::ggdraw(xlim = c(0, 1), ylim = c(0, 1.15)) + 
  cowplot::draw_plot(combined_plots + theme(legend.position = "none"),
                     x = 0, y = 0, width = 1, height = 1) +
  cowplot::draw_plot(legend_grob,
                     x = 0.05, y = 1.1, width = 0.96, height = 0.07,
                     hjust = 0, vjust = 1)


# Export figure
ggexport(final_plot, filename = "./results/figure_s3_leaf_loss.png",  
         width = 2000,
         height = 1000,
         pointsize = 12,
         res = 300)
