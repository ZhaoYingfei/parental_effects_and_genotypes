# Cleaning memory
cat("\014")
rm(list = ls())
gc()

# Load necessary packages
library(broom)
library(readxl)
library(tidyverse)
library(skimr)
library(permuco)

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

data_02gen$fw1st_g <- data_02gen$fw1st_mg / 1000

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

## =========================================================================
# permutation ANOVA: models with covariates
# ==========================================================================
# Function
perform_aovp <- function(var, data) {
  fml <- as.formula(paste(var, "~ fw1st_g + gene + treatment1st + gene:treatment1st + Error(nest.id/(gene))"))
  subset_data <- data[!is.na(data[[var]]), , drop = FALSE]
  set.seed(2024)
  fit <- permuco::aovperm(fml, data = subset_data, np = 999)
  summary(fit)  
}

aovp_results <- setNames(
  purrr::map(variables, ~ perform_aovp(.x, data_02gen)),
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
# Save results to CSV file
write_csv(aovp_results_df, "./results/offspring.performance.anova.withfw.csv")

