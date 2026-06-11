# ============================================================================
# Quantifying L2 vowel categories in acoustic space
# Distinctiveness and nativelikeness metrics in R
# ============================================================================
#
# This script reads the experimental dataset, describes the structure of the
# stimuli, visualises vowel clouds, computes vowel-distance metrics with
# `voweldist()`, explores covariance size, normalises selected output measures,
# and fits one example mixed-effects model.
#
# Optional setup:
# rm(list = ls())
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# Run the setwd() line only if you want the working directory to be the folder
# containing this script.

# ============================================================================
# 1. Packages and data import
# ============================================================================

suppressPackageStartupMessages({
  library(glmmTMB)
  library(bestNormalize)
  library(tidyverse)
  library(readxl)
})

source("https://raw.githubusercontent.com/ebrenc/rstats/refs/heads/main/voweldist.R")
source("https://raw.githubusercontent.com/ebrenc/rstats/refs/heads/main/detcov.R")
source("https://raw.githubusercontent.com/ebrenc/rstats/refs/heads/main/mahalanobis_outliers.R")
source("https://raw.githubusercontent.com/ebrenc/rstats/refs/heads/main/glmmTable.R")

# Read the data and set the intended order of the main categorical variables.
# The order matters for tables, plots, and model reference levels.
df <- read_excel("data.xlsx", guess_max = 15000) %>%
  tibble() %>%
  mutate(
    across(where(is.character), ~ na_if(str_squish(.x), "")),
    Task = factor(Task, levels = c("st", "sp", "wr")),
    WordAge = factor(WordAge, levels = c("Old", "New")),
    WordType = factor(WordType, levels = c("Nonword", "Word")),
    Trained = factor(Trained, levels = c("Untrained", "Trained")),
    ItemPair = factor(ItemPair)
  )

# ============================================================================
# 2. Descriptive checks
# ============================================================================

# ---------------------------------------------------------------------------
# 2.1 Number of participants by task and group/condition
# ---------------------------------------------------------------------------
# One speaker may appear in more than one task, so the count is calculated
# separately for each Task × Group × Condition combination.
df %>%
  distinct(Task, Group, Condition, Speaker) %>%
  group_by(Task, Group, Condition) %>%
  count() %>%
  pivot_wider(names_from = Task, values_from = n) %>%
  print(n = Inf)

# ---------------------------------------------------------------------------
# 2.2 Words by task and experimental condition
# ---------------------------------------------------------------------------
# This table lists the lexical items used in each combination of Task, Vowel,
# WordAge, Trained, and WordType, excluding native speakers from the stimulus
# summary. The two vowels are displayed in separate columns.
df %>%
  filter(Group != "NS") %>%
  distinct(Task, Vowel, Trained, WordAge, WordType, Word) %>%
  group_by(Task, Vowel, Trained, WordAge, WordType) %>%
  summarise(
    Words = paste(sort(unique(Word)), collapse = ", "),
    n_words = n_distinct(Word),
    .groups = "drop"
  ) %>%
  arrange(Task, Vowel, Trained, WordAge, WordType) %>%
  pivot_wider(
    names_from = Vowel,
    values_from = c(Words, n_words),
    names_glue = "{Vowel}_{.value}"
  ) %>%
  mutate(
    n_words = if_else(æ_n_words == ʌ_n_words, æ_n_words, NA_integer_)
  ) %>%
  select(Task, Trained, WordAge, WordType, n_words, æ_Words, ʌ_Words) %>%
  print(n = Inf)

# ---------------------------------------------------------------------------
# 2.3 Contrasts represented in the data
# ---------------------------------------------------------------------------
# st, sp, wr -> Time
# sp, wr     -> within Word: WordAge (Old vs New)
# wr         -> within Old Word: Trained (Trained vs Untrained)
# wr         -> WordType (New -Untrained- Word vs Nonword)

# ============================================================================
# 3. Visual inspection of vowel clouds
# ============================================================================
# The plots show 50% normal ellipses for the two vowels in the acoustic space
# defined by B2B1 on the x-axis and B1B0 on the y-axis.

# Task st: vowel clouds by Condition and Time.
df %>%
  filter(Task == "st") %>%
  ggplot(aes(x = B2B1, y = B1B0, colour = Vowel)) +
  facet_grid(Condition ~ Time) +
  stat_ellipse(type = "norm", level = 0.5, linewidth = 0.9) +
  coord_equal()

# Task sp: vowel clouds by Condition, WordAge, and Time.
df %>%
  filter(Task == "sp") %>%
  ggplot(aes(x = B2B1, y = B1B0, colour = Vowel)) +
  facet_grid(Condition ~ WordAge * Time) +
  stat_ellipse(type = "norm", level = 0.5, linewidth = 0.9) +
  coord_equal()

# Task wr: vowel clouds by Condition and lexical/training factors.
df %>%
  filter(Task == "wr") %>%
  ggplot(aes(x = B2B1, y = B1B0, colour = Vowel)) +
  facet_grid(Condition ~ WordType * WordAge * Trained * Time) +
  stat_ellipse(type = "norm", level = 0.5, linewidth = 0.9) +
  coord_equal()

# ============================================================================
# 4. Split the dataset by task
# ============================================================================
# This creates three objects in the global environment named st, sp, and wr.
# Each object contains the rows for one task.

st = df %>% filter(Task == "st")
sp = df %>% filter(Task == "sp")

# ============================================================================
# 5. Distance metrics for task st
# ============================================================================
# Two types of distances are calculated:
# - within-group distances: inter_group = FALSE
# - between-group distances to the Native reference group: inter_group = TRUE

st_wit <- st %>%
  voweldist(
    dependent_vars = c(B1B0, B2B1),
    inter_group = FALSE,
    contrast_var = Vowel,
    condition_vars = Time,
    speaker_var = Speaker,
    compute_euclidean = TRUE,
    compute_mahalanobis = TRUE,
    compute_pillai = TRUE,
    compute_bhatt = TRUE
  )

st_bet <- st %>%
  voweldist(
    dependent_vars = c(B1B0, B2B1),
    inter_group = TRUE,
    contrast_var = Vowel,
    condition_vars = Time,
    speaker_var = Speaker,
    group_var = Condition,
    reference_group = "Native",
    compute_euclidean = TRUE,
    compute_mahalanobis = TRUE,
    compute_pillai = TRUE,
    compute_bhatt = TRUE
  )

# Join within- and between-group outputs for st.
st <- left_join(st_wit, st_bet)
rm(list = ls(pattern = "^st_(wit|bet)$"))

# Summarise the joined st distances.
st %>%
  group_by(Condition, Time, Vowel) %>%
  mutate(n_Speaker = n_distinct(Speaker), .groups = "drop") %>%
  group_by(Condition, n_Speaker, Time, Vowel) %>%
  summarise(
    n_Items = unique(n() / n_Speaker),
    cloud_wit = unique(cloud_wit),
    cloud_bet = unique(cloud_bet),
    across(matches("^dist_|^detcov") & !matches("_cloud"), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename_with(~ str_remove(.x, "^dist_"), starts_with("dist_")) %>%
  arrange(Condition, Time, Vowel) %>%
  print(n = Inf)

# rm(st)

# ============================================================================
# 6. Distance metrics for task sp
# ============================================================================

sp_wit <- sp %>%
  voweldist(
    dependent_vars = c(B1B0, B2B1),
    inter_group = FALSE,
    contrast_var = Vowel,
    condition_vars = c(Time, WordAge),
    speaker_var = Speaker,
    compute_euclidean = TRUE,
    compute_mahalanobis = TRUE,
    compute_pillai = TRUE,
    compute_bhatt = TRUE
  )

sp_wit %>%
  group_by(Condition, WordAge, cloud_wit) %>%
  count() %>%
  print(n = Inf)

# set.seed(123)
# sample_n <- 2500
# spSample <- sp %>%
#   mutate(.row_id = row_number()) %>%
#   slice_sample(n = sample_n) %>%
#   arrange(.row_id) %>%
#   select(-.row_id)
# spSample_wit <- spSample %>%
#   voweldist(
#     dependent_vars = c(B1B0, B2B1),
#     inter_group = FALSE,
#     contrast_var = Vowel,
#     condition_vars = c(Time, WordAge),
#     speaker_var = Speaker,
#     compute_euclidean = TRUE,
#     compute_mahalanobis = TRUE,
#     compute_pillai = TRUE,
#     compute_bhatt = TRUE
#   )
# 
# spSample_wit %>%
#   group_by(Condition, WordAge, cloud_wit) %>%
#   count() %>%
#   print(n = Inf)

sp_bet <- sp %>%
  voweldist(
    dependent_vars = c(B1B0, B2B1),
    inter_group = TRUE,
    contrast_var = Vowel,
    condition_vars = c(Time, WordAge),
    speaker_var = Speaker,
    group_var = Condition,
    reference_group = "Native",
    compute_euclidean = TRUE,
    compute_mahalanobis = TRUE,
    compute_pillai = TRUE,
    compute_bhatt = TRUE
  )

sp <- left_join(sp_wit, sp_bet)
rm(list = ls(pattern = "^sp_(wit|bet)$"))

# Summarise task-sp distances by Condition, WordAge, Time, and Vowel.
sp %>%
  group_by(Condition, WordAge, Time, Vowel) %>%
  mutate(n_Speaker = n_distinct(Speaker), .groups = "drop") %>%
  group_by(Condition, n_Speaker, WordAge, Time, Vowel) %>%
  summarise(
    n_Items = unique(n() / n_Speaker),
    cloud = unique(cloud_wit),
    across(matches("^dist_|^detcov") & !matches("_cloud"), \(x) mean(x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  rename_with(~ str_remove(.x, "^dist_"), starts_with("dist_")) %>%
  arrange(Condition, WordAge, Time, Vowel) %>%
  print(n = Inf)

# ============================================================================
# 8. Determinant of covariance for task sp
# ============================================================================
# The determinant of the covariance matrix is used here as a measure of the
# size/spread of the vowel cloud. The original script calculates it for æ.

sp_d <- sp %>%
  filter(Vowel == "æ") %>%
  group_by(Condition, Speaker, Vowel, Time, WordAge) %>%
  summarise(d = detcov(pick(B1B0, B2B1)), .groups = "drop")

# Mean determinant by Condition, Vowel, Time, and WordAge.
sp_d %>%
  group_by(Condition, Vowel, Time, WordAge) %>%
  summarise(d = mean(d, na.rm = TRUE), .groups = "drop") %>%
  arrange(d) %>%
  print(n = Inf)

# Plot only the minimum and maximum determinant cases.
df %>%
  filter(Task == "sp") %>%
  inner_join(
    sp_d %>%
      group_by(Condition, Vowel, Time, WordAge) %>%
      summarise(
        d = mean(d, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      {
        bind_rows(
          slice_min(., d, n = 1, with_ties = TRUE) %>%
            mutate(extreme = "min"),
          slice_max(., d, n = 1, with_ties = TRUE) %>%
            mutate(extreme = "max")
        )
      } %>%
      mutate(extreme = factor(extreme, levels = c("min", "max"))) %>%
      arrange(d) %>%
      select(Condition, Vowel, Time, WordAge, extreme, d),
    by = c("Condition", "Vowel", "Time", "WordAge")
  ) %>%
  ggplot(aes(x = B2B1, y = B1B0)) +
  facet_grid(. ~ extreme) +
  geom_point() +
  stat_ellipse(type = "norm", level = 0.5, linewidth = 0.9) +
  coord_equal()

# rm(list = ls(pattern = "^sp_d"))

# ============================================================================
# 9. Normalisation and exploratory plots
# ============================================================================

# library(ggpubr)

# Ordered quantile normalisation. The original script applies it to all
# numeric variables after converting Item and ItemN to character.
ordernormal <- function(var) {
  return(orderNorm(var)$x.t)
}

sp_norm <- sp %>%
  mutate(Item = as.character(Item)) %>%
  mutate(ItemN = as.character(ItemN)) %>%
  mutate(across(where(is.numeric), ~ ordernormal(.x)))

# Example boxplot for a within-group distance.
sp_norm %>%
  ggplot(aes(Condition, dist_wit_bha, fill = Time)) +
  geom_boxplot()

# ============================================================================
# 10. Example mixed model
# ============================================================================
# The model below tests whether the normalised within-group Bhattacharyya
# distance varies as a function of Condition, WordAge, and Time.

model <- sp_norm %>%
  filter(Group == "ES") %>%
  mutate(Condition = Condition %>% factor(levels = c("Ctrl", "Exp1", "Exp2"))) %>%
  glmmTMB(dist_wit_bha ~ Condition * WordAge * Time + (1 + WordAge | Speaker) + (1 | ItemPair), data = .)

model %>%
  glmmTable("mixed_anova_table.html")
