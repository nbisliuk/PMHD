
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(lme4)
library(lmerTest)
library(emmeans)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(scales)

options(dplyr.summarise.inform = FALSE)
options(scipen = 999)

set.seed(2026)

pilot <- read_csv("G16.pilot.data.csv", show_col_types = FALSE) %>%
  rename(quality_idx = quality_index) %>%
  mutate(
    species = factor(species),
    location = factor(location),
    ecotr_id = factor(ecotr_id),
    tree_id = factor(tree_id)
  )

data <- read_csv("G16.outcome.data.csv", show_col_types = FALSE)

# Harmonize column name only if needed
if ("quality_index" %in% names(data)) {
  data <- data %>% rename(quality_idx = quality_index)
}

data <- data %>%
  mutate(
    quality_binary = if_else(quality_idx >= 65, 1L, 0L),
    climate = factor(climate),
    species = factor(species),
    tree_id = factor(tree_id),
    ecotr_id = factor(ecotr_id),
    location = factor(location)
  )

size <- read_csv("G16.size.data.csv", show_col_types = FALSE) %>%
  mutate(
    climate = factor(climate),
    species = factor(species),
    tree_id = factor(tree_id),
    ecotr_id = factor(ecotr_id)
  )

soil <- read_csv("G16.soil.data.csv", show_col_types = FALSE)



cat("Pears per climate\n")
print(table(data$climate))

cat("\nPears per tree summary\n")
print(summary(count(data, tree_id)$n))

quality_summary <- data %>%
  group_by(climate) %>%
  summarise(
    mean = mean(quality_idx, na.rm = TRUE),
    median = median(quality_idx, na.rm = TRUE),
    sd = sd(quality_idx, na.rm = TRUE)
  )

print(quality_summary)

par(mfrow = c(2, 3))

hist(count(data, tree_id)$n, main = "Histogram of Pears per tree", xlab = "")
hist(table(data$ecotr_id), main = "Histogram of Pears per ecotron", xlab = "")
hist(data$quality_idx, breaks = 20, xlab = "", main = "Histogram of Quality Index")

boxplot(quality_idx ~ climate, data, xlab = "", ylab = "Quality Index",
        main = "Quality Index per climate")

boxplot(quality_idx ~ species, data, xlab = "", ylab = "Quality Index",
        main = "Quality Index per species")

boxplot(quality_idx ~ ecotr_id, data, xlab = "", ylab = "Quality Index",
        main = "Quality Index per ecotron")

par(mfrow = c(1, 1))

# Optional exploratory check:
# Are pear counts systematically different between climates?
ecotron_pears <- data %>%
  group_by(ecotr_id) %>%
  summarise(
    climate = first(climate),
    pears = n()
  )

print(kruskal.test(pears ~ climate, data = ecotron_pears))
print(pairwise.wilcox.test(ecotron_pears$pears, ecotron_pears$climate, p.adjust.method = "BH"))



pilot_m <- lmer(quality_idx ~ species + location + (1 | ecotr_id) + (1 | tree_id),
                data = pilot, REML = TRUE)

print(summary(pilot_m))

pilot_vc <- as.data.frame(VarCorr(pilot_m))

sd_ecotron <- pilot_vc$sdcor[pilot_vc$grp == "ecotr_id"]
sd_tree <- pilot_vc$sdcor[pilot_vc$grp == "tree_id"]
sd_resid <- sigma(pilot_m)

if (length(sd_ecotron) == 0 || is.na(sd_ecotron)) sd_ecotron <- 0
if (length(sd_tree) == 0 || is.na(sd_tree)) sd_tree <- 0

mean_pears <- pilot %>%
  count(tree_id) %>%
  summarise(m = mean(n)) %>%
  pull(m)

sd_pears <- pilot %>%
  count(tree_id) %>%
  summarise(s = sd(n)) %>%
  pull(s)

if (is.na(sd_pears) || sd_pears == 0) sd_pears <- max(1, mean_pears * 0.10)

beta_intercept <- fixef(pilot_m)[["(Intercept)"]]
beta_species <- if ("speciesDoyenne" %in% names(fixef(pilot_m))) fixef(pilot_m)[["speciesDoyenne"]] else 0
beta_location <- if ("locationFR" %in% names(fixef(pilot_m))) fixef(pilot_m)[["locationFR"]] else 0

simulate_dataset <- function(n_ecotrons = 12,
                             trees_per_ecotron = 6,
                             climate_effects = c(0, 10, 10, 10)) {
  ecotron_df <- tibble(
    ecotr_id = factor(seq_len(n_ecotrons)),
    climate = factor(rep(1:4, length.out = n_ecotrons)),
    location = factor("BE"),
    ecotron_re = rnorm(n_ecotrons, 0, sd_ecotron)
  )

  tree_df <- ecotron_df %>%
    uncount(weights = trees_per_ecotron) %>%
    group_by(ecotr_id) %>%
    mutate(
      tree_id = factor(row_number()),
      species = factor(rep(c("Conference", "Doyenne"), length.out = n()),
                       levels = c("Conference", "Doyenne")),
      tree_re = rnorm(n(), 0, sd_tree),
      n_pears = pmax(1L, round(rnorm(n(), mean_pears, sd_pears)))
    ) %>%
    ungroup()

  pear_df <- tree_df %>%
    uncount(weights = n_pears) %>%
    group_by(ecotr_id, tree_id) %>%
    mutate(pear_id = row_number()) %>%
    ungroup() %>%
    mutate(
      quality_idx = beta_intercept +
        ifelse(species == "Doyenne", beta_species, 0) +
        ifelse(location == "FR", beta_location, 0) +
        climate_effects[as.integer(as.character(climate))] +
        ecotron_re +
        tree_re +
        rnorm(n(), 0, sd_resid)
    )

  pear_df
}

test_one_simulation <- function(dat, alpha = 0.05) {
  mod <- lmer(quality_idx ~ climate * species + (1 | ecotr_id) + (1 | tree_id), data = dat)
  cf <- summary(mod)$coefficients
  climate_rows <- grep("^climate", rownames(cf))
  if (length(climate_rows) == 0) return(FALSE)

  p_two <- cf[climate_rows, "Pr(>|t|)"]
  t_val <- cf[climate_rows, "t value"]

  # One-sided testing in the beneficial direction, with BH correction
  p_one <- ifelse(t_val > 0, p_two / 2, 1 - p_two / 2)
  p_adj <- p.adjust(p_one, method = "BH")

  any(p_adj < alpha)
}

estimate_power <- function(n_ecotrons_grid = seq(4, 20, by = 4),
                           trees_per_ecotron_grid = c(2, 4, 6),
                           nsim = 200,
                           alpha = 0.05,
                           climate_effects = c(0, 10, 10, 10)) {
  design_grid <- tidyr::crossing(
    n_ecotrons = n_ecotrons_grid,
    trees_per_ecotron = trees_per_ecotron_grid
  )

  design_grid$power <- mapply(
    function(n_ec, n_tree) {
      mean(replicate(nsim, {
        dat_sim <- simulate_dataset(
          n_ecotrons = n_ec,
          trees_per_ecotron = n_tree,
          climate_effects = climate_effects
        )
        test_one_simulation(dat_sim, alpha = alpha)
      }))
    },
    design_grid$n_ecotrons,
    design_grid$trees_per_ecotron
  )

  design_grid %>%
    mutate(
      total_trees = n_ecotrons * trees_per_ecotron,
      trees_with_10pct_dropout = ceiling(total_trees / 0.90),
      meets_80pct_power = power >= 0.80
    )
}

plot_power_curve <- function(power_tbl) {
  ggplot(power_tbl, aes(total_trees, power, color = factor(trees_per_ecotron))) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.80, linetype = 2) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
      x = "Total trees",
      y = "Estimated power",
      color = "Trees per ecotron",
      title = "Simulation-based power for quality-index analysis"
    ) +
    theme_minimal()
}

# Quality index

m1 <- lmer(quality_idx ~ climate * species + (1 | ecotr_id) + (1 | tree_id), data = data)
print(summary(m1))

m1_emm <- emmeans(m1, ~ climate | species)
print(m1_emm)
print(contrast(m1_emm, method = "trt.vs.ctrl", ref = 1))

ggplot(data, aes(climate, quality_idx)) +
  geom_boxplot(outlier.alpha = 0.25) +
  geom_jitter(width = 0.15, alpha = 0.15) +
  facet_wrap(~ species) +
  theme_minimal() +
  labs(title = "Quality index by climate and species", x = NULL, y = "Quality index")

# Binary quality outcome


m2 <- glmer(
  quality_binary ~ climate * species + (1 | ecotr_id) + (1 | tree_id),
  data = data,
  family = binomial
)

print(summary(m2))
print(exp(fixef(m2)))

m2_emm <- emmeans(m2, ~ climate | species, type = "response")
print(m2_emm)
print(contrast(m2_emm, method = "trt.vs.ctrl", ref = 1))

# Pear size over time

if (all(c("time", "size") %in% names(size))) {
  size_long <- size %>%
    transmute(
      pear_id = if ("pear_id" %in% names(size)) pear_id else row_number(),
      tree_id = tree_id,
      ecotr_id = ecotr_id,
      species = species,
      climate = climate,
      week = as.numeric(time),
      size_val = as.numeric(size)
    )
} else {
  size_long <- size %>%
    pivot_longer(
      cols = starts_with("week_"),
      names_to = "week",
      names_prefix = "week_",
      values_to = "size_val"
    ) %>%
    mutate(
      week = as.numeric(week),
      size_val = as.numeric(size_val)
    )
}

ggplot(size_long, aes(x = week, y = size_val, color = climate, group = interaction(tree_id, pear_id))) +
  geom_line(alpha = 0.25, na.rm = TRUE) +
  facet_wrap(~ climate) +
  theme_minimal() +
  labs(
    title = "Pear Size Growth by Climate",
    x = "Week Number",
    y = "Size",
    color = "Climate"
  )

m3 <- lmer(
  size_val ~ week * climate + species + (1 | ecotr_id) + (1 | tree_id / pear_id),
  data = size_long
)

print(summary(m3))

# ======================================================
# Soil PCA

soil_feature_map <- c(
  bulk_density = "Bulk Density",
  infiltration = "Infiltration",
  soil_porosity = "Soil Porosity",
  soil_depth = "Soil Depth",
  water_holding_capacity = "Water Holding Capacity",
  soil_nitrate = "Soil Nitrate",
  soil_ph = "Soil pH",
  phosphorus = "Phosphorus",
  potassium = "Potassium",
  earthworms = "Earthworms",
  microbial_biomass_c = "Microbial Biomass C",
  microbial_biomass_n = "Microbial Biomass N",
  particulate_organic_matter = "Particulate Organic Matter",
  mineralizable_n = "Mineralizable N",
  soil_enzymes = "Soil Enzymes",
  soil_respiration = "Soil Respiration",
  total_organic_carbon = "Total Organic Carbon"
)

names(soil) <- names(soil) |>
  gsub("[^A-Za-z0-9]+", "_", x = _) |>
  tolower()

soil$climate <- factor(soil$climate)

soil_numeric <- soil %>%
  select(any_of(names(soil_feature_map))) %>%
  mutate(across(everything(), as.numeric))

soil_numeric <- soil_numeric[, sapply(soil_numeric, function(x) sd(x, na.rm = TRUE) > 0), drop = FALSE]

if (ncol(soil_numeric) >= 2) {
  res.pca <- PCA(soil_numeric, graph = FALSE, scale.unit = TRUE)

  fviz_eig(res.pca, addlabels = TRUE)
  fviz_cos2(res.pca, choice = "var", axes = 1:2)

  var_pc_corr <- get_pca_var(res.pca)$cor
  var_pc_corr <- var_pc_corr[order(var_pc_corr[, 1], decreasing = TRUE), , drop = FALSE]
  corrplot(var_pc_corr)
  mtext("Correlations of variables with PCs", at = -3.5, line = 2.5)

  fviz_contrib(res.pca, choice = "var", axes = 1:2)

  var_pc_contrib <- get_pca_var(res.pca)$contrib
  var_pc_contrib <- var_pc_contrib[order(var_pc_contrib[, 1], decreasing = TRUE), , drop = FALSE]
  corrplot(var_pc_contrib, is.corr = FALSE)
  mtext("Contributions of variables to PCs", at = -4.5, line = 2.5)

  fviz_pca_ind(
    res.pca,
    geom.ind = "point",
    col.ind = soil$climate,
    addEllipses = TRUE,
    ellipse.type = "convex",
    legend.title = "Climate"
  )

  fviz_pca_var(
    res.pca,
    col.var = "contrib",
    repel = TRUE
  )

  fviz_pca_biplot(res.pca, repel = TRUE, col.ind = soil$climate)
} else {
  message("Too few non-constant soil variables available for PCA.")
}


results <- list(
  pilot_model = pilot_m,
  main_quality_model = m1,
  binary_quality_model = m2,
  size_model = m3
)
