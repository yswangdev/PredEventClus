# ==============================================================================
# Visualization: General Plotting Functions
# ==============================================================================
#
# This script contains code for generating visualizations including:
#   - Kaplan-Meier survival curves by cluster
#   - MDS (Multi-Dimensional Scaling) plots for cluster visualization
#   - Covariate distribution plots
#
# Usage:
#   - Load simulation or real-world results
#   - Generate KM curves and MDS plots
#   - Customize plots for publication
#
# ==============================================================================

# Required libraries
library(survival)
library(dplyr)
library(survminer)
library(cluster)
library(ggplot2)

# ------------------------------------------------------------------------------
# User configuration
# ------------------------------------------------------------------------------

# Examples:
# results_rds_path <- "all_results_no_clus_s1_n_500.rds"       # Scenario 1
# results_rds_path <- "all_results_s2_n_500.rds"               # Scenario 2
# results_rds_path <- "all_results_s3_less_n_500.rds"          # Scenario 3
# results_rds_path <- "all_results_s4_n_500.rds"               # Scenario 4
results_rds_path <- "all_results_no_clus_s1_n_500.rds"
scenario_label <- ifelse(
  grepl("^all_results_no_clus_s1_n_", results_rds_path),
  "Scenario 1",
  ifelse(
    grepl("^all_results_s2_n_", results_rds_path),
    "Scenario 2",
    ifelse(
      grepl("^all_results_s3_", results_rds_path),
      "Scenario 3",
      ifelse(
        grepl("^all_results_s4_n_", results_rds_path),
        "Scenario 4",
        "Simulation Scenario"
      )
    )
  )
)
sim_index <- 1

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

extract_full_data <- function(all_results, sim_index = 1) {
  if (!is.list(all_results) || length(all_results) < sim_index) {
    stop("Invalid all_results object or sim_index out of range.")
  }
  if (is.null(all_results[[sim_index]]$full_data)) {
    stop("full_data not found in selected simulation result.")
  }
  all_results[[sim_index]]$full_data
}

make_cluster_labels <- function(cluster_vec) {
  levels_vec <- sort(unique(cluster_vec))
  paste("Cluster", levels_vec)
}

prepare_mixed_data <- function(full_data) {
  req_cols <- c("age", "hemoglobin", "tumor_burden", "gender", "race", "kps", "disease_sites")
  missing_cols <- setdiff(req_cols, names(full_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns for MDS:", paste(missing_cols, collapse = ", ")))
  }

  full_data %>%
    dplyr::select(all_of(req_cols)) %>%
    mutate(across(c(gender, race, kps, disease_sites), as.factor))
}

# ------------------------------------------------------------------------------
# Load selected simulation output
# ------------------------------------------------------------------------------

all_results <- readRDS(results_rds_path)
full_data <- extract_full_data(all_results, sim_index = sim_index)

# Ensure cluster is a factor for consistent legends/colors
if (!("cluster" %in% names(full_data))) {
  stop("Column `cluster` not found in full_data.")
}
full_data$cluster <- factor(full_data$cluster)

# ------------------------------------------------------------------------------
# 1) Kaplan-Meier curve by true cluster
# ------------------------------------------------------------------------------

surv_obj <- Surv(time = full_data$time, event = full_data$status)
fit <- survfit(surv_obj ~ cluster, data = full_data)
cluster_labs <- make_cluster_labels(as.integer(as.character(full_data$cluster)))

km_plot <- ggsurvplot(
  fit,
  data = full_data,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  legend.title = "True Cluster",
  legend.labs = cluster_labs,
  palette = "Set1",
  title = scenario_label,
  xlab = "Time",
  ylab = "Survival Probability"
)

print(km_plot)

# ------------------------------------------------------------------------------
# 2) MDS plot for mixed covariates
# ------------------------------------------------------------------------------

mixed_data <- prepare_mixed_data(full_data)
gower_dist <- daisy(mixed_data, metric = "gower")
mds_coords <- cmdscale(gower_dist, k = 2)

mds_df <- as.data.frame(mds_coords)
colnames(mds_df) <- c("Dim1", "Dim2")
mds_df$cluster <- full_data$cluster

mds_plot <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = paste0(scenario_label, ": MDS of Mixed Covariates"),
    x = "MDS Dimension 1",
    y = "MDS Dimension 2",
    color = "True Cluster"
  ) +
  theme_minimal()

print(mds_plot)
