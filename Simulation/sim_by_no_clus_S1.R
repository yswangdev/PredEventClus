# ==============================================================================
# Simulation Study: Scenario 1 (No True Clustering, S1-Matched Parameters)
# ==============================================================================
#
# This script evaluates milestone event time prediction when there is no true
# cluster structure. A single cohort (one true cluster) is generated, then
# "None" (original method) is compared against clustering-based methods.
#
# Output:
#   - all_results_no_clus_s1_n_[n_total].rds
#   - predict_timepoints_no_clus_s1_n_[n_total].csv
#   - predict_if_no_clus_s1_n_[n_total].csv
#   - compare_vs_none_timepoints_no_clus_s1_n_[n_total].csv
#   - compare_vs_none_if_no_clus_s1_n_[n_total].csv
#
# ==============================================================================

library(simsurv)
library(flexsurv)
library(dplyr)
library(cluster)
library(clustMixType)
library(kamila)
library(parallel)
library(future.apply)
library(MASS)

# Source shared functions
source("simulation_functions.R")

# ==============================================================================
# Simulation Parameters
# ==============================================================================

n_total <- 500
n_simulations <- 1
FU_time <- 18
target_event <- 0.7 * n_total
n_cores <- detectCores()

# "None" is the original method (no clustering); others are clustering methods.
clustering_methods <- c("None", "PAM", "Hierarchical", "K-prototypes", "Kamila")
model_types <- c("Exponential", "Weibull")

# Scenario 1 parameters for a single homogeneous cohort.
betas_common <- c(age = 0.03, gender = 0.2, hemoglobin = -0.3, tumor_burden = 0.03)
lambda_common <- 0.15
gamma_common <- 1

# Generate patient covariates from one distribution (single true cluster).
simulate_covariates <- function(n) {
  data.frame(
    id = seq_len(n),
    cluster = 1,
    age = rnorm(n, mean = 65, sd = 5),
    gender = rbinom(n, 1, 0.7),
    race = sample(c("White", "Black", "Other"), n, replace = TRUE, prob = c(0.8, 0.1, 0.1)),
    kps = sample(1:4, n, replace = TRUE, prob = c(0.4, 0.4, 0.1, 0.1)),
    hemoglobin = rnorm(n, mean = 11, sd = 1),
    disease_sites = sample(1:5, n, replace = TRUE, prob = c(0.1, 0.1, 0.2, 0.2, 0.4)),
    tumor_burden = rnorm(n, mean = 70, sd = 8)
  )
}

# Main simulation function
sim_func <- function(i, clustering_methods, model_types, targetev, sims, n_clus) {
  covs_all <- simulate_covariates(n_total)

  # Simulate event times for one homogeneous cohort.
  sim_all <- simsurv(dist = "exponential", lambdas = lambda_common, betas = betas_common, x = covs_all)

  # Exponential dropout censoring
  drop_time <- rexp(n_total, 0.01)
  sim_all$dropout <- ifelse(drop_time < sim_all$eventtime, 1, 0)
  sim_all$status <- ifelse(sim_all$dropout == 1, 0, sim_all$status)
  sim_all$eventtime <- pmin(sim_all$eventtime, drop_time)
  sim_all$id <- 1:n_total
  colnames(sim_all)[2] <- "time"

  # Staggered entry + administrative censoring
  sim_all$entry_time <- runif(n_total, 0, 2)
  sim_all$ttf <- sim_all$entry_time + sim_all$time
  sim_all$status <- ifelse(sim_all$ttf > FU_time, 0, sim_all$status)
  sim_all$ttf <- ifelse(sim_all$ttf < FU_time, sim_all$ttf, FU_time)
  sim_all$time <- ifelse(sim_all$ttf < FU_time, sim_all$time, FU_time)

  full_data <- left_join(sim_all, covs_all, by = "id") %>%
    mutate(
      gender = factor(gender, labels = c("Female", "Male")),
      race = factor(race),
      kps = factor(kps, ordered = TRUE),
      disease_sites = factor(disease_sites, ordered = TRUE)
    )
  mixed_data <- full_data %>% dplyr::select(age, hemoglobin, tumor_burden, gender, race, kps, disease_sites)
  full_data$sim_id <- i

  rolling_list <- list()
  rolling_predtime_list <- list()
  if_list <- list()
  if_predtimes_list <- list()

  for (cm in clustering_methods) {
    tryCatch({
      rolling_result <- rolling_milestone_prediction(
        full_data = full_data,
        mixed_data = mixed_data,
        clustering_method = cm,
        modeltypes = model_types,
        criterion = "BIC",
        targetev = targetev,
        time_points = as.integer(seq.int(2, 14, length.out = 6)),
        sims = sims,
        n_clus = n_clus,
        sim_id = i
      )

      info_result <- info_frac_prediction(
        full_data = full_data,
        mixed_data = mixed_data,
        clustering_method = cm,
        modeltypes = model_types,
        criterion = "BIC",
        targetev = targetev,
        info_frac = c(0.4, 0.6, 0.75, 0.8),
        sims = sims,
        n_clus = n_clus,
        sim_id = i
      )

      rolling_list[[length(rolling_list) + 1]] <- rolling_result$pred_results
      rolling_predtime_list[[length(rolling_predtime_list) + 1]] <- rolling_result$pred_times
      if_list[[length(if_list) + 1]] <- info_result$pred_results
      if_predtimes_list[[length(if_predtimes_list) + 1]] <- info_result$pred_times
    }, error = function(e) {
      message(paste("Error in", cm, ":", e$message))
      NULL
    })
  }

  list(
    rolling = do.call(rbind, rolling_list),
    rolling_predtimes = do.call(rbind, rolling_predtime_list),
    info = do.call(rbind, if_list),
    info_predtimes = do.call(rbind, if_predtimes_list),
    full_data = full_data
  )
}

# ==============================================================================
# Run Simulation
# ==============================================================================

plan(multisession, workers = n_cores)
all_results <- future_lapply(
  1:n_simulations,
  sim_func,
  clustering_methods = clustering_methods,
  model_types = model_types,
  targetev = target_event,
  sims = 1000,
  n_clus = 1,
  future.seed = TRUE
)

rolling_results <- bind_rows(lapply(all_results, `[[`, "rolling"))
info_results <- bind_rows(lapply(all_results, `[[`, "info"))

saveRDS(all_results, paste0("all_results_no_clus_s1_n_", n_total, ".rds"))
gc()

# ==============================================================================
# Summaries
# ==============================================================================

summary_rolling <- rolling_results %>%
  mutate(
    mean_rae = as.numeric(mean_rae),
    mean_err = as.numeric(mean_err),
    mean_ae = as.numeric(mean_ae)
  )

inf_sim_df <- summary_rolling %>%
  group_by(Sim) %>%
  summarise(has_inf = any(is.infinite(c_across(c(mean_rae, mean_err, mean_ae)))), .groups = "drop") %>%
  filter(has_inf)

n_total_sim <- summary_rolling$Sim %>% unique() %>% length()
n_inf_sim <- nrow(inf_sim_df)

if (n_inf_sim / n_total_sim < 0.05) {
  summary_rolling <- summary_rolling %>%
    filter(!(Sim %in% inf_sim_df$Sim))
}

summary_by_time <- summary_rolling %>%
  group_by(Landmark_Time, Clustering_Method) %>%
  summarise(
    Mean_RAE = mean(mean_rae, na.rm = TRUE),
    Mean_ERR = mean(mean_err, na.rm = TRUE),
    MEAN_AE = mean(mean_ae, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_by_time, paste0("predict_timepoints_no_clus_s1_n_", n_total, ".csv"))

summary_if <- info_results %>%
  mutate(
    mean_rae = as.numeric(mean_rae),
    mean_err = as.numeric(mean_err),
    mean_ae = as.numeric(mean_ae)
  )

inf_sim_df <- summary_if %>%
  group_by(Sim) %>%
  summarise(has_inf = any(is.infinite(c_across(c(mean_rae, mean_err, mean_ae)))), .groups = "drop") %>%
  filter(has_inf)

n_total_sim <- summary_if$Sim %>% unique() %>% length()
n_inf_sim <- nrow(inf_sim_df)

if (n_inf_sim / n_total_sim < 0.05) {
  summary_if <- summary_if %>%
    filter(!(Sim %in% inf_sim_df$Sim))
}

summary_by_IF <- summary_if %>%
  group_by(IF, Clustering_Method) %>%
  summarise(
    Mean_RAE = mean(mean_rae, na.rm = TRUE),
    Mean_ERR = mean(mean_err, na.rm = TRUE),
    MEAN_AE = mean(mean_ae, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_by_IF, paste0("predict_if_no_clus_s1_n_", n_total, ".csv"))

# ==============================================================================
# Explicit comparison: original (None) vs clustering methods
# ==============================================================================

none_by_time <- summary_by_time %>%
  filter(Clustering_Method == "None") %>%
  dplyr::select(Landmark_Time, None_Mean_RAE = Mean_RAE, None_Mean_ERR = Mean_ERR, None_MEAN_AE = MEAN_AE)

compare_vs_none_time <- summary_by_time %>%
  filter(Clustering_Method != "None") %>%
  left_join(none_by_time, by = "Landmark_Time") %>%
  mutate(
    Delta_RAE_vs_None = Mean_RAE - None_Mean_RAE,
    Delta_ERR_vs_None = Mean_ERR - None_Mean_ERR,
    Delta_AE_vs_None = MEAN_AE - None_MEAN_AE
  )

write.csv(compare_vs_none_time, paste0("compare_vs_none_timepoints_no_clus_s1_n_", n_total, ".csv"))

none_by_if <- summary_by_IF %>%
  filter(Clustering_Method == "None") %>%
  dplyr::select(IF, None_Mean_RAE = Mean_RAE, None_Mean_ERR = Mean_ERR, None_MEAN_AE = MEAN_AE)

compare_vs_none_if <- summary_by_IF %>%
  filter(Clustering_Method != "None") %>%
  left_join(none_by_if, by = "IF") %>%
  mutate(
    Delta_RAE_vs_None = Mean_RAE - None_Mean_RAE,
    Delta_ERR_vs_None = Mean_ERR - None_Mean_ERR,
    Delta_AE_vs_None = MEAN_AE - None_MEAN_AE
  )

write.csv(compare_vs_none_if, paste0("compare_vs_none_if_no_clus_s1_n_", n_total, ".csv"))
