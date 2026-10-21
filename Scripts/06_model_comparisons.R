rm(list=ls())
library(dplyr)
library(lme4)
library(car)
library(performance)
library(stringr)
library(lmerTest)
library(broom.mixed)
library(purrr)
library(emmeans)
library(multcomp)
library(here)
library(sjPlot)
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))

df_modeling <- meta_filteredK %>%
  filter(!is.na(ED50), !is.na(Site), !is.na(Species)) %>%
  mutate(
    YY      = ED50,
    Species = droplevels(as.factor(Species)),
    Site    = droplevels(as.factor(Site)),
    across(starts_with("K_"), ~ droplevels(as.factor(.)))
  )


#test which model is best for ed50 (K2:15)
k_ms_models <- list()
k_range <- 2:15

for (k in k_range) {
  k_col <- paste0("K_", k)
  if (k_col %in% names(df_modeling)) {
    f_kms <- as.formula(paste0("YY ~ Species + ", k_col, " + (1 | Site)"))
    
    k_ms_models[[paste0("k", k, "ms")]] <- tryCatch(
      lmer(f_kms, data = df_modeling, REML = FALSE),
      error = function(e) NULL
    )
  }
}

k_ms_models <- compact(k_ms_models)

# Build Supp Table 1 (Evaluating optimal K)
table_1_kms <- imap_dfr(k_ms_models, function(model_obj, model_name) {
  perf    <- model_performance(model_obj, estimator = "ML")
  r2_vals <- tryCatch(r2(model_obj), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
  k_val   <- str_extract(model_name, "\\d+")
  
  data.frame(
    Model_Name      = model_name,
    K               = k_val,
    Model_Type      = "Species + K (Additive)",
    AIC             = round(perf$AIC, 1),
    BIC             = round(perf$BIC, 1),
    Marginal_R2     = round(r2_vals$R2_marginal, 3),
    Conditional_R2  = round(r2_vals$R2_conditional, 3),
    stringsAsFactors = FALSE
  )
}) %>% arrange(BIC)

cat("\n--- STEP 2: TABLE 1 - ADDITIVE MODEL SELECTION (BEST K) ---\n")
print(head(table_1_kms))
write.csv(table_1_kms, here::here("Outputs", "Table_K_selection.csv"), row.names = FALSE)

evaluate_k_models <- function(metric_name, data = meta_filteredK) {
  
  # Prepare modeling subset ensuring consistent data points
  df_sub <- data %>%
    filter(!is.na(.data[[metric_name]]), !is.na(Site), !is.na(Species)) %>%
    mutate(
      YY      = .data[[metric_name]],
      Species = droplevels(as.factor(Species)),
      Site    = droplevels(as.factor(Site)),
      across(starts_with("K_"), ~ droplevels(as.factor(.)))
    )
  
  k_ms_models <- list()
  k_range <- 2:15
  
  # Base Model (K=0: Taxonomy only) for baseline comparison
  k_ms_models[["K_0"]] <- tryCatch(
    lmer(YY ~ Species + (1 | Site), data = df_sub, REML = FALSE),
    error = function(e) NULL
  )
  
  # Candidate K Models (K=2 to 15)
  for (k in k_range) {
    k_col <- paste0("K_", k)
    if (k_col %in% names(df_sub)) {
      f_kms <- as.formula(paste0("YY ~ Species + ", k_col, " + (1 | Site)"))
      k_ms_models[[paste0("K_", k)]] <- tryCatch(
        lmer(f_kms, data = df_sub, REML = FALSE),
        error = function(e) NULL
      )
    }
  }
  
  k_ms_models <- compact(k_ms_models)
  
  # Extract statistics
  imap_dfr(k_ms_models, function(model_obj, model_name) {
    perf    <- model_performance(model_obj, estimator = "ML")
    r2_vals <- tryCatch(r2(model_obj), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
    
    data.frame(
      K              = model_name,
      AIC            = round(perf$AIC, 1),
      BIC            = round(perf$BIC, 1),
      Marginal_R2    = sprintf("%.3f", r2_vals$R2_marginal),
      Conditional_R2 = sprintf("%.3f", r2_vals$R2_conditional),
      stringsAsFactors = FALSE
    )
  })
}

# -----------------------------------------------------------------------------
# 2. RUN EVALUATIONS FOR ED50 AND COLOR SLOPE
# -----------------------------------------------------------------------------
res_ed50  <- evaluate_k_models("ED50")
res_color <- evaluate_k_models("color_slope") # Adjust column name if needed

write.csv(
  res_ed50, 
  file = here("Outputs", "Supplemental_Table_K_Selection_ED50.csv"), 
  row.names = FALSE
)

write.csv(
  res_color, 
  file = here("Outputs", "Supplemental_Table_K_Selection_color.csv"), 
  row.names = FALSE
)

#Table 1: Species, K_3, or Species+K_3 model evaluations
# Function to run Table 1 Model Selection across any dose-response metric
run_dose_response_selection <- function(metric_col_name, data = meta_filteredK) {
  
  cat("\n=======================================================\n")
  cat(" RUNNING MODEL SELECTION FOR:", metric_col_name, "\n")
  cat("=======================================================\n")
  
  # 1. Prepare data (Filter NAs for target metric and ensure identical N)
  df_sub <- data %>%
    filter(
      !is.na(.data[[metric_col_name]]), 
      !is.na(Site), 
      !is.na(Species),
      !is.na(K_3)
    ) %>%
    mutate(
      YY      = .data[[metric_col_name]],
      Species = droplevels(as.factor(Species)),
      Site    = droplevels(as.factor(Site)),
      K_3     = droplevels(as.factor(K_3))
    )
  
  # 2. Fit Candidate Models (REML = FALSE for ML comparison)
  m_spp  <- lmer(YY ~ Species + (1 | Site), data = df_sub, REML = FALSE)
  m_k3   <- lmer(YY ~ K_3 + (1 | Site), data = df_sub, REML = FALSE)
  m_full <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_sub, REML = FALSE)
  
  model_list <- list(
    "Species Only"         = m_spp,
    "Cluster (K=3) Only"   = m_k3,
    "Species + K_3 (Full)" = m_full
  )
  
  # 3. Build Model Selection Table
  table_out <- imap_dfr(model_list, function(mod, name) {
    perf    <- model_performance(mod, estimator = "ML")
    r2_vals <- tryCatch(r2(mod), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
    
    # Variance components
    var_comp  <- as.data.frame(VarCorr(mod))
    site_var  <- var_comp$vcov[var_comp$grp == "Site"]
    resid_var <- var_comp$vcov[var_comp$grp == "Residual"]
    tot_var   <- site_var + resid_var
    
    anova_res <- as.data.frame(anova(mod, type = 2))
    
    # Species F-test
    f_spp <- if ("Species" %in% rownames(anova_res)) {
      sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
              anova_res["Species", "NumDF"], anova_res["Species", "DenDF"], 
              anova_res["Species", "F value"], anova_res["Species", "Pr(>F)"])
    } else "NA"
    
    # K3 F-test
    f_k3 <- if ("K_3" %in% rownames(anova_res)) {
      sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
              anova_res["K_3", "NumDF"], anova_res["K_3", "DenDF"], 
              anova_res["K_3", "F value"], anova_res["K_3", "Pr(>F)"])
    } else "NA"
    
    data.frame(
      Metric            = metric_col_name,
      Model             = name,
      AIC               = round(perf$AIC, 1),
      BIC               = round(perf$BIC, 1),
      Marginal_R2       = sprintf("%.3f (%.1f%%)", r2_vals$R2_marginal, r2_vals$R2_marginal * 100),
      Conditional_R2    = sprintf("%.3f (%.1f%%)", r2_vals$R2_conditional, r2_vals$R2_conditional * 100),
      Site_Var_Pct      = sprintf("%.1f%%", (site_var / tot_var) * 100),
      Species_Effect    = f_spp,
      K3_Cluster_Effect = f_k3,
      stringsAsFactors  = FALSE
    )
  })
  
  return(table_out)
}

# Execute model selection across the complete physiological and optical suite
tbl_ed5    <- run_dose_response_selection("ED5")
tbl_ed50   <- run_dose_response_selection("ED50")
tbl_ed95   <- run_dose_response_selection("ED95")
tbl_dw     <- run_dose_response_selection("DW")
tbl_slope  <- run_dose_response_selection("color_slope") # or "slope" / "Paling"

# Combine into a single comprehensive output table
master_trait_table <- bind_rows(
  tbl_ed5, 
  tbl_ed50, 
  tbl_ed95, 
  tbl_dw, 
  tbl_slope
)

# Export or inspect table
print(master_trait_table)
write.csv(master_trait_table, here::here("Outputs", "Table_S1_full_response_model_table.csv"), row.names = FALSE)
##################
#here want to see anova 
library(dplyr)
library(lme4)
library(lmerTest) # Loading lmerTest enables Satterthwaite df for stats::anova()
library(performance)
library(purrr)

run_dose_response_selection <- function(data, metric_col_name) {
  
  # 1. Prepare data (Filter NAs for target metric and ensure identical N)
  df_sub <- data %>%
    filter(
      !is.na(.data[[metric_col_name]]), 
      !is.na(Site), 
      !is.na(Species),
      !is.na(K_3)
    ) %>%
    mutate(
      YY      = .data[[metric_col_name]],
      Species = droplevels(as.factor(Species)),
      Site    = droplevels(as.factor(Site)),
      K_3     = droplevels(as.factor(K_3))
    )
  
  # 2. Fit Candidate Models for Selection (ML / REML = FALSE)
  m_spp_ml    <- lmer(YY ~ Species + (1 | Site), data = df_sub, REML = FALSE)
  m_k3_ml     <- lmer(YY ~ K_3 + (1 | Site), data = df_sub, REML = FALSE)
  m_full_ml   <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_sub, REML = FALSE)
  m_k3sppr_ml <- lmer(YY ~ K_3 + (1 | Species) + (1 | Site), data = df_sub, REML = FALSE)
  
  # 3. Fit Identical Models for ANOVA Inference (REML = TRUE)
  m_spp_reml    <- lmer(YY ~ Species + (1 | Site), data = df_sub, REML = TRUE)
  m_k3_reml     <- lmer(YY ~ K_3 + (1 | Site), data = df_sub, REML = TRUE)
  m_full_reml   <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_sub, REML = TRUE)
  m_k3sppr_reml <- lmer(YY ~ K_3 + (1 | Species) + (1 | Site), data = df_sub, REML = TRUE)
  
  model_pairs <- list(
    "Species Only"               = list(ml = m_spp_ml,    reml = m_spp_reml),
    "Cluster (K=3) Only"         = list(ml = m_k3_ml,     reml = m_k3_reml),
    "Species + K_3 (Full)"       = list(ml = m_full_ml,   reml = m_full_reml),
    "K_3 + (1|Species) (Random)" = list(ml = m_k3sppr_ml, reml = m_k3sppr_reml)
  )
  
  # 4. Extract Performance Metrics & ANOVA F-Statistics
  table_out <- imap_dfr(model_pairs, function(mod_pair, name) {
    
    # Selection metrics from ML model
    perf    <- model_performance(mod_pair$ml, estimator = "ML")
    r2_vals <- tryCatch(r2(mod_pair$reml), error = function(e) list(R2_marginal = NA, R2_conditional = NA))
    
    # Variance components from REML model
    var_comp  <- as.data.frame(VarCorr(mod_pair$reml))
    site_var  <- ifelse("Site" %in% var_comp$grp, var_comp$vcov[var_comp$grp == "Site"], 0)
    resid_var <- var_comp$vcov[var_comp$grp == "Residual"]
    tot_var   <- sum(var_comp$vcov)
    
    # Type II Satterthwaite ANOVA using stats::anova with lmerTest active
    anova_res <- as.data.frame(anova(mod_pair$reml, type = 2))
    
    # Helper to format F-statistic strings
    format_f <- function(term) {
      if (term %in% rownames(anova_res)) {
        sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
                anova_res[term, "NumDF"], anova_res[term, "DenDF"], 
                anova_res[term, "F value"], anova_res[term, "Pr(>F)"])
      } else "—"
    }
    
    data.frame(
      Metric            = metric_col_name,
      Model             = name,
      AIC               = round(perf$AIC, 1),
      BIC               = round(perf$BIC, 1),
      Marginal_R2       = sprintf("%.3f (%.1f%%)", r2_vals$R2_marginal, r2_vals$R2_marginal * 100),
      Conditional_R2    = sprintf("%.3f (%.1f%%)", r2_vals$R2_conditional, r2_vals$R2_conditional * 100),
      Site_Var_Pct      = sprintf("%.1f%%", (site_var / tot_var) * 100),
      Species_Effect    = format_f("Species"),
      K3_Cluster_Effect = format_f("K_3"),
      stringsAsFactors  = FALSE
    )
  })
  
  return(table_out)
}

# Run selection across traits
metrics_to_test <- c("ED5", "ED50", "ED95", "DW", "color_slope")

master_anova_selection <- map_dfr(metrics_to_test, ~run_dose_response_selection(data = df_modeling, metric_col_name = .x))

# Save master selection table
write.csv(master_anova_selection, here::here("Outputs", "Table_DoseResponse_Model_Selection.csv"), row.names = FALSE)

# Inspect output
print(master_anova_selection, width = Inf)

# ==============================================================================
# 1. FIT FULL ADDITIVE REML MODELS FOR ED50 AND ED95
# ==============================================================================
m_ed50_full <- lmer(ED50 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)
m_ed95_full <- lmer(ED95 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)

# ==============================================================================
# 2. EXTRACT POST-HOC TUKEY HSD PAIRWISE COMPARISONS
# ==============================================================================
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(purrr)
library(here)

# 1. Clean data and factors
df_modeling <- df_modeling %>%
  filter(!is.na(Site), !is.na(Species), !is.na(K_3)) %>%
  mutate(
    Species = droplevels(as.factor(Species)),
    K_3     = droplevels(as.factor(K_3)),
    Site    = droplevels(as.factor(Site))
  )

# 2. Helper function to extract and format Tukey pairwise post-hoc contrasts
get_sig_posthoc <- function(model, metric_name, factor_name) {
  
  # Calculate EMMs and pairwise Tukey contrasts
  emm <- emmeans(model, specs = as.formula(paste0("~ ", factor_name)))
  pw  <- contrast(emm, method = "pairwise", adjust = "tukey")
  
  # Format output table
  as.data.frame(pw) %>%
    transmute(
      Metric     = metric_name,
      Factor     = factor_name,
      Contrast   = contrast,
      Estimate   = round(estimate, 3),
      SE         = round(SE, 3),
      df         = round(df, 1),
      t_ratio    = round(t.ratio, 2),
      p_value    = ifelse(p.value < 0.001, "< 0.001", sprintf("%.4f", p.value)),
      Signif     = case_when(
        p.value < 0.001 ~ "***",
        p.value < 0.01  ~ "**",
        p.value < 0.05  ~ "*",
        TRUE            ~ "ns"
      )
    )
}

# 3. Fit REML models and extract significant post-hocs

# --- ED50 (Both Species & K3 Significant) ---
m_ed50 <- lmer(ED50 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)
ph_ed50_spp <- get_sig_posthoc(m_ed50, "ED50", "Species")
ph_ed50_k3  <- get_sig_posthoc(m_ed50, "ED50", "K_3")

# --- ED95 (Both Species & K3 Significant) ---
m_ed95 <- lmer(ED95 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)
ph_ed95_spp <- get_sig_posthoc(m_ed95, "ED95", "Species")
ph_ed95_k3  <- get_sig_posthoc(m_ed95, "ED95", "K_3")

# --- ED5 (K3 Only Significant) ---
m_ed5  <- lmer(ED5 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)
ph_ed5_k3   <- get_sig_posthoc(m_ed5, "ED5", "K_3")

# --- Dynamic Width (DW) (K3 Only Significant) ---
m_dw   <- lmer(DW ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)
ph_dw_k3    <- get_sig_posthoc(m_dw, "DW", "K_3")

# 4. Combine into Master Post-Hoc Table
master_posthoc_table <- bind_rows(
  ph_ed50_spp,
  ph_ed50_k3,
  ph_ed95_spp,
  ph_ed95_k3,
  ph_ed5_k3,
  ph_dw_k3
)

cat("\n--- POST-HOC TUKEY HSD PAIRWISE COMPARISONS (SIGNIFICANT FULL-MODEL TERMS) ---\n")
print(master_posthoc_table, row.names = FALSE)

# 5. Export for Supplementary Data
write.csv(
  master_posthoc_table, 
  here("Outputs", "Table_S3_PostHoc_Pairwise_Summary.csv"), 
  row.names = FALSE
)
library(partR2)
library(lme4)

# Fit REML full model
m_ed50_full <- lmer(ED50 ~ Species + K_3 + (1 | Site), data = df_modeling, REML = TRUE)

# Calculate semi-partial R2 for each fixed term
part_r2_ed50 <- partR2(
  m_ed50_full, 
  partvars = c("Species", "K_3"), 
  R2_type = "marginal",
  nboot = 100 # Set higher (e.g. 1000) for manuscript CIs
)

summary(part_r2_ed50)
###########################
library(dplyr)
library(lme4)
library(lmerTest)
library(performance)
library(purrr)

get_winning_model_outputs <- function(data, metric_col_name) {
  
  # 1. Subset data for target metric
  df_sub <- data %>%
    filter(
      !is.na(.data[[metric_col_name]]), 
      !is.na(Site), 
      !is.na(Species),
      !is.na(K_3)
    ) %>%
    mutate(
      YY      = .data[[metric_col_name]],
      Species = droplevels(as.factor(Species)),
      Site    = droplevels(as.factor(Site)),
      K_3     = droplevels(as.factor(K_3))
    )
  
  # 2. Fit Final Full Additive Model (REML = TRUE)
  m_full <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_sub, REML = TRUE)
  
  # 3. Extract R2 Metrics safely
  r2_vals <- tryCatch(
    suppressWarnings(r2(m_full)), 
    error = function(e) list(R2_marginal = NA, R2_conditional = NA)
  )
  r2_m <- if (!is.null(r2_vals$R2_marginal)) r2_vals$R2_marginal else NA
  r2_c <- if (!is.null(r2_vals$R2_conditional)) r2_vals$R2_conditional else r2_m
  
  # 4. Extract Site Variance (%)
  var_comp  <- as.data.frame(VarCorr(m_full))
  site_var  <- ifelse("Site" %in% var_comp$grp, var_comp$vcov[var_comp$grp == "Site"], 0)
  tot_var   <- sum(var_comp$vcov)
  site_pct  <- (site_var / tot_var) * 100
  
  # 5. Extract Type II Satterthwaite ANOVA Results
  anova_res <- as.data.frame(anova(m_full, type = 2))
  
  format_f <- function(term) {
    if (term %in% rownames(anova_res)) {
      sprintf("F(%.0f, %.1f) = %.2f, p = %.4f", 
              anova_res[term, "NumDF"], anova_res[term, "DenDF"], 
              anova_res[term, "F value"], anova_res[term, "Pr(>F)"])
    } else "—"
  }
  
  # 6. Format Final Row
  data.frame(
    Metric            = metric_col_name,
    Model_Formula     = "Species + K_3 + (1|Site)",
    Marginal_R2       = ifelse(is.na(r2_m), "NA", sprintf("%.3f (%.1f%%)", r2_m, r2_m * 100)),
    Conditional_R2    = ifelse(is.na(r2_c), "NA", sprintf("%.3f (%.1f%%)", r2_c, r2_c * 100)),
    Site_Var_Pct      = sprintf("%.1f%%", site_pct),
    Species_Effect    = format_f("Species"),
    K3_Cluster_Effect = format_f("K_3"),
    stringsAsFactors  = FALSE
  )
}

# ==============================================================================
# EXECUTE ACROSS TRAIT SUITE
# ==============================================================================
metrics_to_run <- c("ED5", "ED50", "ED95", "DW", "color_slope")

table_winning_models <- map_dfr(
  metrics_to_run, 
  ~get_winning_model_outputs(data = df_modeling, metric_col_name = .x)
)

cat("\n--- SUMMARY OF FULL ADDITIVE MODELS (SPECIES + K_3 + (1|SITE)) ---\n")
print(table_winning_models, width = Inf)

# Save to CSV
write.csv(
  table_winning_models, 
  here::here("Outputs", "Table_FullModel_Performance_Effects.csv"), 
  row.names = FALSE
)


##########################
#only looking here at ed50~spp+K_3+(1|Site)
k3ms   <- lmer(YY ~ K_3 + Species+ (1 | Site), data = df_modeling, REML = TRUE)
# coefs
coefs <- broom.mixed::tidy(k3ms, effects = "fixed", conf.int = TRUE) %>%
  transmute(
    Term          = term,
    Estimate      = round(estimate, 3),
    `Std. Error`  = round(std.error, 3),
    `2.5% CI`     = round(conf.low, 3),
    `97.5% CI`    = round(conf.high, 3),
    `Statistic`   = round(statistic, 2),
    `p-value`     = ifelse(p.value < 0.001, "< 0.001", sprintf("%.4f", p.value))
  );coefs

# random
random_vars <- broom.mixed::tidy(k3ms, effects = "ran_pars") %>%
  transmute(
    Group    = group,
    Term     = term,
    StdDev   = round(estimate, 3),
    Variance = round(estimate^2, 3)
  );random_vars

plot_model(k3ms, type = "re")

# ==============================================================================
# Diagnostics
# ==============================================================================
# --- Residual Diagnostics ---
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(resid(k3ms) ~ fitted(k3ms), main = "Residuals vs Fitted (k3msr)", xlab = "Fitted", ylab = "Residuals")
qqnorm(resid(k3ms), main = "Normal Q-Q Plot (k3ms)")
qqline(resid(k3ms), col = "red")
par(mfrow = c(1, 1))

# ==============================================================================
# 2. FOR SPP: Extract EMMs & Compact Letter Display (CLD)
# ==============================================================================
k3_emmeans_spp <- emmeans(k3ms, pairwise ~ Species, adjust = "tukey")
k3_emm_df       <- as.data.frame(k3_emmeans_spp$emmeans)
k3_contrasts_df <- as.data.frame(k3_emmeans_spp$contrasts)

cat("\n--- Estimated Marginal Means (ED50) for K_3 Clusters ---\n")
print(k3_emm_df %>% dplyr::select(Species, emmean, SE, df, lower.CL, upper.CL))

cat("\n--- Pairwise Tukey Contrasts Between K_3 Clusters ---\n")
print(k3_contrasts_df %>% dplyr::select(contrast, estimate, SE, p.value))

# Compute Tukey CLD Letters
cld_spp   <- cld(k3_emmeans_spp$emmeans, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group))

# # ==============================================================================
# # 3. Align Cluster Labels & Factor Ordering
# # ==============================================================================
# Verify counts in spp
cat("Raw Data Counts per spp Label:\n")
print(table(df_modeling$Species, useNA = "always"))

# ==============================================================================
# 4. Publication-Ready Plot
# ==============================================================================
y_min_val <- min(df_modeling$ED50, na.rm = TRUE) - 0.5
y_max_val <- max(cld_spp$upper.CL, na.rm = TRUE) + 0.8

plot_k3spp <- ggplot() +
  # Layer A: Raw data points
  # geom_jitter(
  #   data = df_modeling, 
  #   aes(x = Species, y = ED50, color = Species),
  #   width = 0.12, 
  #   alpha = 0.5, 
  #   size = 2
  # ) +
  # Layer B: Model Estimated Marginal Means & 95% CIs
  geom_errorbar(
    data = cld_spp,
    aes(x = Species, ymin = lower.CL, ymax = upper.CL),
    width = 0.15, 
    linewidth = 0.8, 
    color = "black"
  ) +
  geom_point(
    data = cld_spp,
    aes(x = Species, y = emmean, fill = Species, shape = Species),
    size = 5,
    #fill = Species,
    color = "black",
    #shape = 21
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_spp,
    aes(x = Species, y = upper.CL + 0.5, label = .group),
    #fontface = "bold", 
    hjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Custom Aesthetics & Formatting
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    y = expression("ED50" ~ "(°C)"),
    color = "Species"
  ) +
  
  # FORCES X-AXIS ORDER & ITALIC SPECIES LABELS
  scale_x_discrete(
    limits = species_order, 
    labels = species_labels
  ) +
  
  scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_fill_manual(values = species_colors, limits = species_order, labels = species_labels) +
  scale_shape_manual(values = species_shapes, limits = species_order, labels = species_labels)+
  
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title      = element_text(face = "bold", size = 13, color = "black"),
    axis.text       = element_text(size = 11,  color = "black"),
    axis.text.x     = element_text(face = "italic", vjust = 1, hjust = 1, angle = 30,size = 11, color = "black"), # Added italic for species names
    axis.line       = element_line(linewidth = 0.6, color = "black"),
    axis.ticks      = element_line(linewidth = 0.6, color = "black"),
    legend.position = "none",
    # legend.title    = element_text(face = "bold", size = 11),
    # legend.text     = element_text(size = 11, face = "italic"),
    plot.margin     = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  );plot_k3spp

saveRDS(plot_k3spp, "Outputs/plot_ed50K_spp.rds")

# ==============================================================================
# 3. For Cluster: Extract EMMs & Compact Letter Display (CLD)
# ==============================================================================
k3_emmeans_obj <- emmeans(k3ms, pairwise ~ K_3, adjust = "tukey")
k3_emm_df       <- as.data.frame(k3_emmeans_obj$emmeans)
k3_contrasts_df <- as.data.frame(k3_emmeans_obj$contrasts)

cat("\n--- Estimated Marginal Means (ED50) for K_3 Clusters ---\n")
print(k3_emm_df %>% dplyr::select(K_3, emmean, SE, df, lower.CL, upper.CL))

cat("\n--- Pairwise Tukey Contrasts Between K_3 Clusters ---\n")
print(k3_contrasts_df %>% dplyr::select(contrast, estimate, SE, p.value))

# Compute Tukey CLD Letters
cld_clusters   <- cld(k3_emmeans_obj$emmeans, adjust = "tukey", Letters = letters) %>%
  as.data.frame() %>%
  mutate(.group = trimws(.group))

# # ==============================================================================
# # 3. Align Cluster Labels & Factor Ordering
# # ==============================================================================
df_modeling <- df_modeling %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

cld_clusters <- cld_clusters %>%
  mutate(
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

# Verify counts in console
cat("Raw Data Counts per Cluster Label:\n")
print(table(df_modeling$Cluster_Label, useNA = "always"))

cat("\nEMM Summary Counts per Cluster Label:\n")
print(table(cld_clusters$Cluster_Label, useNA = "always"))
# ==============================================================================
# 4. Publication-Ready Plot
# ==============================================================================
y_min_val <- min(df_modeling$ED50, na.rm = TRUE) - 0.5
y_max_val <- max(cld_clusters$upper.CL, na.rm = TRUE) + 0.8

plot_k3msr <- ggplot() +
  # Layer A: Raw data points
  # geom_jitter(
  #   data = df_modeling, 
  #   aes(x = Cluster_Label, y = ED50, color = Species),
  #   width = 0.12, 
  #   alpha = 0.8, 
  #   size = 3
  # ) +
  # Layer B: Model Estimated Marginal Means & 95% CIs
  geom_errorbar(
    data = cld_clusters,
    aes(x = Cluster_Label, ymin = lower.CL, ymax = upper.CL),
    width = 0.15, 
    linewidth = 0.8, 
    color = "black"
  ) +
  geom_point(
    data = cld_clusters,
    aes(x = Cluster_Label, y = emmean, fill = Cluster_Label),
    size = 5,
    shape=21
    # fill = "black",
    # color = "black"
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_clusters,
    aes(x = Cluster_Label, y = upper.CL + 0.5, label = .group),
    #fontface = "bold", 
    hjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Custom Aesthetics & Formatting
  #scale_color_manual(values = species_colors, limits = species_order, labels = species_labels) +
 # scale_fill_manual(values = species_colors, limits = species_order, labels = species_labels) +
 # scale_color_manual(values = cluster_colors, limits = roman_map) +
  #scale_fill_manual(values = cluster_colors, limits = roman_map) +
  # Discrete Fill Scale using target_levels for limits
  scale_color_manual(name = "Cluster", values = cluster_colors, limits = target_levels) +
  scale_fill_manual(name = "Cluster", values = cluster_colors, limits = target_levels) +
  #scale_shape_manual(values = custom_shapes)+
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    x = "Photophysiological Cluster (K = 3)",
    y = expression("ED50" ~ "(°C)"),
    color = "Cluster",
    #shape = "Site"
  ) +
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III", "Cluster IV" = "IV")
  ) +
  #ggrepel::geom_text_repel()+
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title      = element_text(face = "bold", size = 13, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    axis.text.x     = element_text(face = "bold", size = 11, color = "black"),
    axis.line       = element_line(linewidth = 0.6, color = "black"),
    axis.ticks      = element_line(linewidth = 0.6, color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11),
    plot.margin = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  );plot_k3msr

saveRDS(plot_k3msr, "Outputs/plot_ed50K_3.rds")

plot_k3msr <- ggplot() +
  # Layer B: Model Estimated Marginal Means & 95% CIs
  geom_errorbar(
    data = cld_clusters,
    aes(x = Cluster_Label, ymin = lower.CL, ymax = upper.CL),
    width = 0.15, 
    linewidth = 0.8, 
    color = "black"
  ) +
  geom_point(
    data = cld_clusters,
    aes(x = Cluster_Label, y = emmean, fill = Cluster_Label, shape = species_shapes),
    size = 5,
    #shape = 21,       # Required so points accept 'fill' colors with a black outline
    color = "black"
  ) +
  # Layer C: Significance Letters (CLD)
  geom_text(
    data = cld_clusters,
    aes(x = Cluster_Label, y = upper.CL + 0.7, label = .group),
    fontface = "bold", 
    hjust = 0.5,
    size = 6,
    color = "black",
    inherit.aes = FALSE
  ) +
  # Enforce Cluster Order on X-Axis
  scale_x_discrete(limits = target_levels) +
  
  # Discrete Fill Scale using target_levels for limits
  scale_fill_manual(
    name = "Cluster",
    values = cluster_colors, 
    limits = target_levels
  ) +
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")
  ) +
  
  coord_cartesian(ylim = c(y_min_val, y_max_val)) +
  labs(
    x = "Photophysiological Cluster (K = 3)",
    y = expression("ED50" ~ "(°C)")
  ) +
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    axis.title      = element_text(face = "bold", size = 13, color = "black"),
    axis.text       = element_text(size = 11, color = "black"),
    axis.text.x     = element_text( size = 11, color = "black"), # Non-italic for clusters
    #axis.title.y     = element_text(face = "bold", size = 11, color = "black"),
    axis.line       = element_line(linewidth = 0.6, color = "black"),
    axis.ticks      = element_line(linewidth = 0.6, color = "black"),
    legend.position = "right",
    legend.title    = element_text(face = "bold", size = 11),
    legend.text     = element_text(size = 11),
    plot.margin     = margin(t = 10, r = 10, b = 5, l = 10, unit = "pt")
  )

print(plot_k3msr)
saveRDS(plot_k3msr, "Outputs/plot_ed50K_3.rds")

# ==============================================================================
# FACETED DOSE-RESPONSE TRAIT PLOTTING SCRIPT
# ==============================================================================
library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(multcomp)
library(ggplot2)
library(purrr)
library(tidyr)
library(here)

# Load helper variables (species_colors, cluster_colors, target_levels, roman_map, etc.)
source(here("Scripts/00_visualization_prep.R"))
load(here("Outputs", "photophys_preprocessed_data.RData"))

# 1. Clean data & ensure factor definitions
df_modeling <- meta_filteredK %>%
  filter(!is.na(Site), !is.na(Species), !is.na(K_3)) %>%
  mutate(
    Species       = droplevels(as.factor(Species)),
    K_3           = droplevels(as.factor(K_3)),
    Site          = droplevels(as.factor(Site)),
    Cluster_Label = roman_map[as.character(K_3)],
    Cluster_Label = factor(Cluster_Label, levels = target_levels)
  )

# Define dose-response metrics to facet
dr_metrics <- c("ED5", "ED50", "ED95", "DW")

# Clean metric display labels for facets
metric_facet_labels <- c(
  "ED5"  = "ED[5]~('°C')",
  "ED50" = "ED[50]~('°C')",
  "ED95" = "ED[95]~('°C')",
  "DW"   = "Decline~Width~('°C')"
)

# ==============================================================================
# 2. ITERATE OVER METRICS: FIT MODELS, EXTRACT EMMEANS & CLDs
# ==============================================================================
all_cld_list <- list()

for (m in dr_metrics) {
  
  # Prepare temporary dataset for target metric
  df_temp <- df_modeling %>%
    filter(!is.na(.data[[m]])) %>%
    mutate(YY = .data[[m]])
  
  # Fit winning additive REML model
  m_reml <- lmer(YY ~ Species + K_3 + (1 | Site), data = df_temp, REML = TRUE)
  
  # Extract EMMs and calculate Tukey CLD letters across Cluster
  emm_obj <- emmeans(m_reml, pairwise ~ K_3, adjust = "tukey")
  
  cld_df <- cld(emm_obj$emmeans, adjust = "tukey", Letters = letters) %>%
    as.data.frame() %>%
    mutate(
      Metric        = m,
      Cluster_Label = roman_map[as.character(K_3)],
      Cluster_Label = factor(Cluster_Label, levels = target_levels),
      .group        = trimws(.group)
    )
  
  all_cld_list[[m]] <- cld_df
}

# Master data frame containing EMMs, CIs, and CLDs for all facets
master_cld_df <- bind_rows(all_cld_list) %>%
  mutate(Metric = factor(Metric, levels = dr_metrics))

# ==============================================================================
# 3. BUILD FACETED PUBLICATION-READY FIGURE (Clusters I, II, III)
# ==============================================================================
plot_faceted_dr <- ggplot(data = master_cld_df, aes(x = Cluster_Label, y = emmean)) +
  
  # Layer A: 95% Confidence Intervals
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.18,
    linewidth = 0.7,
    color = "black"
  ) +
  
  # Layer B: Estimated Marginal Means (Points filled by Cluster)
  geom_point(
    aes(fill = Cluster_Label),
    size = 4,
    shape = 21,
    color = "black"
  ) +
  
  # Layer C: Tukey Significance Letters
  geom_text(
    aes(y = upper.CL + (max(upper.CL) - min(lower.CL)) * 0.08, label = .group),
    fontface = "bold",
    size = 4.5,
    color = "black",
    vjust = 0
  ) +
  
  # Facet Panel Layout (2x2 grid with free y-axes)
  facet_wrap(
    ~ Metric, 
    scales = "free_y", 
    ncol = 2,
    labeller = as_labeller(metric_facet_labels, default = label_parsed)
  ) +
  
  # Manual Aesthetics & Scales
  scale_fill_manual(
    name   = "Cluster",
    values = cluster_colors,
    limits = target_levels
  ) +
  scale_x_discrete(
    limits = target_levels,
    labels = c("Cluster I" = "I", "Cluster II" = "II", "Cluster III" = "III")
  ) +
  
  # Axes & Labels
  labs(
    x = "Photophysiological Cluster (K = 3)",
    y = "Estimated Marginal Mean (°C)"
  ) +
  
  # Publication-Ready Theme
  theme_classic(base_size = 12, base_family = "sans") +
  theme(
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(face = "bold", size = 12, color = "black"),
    axis.title       = element_text(face = "bold", size = 13, color = "black"),
    axis.text        = element_text(size = 11, color = "black"),
    axis.text.x      = element_text(face = "bold", size = 11, color = "black"),
    axis.line        = element_line(linewidth = 0.6, color = "black"),
    axis.ticks       = element_line(linewidth = 0.6, color = "black"),
    panel.spacing    = unit(1.2, "lines"),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 11)
  )

# Display plot
print(plot_faceted_dr)

# Save high-resolution outputs
ggsave(
  filename = here("Outputs", "Figure_Faceted_DoseResponse_Clusters.png"),
  plot     = plot_faceted_dr,
  width    = 8.5,
  height   = 7.0,
  dpi      = 300
)

saveRDS(plot_faceted_dr, here("Outputs", "plot_faceted_dr.rds"))
