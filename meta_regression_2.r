library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(metafor)
library(writexl)

# version-safe robust()
robust_safe <- function(model, cluster_vec) {
  rb_formals <- names(formals(metafor::robust))
  if ("small" %in% rb_formals) robust(model, cluster = cluster_vec, small = TRUE) else
    if ("adjust" %in% rb_formals) robust(model, cluster = cluster_vec, adjust = TRUE) else
      robust(model, cluster = cluster_vec)
}

# CI extractor that works across versions/shapes
extract_ci <- function(fit) {
  ct <- suppressWarnings(confint(fit))
  if ((is.matrix(ct) || is.data.frame(ct)) && all(c("ci.lb","ci.ub") %in% colnames(ct))) {
    return(c(lower = as.numeric(ct[1, "ci.lb"]), upper = as.numeric(ct[1, "ci.ub"])))
  }
  if (is.list(ct) && !is.null(ct$beta) && all(c("ci.lb","ci.ub") %in% colnames(ct$beta))) {
    return(c(lower = as.numeric(ct$beta[1, "ci.lb"]), upper = as.numeric(ct$beta[1, "ci.ub"])))
  }
  if (is.numeric(ct) && all(c("ci.lb","ci.ub") %in% names(ct))) {
    return(c(lower = as.numeric(ct[["ci.lb"]]), upper = as.numeric(ct[["ci.ub"]])))
  }
  est <- unname(coef(fit)[1]); se <- sqrt(vcov(fit)[1,1])
  c(lower = est - qnorm(0.975)*se, upper = est + qnorm(0.975)*se)
}

# helper: predict with version safety (works for rma/rma.mv)
# Safer predict that checks if the fitted model has moderators
predict_safe <- function(fit, newmods_mat = NULL) {
  has_mods <- isTRUE(!is.null(fit$X) && NCOL(fit$X) > 1L) || isTRUE(fit$p > 1L)
  if (!has_mods) {
    # Intercept-only: newmods must be NULL
    preds <- predict(fit)
  } else {
    if (is.null(newmods_mat)) stop("newmods must be provided for models with moderators")
    preds <- predict(fit, newmods = newmods_mat)
  }
  c(pred = as.numeric(preds$pred[1]),
    lci  = as.numeric(preds$ci.lb[1]),
    uci  = as.numeric(preds$ci.ub[1]))
}


# core analyzer: intercept-only and meta-regression on study_duration_months
analyze_sheet_srm_meta <- function(sheet_name, path = excel_path, standard_months = 12, scale_months = FALSE) {
  raw <- readxl::read_excel(path, sheet = sheet_name)
  
  req <- c("Number","SRM_value_final","SRM_SE_final","effect_id","study_duration_months")
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0) stop(sprintf("Sheet '%s' missing: %s", sheet_name, paste(miss, collapse=", ")))
  
  df <- raw %>%
    transmute(
      study_id  = as.character(Number),
      effect_id = as.character(effect_id),
      SRM       = as.numeric(SRM_value_final),
      SE        = as.numeric(SRM_SE_final),
      duration  = as.numeric(study_duration_months)
    ) %>%
    mutate(
      yi = SRM,
      vi = SE^2
    ) %>%
    filter(is.finite(yi), is.finite(vi), is.finite(duration))
  
  if (nrow(df) < 2) {
    return(tibble(
      sheet = sheet_name, k = nrow(df),
      n_studies = dplyr::n_distinct(df$study_id),
      pooled_srm = NA_real_, lci = NA_real_, uci = NA_real_,
      tau2_between = NA_real_, sigma2_within = NA_real_,
      I2_studylevel_perc  = NA_real_, I2_withinstudy_perc = NA_real_,
      used_crve = FALSE,
      meta_reg_coef_duration = NA_real_, meta_reg_lci = NA_real_, meta_reg_uci = NA_real_,
      meta_reg_p = NA_real_,
      standardized_months = standard_months,
      adjusted_srm_at_standard_months = NA_real_,
      adjusted_lci_at_standard_months = NA_real_,
      adjusted_uci_at_standard_months = NA_real_
    ))
  }
  
  # Optionally rescale duration (e.g., per year) to aid interpretability/numerics
  if (isTRUE(scale_months)) {
    df <- df %>% mutate(duration_scaled = duration / 12)
    mod_vec <- ~ duration_scaled
    newmods_std <- matrix(standard_months/12, nrow = 1)
  } else {
    df <- df %>% mutate(duration_scaled = duration)  # keep original months for modeling
    mod_vec <- ~ duration_scaled
    newmods_std <- matrix(standard_months, nrow = 1)
  }
  
  # 1) Intercept-only multilevel model (same as your current analysis)
  m0 <- rma.mv(yi, vi, random = ~ 1 | study_id/effect_id, data = df, method = "REML")  # [web:4]
  n_clusters <- dplyr::n_distinct(df$study_id)
  p_fixed0   <- length(coef(m0))
  if (n_clusters > p_fixed0) {
    m0_use <- robust_safe(m0, cluster_vec = df$study_id)  # [web:9][web:12][web:18]
    used_crve0 <- TRUE
  } else {
    m0_use <- m0
    used_crve0 <- FALSE
  }
  est0 <- unname(coef(m0_use)[1])
  ci0  <- extract_ci(m0_use)
  
  # 2) Meta-regression with study duration as moderator
  m1 <- rma.mv(yi, vi, mods = mod_vec, random = ~ 1 | study_id/effect_id, data = df, method = "REML")  # [web:4]
  p_fixed1 <- length(coef(m1))
  if (n_clusters > p_fixed1) {
    m1_use <- robust_safe(m1, cluster_vec = df$study_id)  # [web:9][web:12][web:18]
    used_crve1 <- TRUE
  } else {
    m1_use <- m1
    used_crve1 <- FALSE
  }
  
  # Extract moderator slope and CI for duration
  b <- coef(m1_use)
  vc <- vcov(m1_use)
  b1 <- as.numeric(b[grep("duration_scaled", names(b), fixed = TRUE)])
  se1 <- sqrt(as.numeric(vc[grep("duration_scaled", names(b), fixed = TRUE),
                            grep("duration_scaled", names(b), fixed = TRUE)]))
  zcrit <- qnorm(0.975)
  b1_lci <- b1 - zcrit*se1
  b1_uci <- b1 + zcrit*se1
  
  # 3) Duration-standardized prediction at standard_months
  # For rma.mv with moderators, use predict(newmods=...) to obtain adjusted estimate at chosen covariate value [web:1][web:11][web:17][web:20]
  adj_pred <- predict_safe(m1, newmods_mat = newmods_std)  # [web:1][web:11][web:20]
  V_sampling   <- mean(df$vi)
  # Collate output
  tibble(
    sheet          = sheet_name,
    k              = nrow(df),
    n_studies      = n_clusters,
    pooled_srm     = est0,
    lci            = ci0["lower"],
    uci            = ci0["upper"],
    tau2_between   = m0$sigma2[1],
    sigma2_within  = m0$sigma2[2],
    I2_studylevel_perc  = 100*m0$sigma2[1]/(m0$sigma2[1]+m0$sigma2[2]+V_sampling),
    I2_withinstudy_perc = 100*m0$sigma2[2]/(m0$sigma2[1]+m0$sigma2[2]+V_sampling),
    used_crve      = used_crve0,
    used_crve_meta_reg = used_crve1,
    meta_reg_coef_duration = b1,
    meta_reg_lci   = b1_lci,
    meta_reg_uci   = b1_uci,
    meta_reg_p     = 2*pnorm(-abs(b1/se1)),
    standardized_months = standard_months,
    adjusted_srm_at_standard_months = adj_pred["pred"],
    adjusted_lci_at_standard_months = adj_pred["lci"],
    adjusted_uci_at_standard_months = adj_pred["uci"]
  )
}

# run across sheets
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 10.xlsx"
sheets <- readxl::excel_sheets(excel_path)
summary_tbl_srm_meta <- purrr::map_dfr(sheets, ~analyze_sheet_srm_meta(.x, path = excel_path, standard_months = 12, scale_months = FALSE))

print(summary_tbl_srm_meta)

# Save CSV
out_csv <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table10_metaregression.csv"
write.csv(summary_tbl_srm_meta, out_csv, row.names = FALSE)
