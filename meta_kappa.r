# packages (same as your CV template)
install.packages(c("readxl", "dplyr", "purrr", "stringr", "tibble", "metafor", "writexl"))
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(metafor)
library(writexl)



# ---- Utilities (reused) ----

robust_safe <- function(model, cluster_vec) {
  rb_formals <- names(formals(metafor::robust))
  if ("small" %in% rb_formals) robust(model, cluster = cluster_vec, small = TRUE)
  else if ("adjust" %in% rb_formals) robust(model, cluster = cluster_vec, adjust = TRUE)
  else robust(model, cluster = cluster_vec)
}

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

# Fisher r-to-z transform and back-transform
fisher_z <- function(r) 0.5 * log((1 + r) / (1 - r))
inv_fisher_z <- function(z) (exp(2*z) - 1) / (exp(2*z) + 1)



# ---- Kappa analyzer ----
# Expected columns (customize if needed):
# Number (study id), effect_id, Kappa_value_final (or 'estimate'), Kappa_SE_final (SE on kappa scale)
# Note: Cohen's kappa is correlation-like in [-1,1]; Fisher r-to-z used here for pooling.
analyze_sheet_kappa <- function(sheet_name, path = excel_path,
                                est_col = "Kappa_value_final",
                                se_col  = "Kappa_SE_final") {
  raw <- readxl::read_excel(path, sheet = sheet_name)
  
  req <- c("Number", est_col, se_col, "effect_id")
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0) stop(sprintf("Sheet '%s' missing: %s", sheet_name, paste(miss, collapse=", ")))
  
  df <- raw %>%
    transmute(
      study_id  = as.character(Number),
      effect_id = as.character(effect_id),
      estimate  = as.numeric(.data[[est_col]]),
      SE        = as.numeric(.data[[se_col]])
    ) %>%
    mutate(
      yi = fisher_z(estimate),
      vi = SE^2
    ) %>%
    filter(is.finite(yi), is.finite(vi))
  
  if (nrow(df) < 2) {
    return(tibble(sheet = sheet_name, k = nrow(df),
                  n_studies = dplyr::n_distinct(df$study_id),
                  pooled     = NA_real_, lci = NA_real_, uci = NA_real_,
                  tau2_between = NA_real_, sigma2_within = NA_real_,
                  used_crve = FALSE,
                  metric = "Kappa"))
  }
  
  m <- rma.mv(yi, vi, random = ~ 1 | study_id/effect_id, data = df, method = "REML")
  
  n_clusters <- dplyr::n_distinct(df$study_id)
  p_fixed    <- length(coef(m))
  
  if (n_clusters > p_fixed) {
    m_use <- robust_safe(m, cluster_vec = df$study_id)
    used_crve <- TRUE
  } else {
    m_use <- m
    used_crve <- FALSE
  }
  
  est_z <- unname(coef(m_use)[1])
  ci_z  <- extract_ci(m_use)
  
  est_r <- inv_fisher_z(est_z)
  lci_r <- inv_fisher_z(ci_z["lower"])
  uci_r <- inv_fisher_z(ci_z["upper"])
  
  tibble(
    sheet          = sheet_name,
    k              = nrow(df),
    n_studies      = n_clusters,
    pooled         = est_r,
    lci            = lci_r,
    uci            = uci_r,
    tau2_between   = m$sigma2[1],
    sigma2_within  = m$sigma2[2],
    I2_studylevel_perc  = 100*m$sigma2[1]/(m$sigma2[1]+m$sigma2[2]+1),
    I2_withinstudy_perc = 100*m$sigma2[2]/(m$sigma2[1]+m$sigma2[2]+1),
    used_crve      = used_crve,
    metric         = "Kappa"
  )
}

# ---- Run over sheets ----
# If you have separate sets of sheets for ICC and Kappa, subset 'sheets' accordingly.
# table 5
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 5.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
kappa_summary <- purrr::map_dfr(sheets, ~analyze_sheet_kappa(.x))


print(kappa_summary)

# Optional: write results
write.csv(kappa_summary, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_table5.csv", row.names = FALSE)


# table 6
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 6.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
kappa_summary <- purrr::map_dfr(sheets, ~analyze_sheet_kappa(.x))
print(kappa_summary)
# Optional: write results
write.csv(kappa_summary, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_table6.csv", row.names = FALSE)

# FNIH
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/FNIH Kappa.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
kappa_summary <- purrr::map_dfr(sheets, ~analyze_sheet_kappa(.x))
print(kappa_summary)
# Optional: write results
write.csv(kappa_summary, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_FNIHKappa.csv", row.names = FALSE)

# FNIH
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/treatment Kappa.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
kappa_summary <- purrr::map_dfr(sheets, ~analyze_sheet_kappa(.x))
print(kappa_summary)
# Optional: write results
write.csv(kappa_summary, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_treatmentKappa.csv", row.names = FALSE)
