install.packages(c("readxl", "dplyr", "purrr", "stringr", "tibble", "metafor", "writexl"))
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(metafor)
library(writexl)
#excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_July18/Table2.xlsx"   # e.g., "/path/to/Table2.xlsx"
#sheets <- readxl::excel_sheets(excel_path)
#sheets



# version-safe robust()
robust_safe <- function(model, cluster_vec) {
  rb_formals <- names(formals(metafor::robust))
  if ("small" %in% rb_formals) robust(model, cluster = cluster_vec, small = TRUE)
  else if ("adjust" %in% rb_formals) robust(model, cluster = cluster_vec, adjust = TRUE)
  else robust(model, cluster = cluster_vec)
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
  # Wald fallback
  est <- unname(coef(fit)[1]); se <- sqrt(vcov(fit)[1,1])
  c(lower = est - qnorm(0.975)*se, upper = est + qnorm(0.975)*se)
}

analyze_sheet <- function(sheet_name, path = excel_path) {
  raw <- readxl::read_excel(path, sheet = sheet_name)
  
  req <- c("Number","CV_value_final","CV_SE_final","effect_id")
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0) stop(sprintf("Sheet '%s' missing: %s", sheet_name, paste(miss, collapse=", ")))
  
  df <- raw %>%
    transmute(
      study_id  = as.character(Number),
      effect_id = as.character(effect_id),
      CV        = as.numeric(CV_value_final),   # 5 means 5%
      SE        = as.numeric(CV_SE_final)       # 0.4 means 0.4%
    ) %>%
    mutate(
      CV_num = CV/100,
      SE_num = SE/100,
      yi     = log(CV_num),
      vi     = SE_num^2
    ) %>%
    filter(is.finite(yi), is.finite(vi))
  
  if (nrow(df) < 2) {
    return(tibble(sheet = sheet_name, k = nrow(df),
                  n_studies = dplyr::n_distinct(df$study_id),
                  pooled_cv = NA_real_, lci = NA_real_, uci = NA_real_,
                  tau2_between = NA_real_, sigma2_within = NA_real_,
                  used_crve = FALSE))
  }
  
  m <- rma.mv(yi, vi, random = ~ 1 | study_id/effect_id, data = df, method = "REML")
  
  n_clusters <- dplyr::n_distinct(df$study_id)
  p_fixed    <- length(coef(m))  # usually 1
  
  if (n_clusters > p_fixed) {
    m_use <- robust_safe(m, cluster_vec = df$study_id)
    used_crve <- TRUE
  } else {
    m_use <- m                      # skip CRVE — too few clusters
    used_crve <- FALSE
  }
  
  est_log <- unname(coef(m_use)[1])
  ci_log  <- extract_ci(m_use)
  
  tibble(
    sheet          = sheet_name,
    k              = nrow(df),
    n_studies      = n_clusters,
    pooled_cv      = exp(est_log) * 100,
    lci            = exp(ci_log["lower"]) * 100,
    uci            = exp(ci_log["upper"]) * 100,
    tau2_between   = m$sigma2[1],
    sigma2_within  = m$sigma2[2],
    I2_studylevel_perc  = 100*m$sigma2[1]/(m$sigma2[1]+m$sigma2[2]+1),
    I2_withinstudy_perc = 100*m$sigma2[2]/(m$sigma2[1]+m$sigma2[2]+1),
    used_crve      = used_crve
  )
}

# run
# table 1
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 1.xlsx"   # e.g., "/path/to/Table2.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl <- purrr::map_dfr(sheets, analyze_sheet)
print(summary_tbl)
write.csv(summary_tbl, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_Table1.csv", row.names = FALSE)

# table 2
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 2.xlsx"   # e.g., "/path/to/Table2.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl <- purrr::map_dfr(sheets, analyze_sheet)
print(summary_tbl)
write.csv(summary_tbl, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_Table2.csv", row.names = FALSE)

