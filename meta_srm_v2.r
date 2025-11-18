# packages (as before)
install.packages(c("readxl", "dplyr", "purrr", "stringr", "tibble", "metafor", "writexl"))
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
  est <- unname(coef(fit)[1]); se <- sqrt(vcov(fit)[1,1])
  c(lower = est - qnorm(0.975)*se, upper = est + qnorm(0.975)*se)
}

analyze_sheet_srm <- function(sheet_name, path = excel_path) {
  raw <- readxl::read_excel(path, sheet = sheet_name)
  
  req <- c("Number","SRM_value_final","SRM_SE_final","effect_id")
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0) stop(sprintf("Sheet '%s' missing: %s", sheet_name, paste(miss, collapse=", ")))
  
  df <- raw %>%
    transmute(
      study_id  = as.character(Number),
      effect_id = as.character(effect_id),
      SRM       = as.numeric(SRM_value_final),  # standardized response mean
      SE        = as.numeric(SRM_SE_final)
    ) %>%
    mutate(
      yi = SRM,
      vi = SE^2
    ) %>%
    filter(is.finite(yi), is.finite(vi))
  
  if (nrow(df) < 2) {
    return(tibble(sheet = sheet_name, k = nrow(df),
                  n_studies = dplyr::n_distinct(df$study_id),
                  pooled_srm = NA_real_, lci = NA_real_, uci = NA_real_,
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
  
  est <- unname(coef(m_use)[1])
  ci  <- extract_ci(m_use)
  V_sampling   <- mean(df$vi)
  tibble(
    sheet          = sheet_name,
    k              = nrow(df),
    n_studies      = n_clusters,
    pooled_srm     = est,
    lci            = ci["lower"],
    uci            = ci["upper"],
    tau2_between   = m$sigma2[1],
    sigma2_within  = m$sigma2[2],
    I2_studylevel_perc  = 100*m$sigma2[1]/(m$sigma2[1]+m$sigma2[2]+V_sampling),
    I2_withinstudy_perc = 100*m$sigma2[2]/(m$sigma2[1]+m$sigma2[2]+V_sampling),
    used_crve      = used_crve
  )
}

# run
# table 7 - n 
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 7quantative negative.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table7quantative negative.csv", row.names = FALSE)

# table 7 - p
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 7quantative positive.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table7quantative positive.csv", row.names = FALSE)

# table 8 - n 
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 8Semi-quantative negative.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table8Semi-quantative negative.csv", row.names = FALSE)


# table 8 - p
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 8Semi-quantative positive.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table8Semi-quantative positive.csv", row.names = FALSE)

# table 9
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 9.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table9.csv", row.names = FALSE)

# table 10
excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 10.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_Table10.csv", row.names = FALSE)

# FNIH

excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/FNIH SRM.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_FNIHSRM.csv", row.names = FALSE)

# treatment

excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/treatment SRM.xlsx"
sheets <- readxl::excel_sheets(excel_path)
sheets
summary_tbl_srm <- purrr::map_dfr(sheets, analyze_sheet_srm)
print(summary_tbl_srm)
write.csv(summary_tbl_srm, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_v2/result_treatmentSRM.csv", row.names = FALSE)

