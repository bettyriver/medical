library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(tibble)
library(metafor)
library(writexl)

analyze_sheet_srm_meta_reg <- function(sheet_name, path = excel_path, time_unit = c("years","months")) {
  time_unit <- match.arg(time_unit)
  raw <- readxl::read_excel(path, sheet = sheet_name)
  
  req <- c("Number","SRM_value_final","SRM_SE_final","effect_id","study_duration_months")
  miss <- setdiff(req, names(raw))
  if (length(miss) > 0) stop(sprintf("Sheet '%s' missing: %s", sheet_name, paste(miss, collapse = ", ")))
  
  df <- raw %>%
    transmute(
      study_id  = as.character(Number),
      effect_id = as.character(effect_id),
      SRM       = as.numeric(SRM_value_final),
      SE        = as.numeric(SRM_SE_final),
      follow_up = as.numeric(study_duration_months)
    ) %>%
    mutate(
      yi = SRM,
      vi = SE^2
    ) %>%
    filter(is.finite(yi), is.finite(vi), is.finite(follow_up))
  
  if (nrow(df) < 2) {
    return(tibble(sheet = sheet_name, k = nrow(df),
                  n_studies = dplyr::n_distinct(df$study_id),
                  pooled_srm = NA_real_, lci = NA_real_, uci = NA_real_,
                  tau2_between = NA_real_, sigma2_within = NA_real_,
                  I2_studylevel_perc = NA_real_, I2_withinstudy_perc = NA_real_,
                  used_crve = FALSE,
                  time_coef = NA_real_, time_lci = NA_real_, time_uci = NA_real_,
                  adj_srm_1yr = NA_real_, adj_lci_1yr = NA_real_, adj_uci_1yr = NA_real_))
  }
  
  # define time moderator in years; center at 1 year so intercept = SRM at 1 year
  time_years <- if (time_unit == "months") df$follow_up / 12 else df$follow_up
  df$time_c <- time_years - 1
  
  m <- rma.mv(yi, vi,
              mods = ~ time_c,
              random = ~ 1 | study_id/effect_id,
              data = df, method = "REML")
  
  # robust small-sample cluster-robust inference by study
  n_clusters <- dplyr::n_distinct(df$study_id)
  p_fixed <- length(coef(m))
  if (n_clusters > p_fixed) {
    m_use <- robust_safe(m, cluster_vec = df$study_id)
    used_crve <- TRUE
  } else {
    m_use <- m
    used_crve <- FALSE
  }
  
  # extract overall intercept (adjusted SRM at 1 year) and slope (change per +1 year)
  b <- coef(m_use)
  vc <- vcov(m_use)
  se_b0 <- sqrt(vc[1,1]); se_b1 <- sqrt(vc[2,2])
  
  # generic 95% CI for coefficients (falls back if confint not available)
  ci_int <- extract_ci(m_use)  # returns CI for first coefficient (intercept at 1 year)
  # get CI for slope explicitly
  z <- qnorm(0.975)
  time_lci <- b[2] - z * se_b1
  time_uci <- b[2] + z * se_b1
  
  # predicted (adjusted) SRM at exactly 1 year and its CI via predict()
  newdat <- data.frame(time_c = 0)
  pr <- predict(m, newmods = model.matrix(~ time_c, newdat), transf = NULL)
  # If robust CI desired for predictions, approximate via delta method using vcov:
  adj_srm_1yr <- as.numeric(b[1])
  adj_lci_1yr <- as.numeric(ci_int["lower"])
  adj_uci_1yr <- as.numeric(ci_int["upper"])
  
  tibble(
    sheet = sheet_name,
    k = nrow(df),
    n_studies = n_clusters,
    # overall random-effects decomposition from the fitted (non-robust) model
    tau2_between = m$sigma2[1],
    sigma2_within = m$sigma2[2],
    I2_studylevel_perc  = 100 * m$sigma2[1] / (m$sigma2[1] + m$sigma2[2] + 1),
    I2_withinstudy_perc = 100 * m$sigma2[2] / (m$sigma2[1] + m$sigma2[2] + 1),
    used_crve = used_crve,
    # moderator results
    time_coef = unname(b[2]),
    time_lci  = unname(time_lci),
    time_uci  = unname(time_uci),
    # adjusted SRM at 1 year (intercept due to centering)
    adj_srm_1yr = adj_srm_1yr,
    adj_lci_1yr = adj_lci_1yr,
    adj_uci_1yr = adj_uci_1yr
  )
}

##
rm(df)              # removes the function/object named df
# or to be safe, restart R and run the script top-to-bottom

dat <- raw %>%
  transmute(
    study_id  = as.character(Number),
    effect_id = as.character(effect_id),
    yi        = as.numeric(SRM_value_final),
    vi        = as.numeric(SRM_SE_final)^2,
    time_c    = (as.numeric(follow_up)/12) - 1
  ) %>%
  filter(is.finite(yi), is.finite(vi), is.finite(time_c))

class(dat)   # expect "data.frame" or "tbl_df"
str(dat[, c("yi","vi","time_c","study_id","effect_id")])
m <- rma.mv(yi, vi, mods = ~ time_c,
            random = ~ 1 | study_id/effect_id,
            data = dat, method = "REML")



#


m_grp <- rma.mv(yi, vi,
                mods = ~ time_c ,             # additive; add interactions if needed
                random = ~ 1 | study_id/effect_id,
                data = df, method = "REML")
m_grp_r <- robust_safe(m_grp, cluster_vec = df$study_id)

# Estimated adjusted SRM at 1 year for each tissue:
mm <- model.matrix(~ time_c , data.frame(time_c = 0))
pred <- predict(m_grp, newmods = mm)
# pred$pred gives adjusted SRM at 1 year for each level; use robust SEs for tests if needed

excel_path <- "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/Table 10.xlsx"
sheets <- readxl::excel_sheets(excel_path)
res <- purrr::map_dfr(sheets, ~ analyze_sheet_srm_meta_reg(.x, path = excel_path, time_unit = "months"))
write.csv(res, "/Users/ymai0110/Documents/medical_data/R/Rtable_Nov02/result_Table10_metaregression.csv", row.names = FALSE)
