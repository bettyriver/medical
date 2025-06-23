#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 22:44:22 2025

@author: Yifan Mai
"""

import numpy as np
import pandas as pd
from statsmodels.stats.meta_analysis import combine_effects

parent_path = '/Users/ymai0110/Documents/medical_data/folder_June18/'
df = pd.read_excel(parent_path+'table2.xlsx', sheet_name='1cartilage')
cv = pd.to_numeric(df['CV_value'], errors='coerce')
cv_low = pd.to_numeric(df['CV_Lowerlimits'], errors='coerce')
cv_high = pd.to_numeric(df['CV_upperlimits'], errors='coerce')
mask = np.isnan(cv)
cv[mask] = (cv_low[mask] + cv_high[mask]) / 2
df['cv_value_num'] = cv

n_iter = 500
pooled_estimates = []
pooled_se = []

for i in range(n_iter):
    # Randomly select one estimate per study
    sampled = df.groupby("Number").sample(n=1, replace=False)

    # obtain CV
    effects = sampled['cv_value_num'].values
    variances = np.ones_like(effects)  # no se available, simply use the same variances

    # note: variance = SE^2, SE = SD/sqrt(n-1)

    
    # use DerSimonian-Laird to perform random effect combine
    
    meta_result = combine_effects(effect=effects, variance=variances, method_re="dl")
    summary_df = meta_result.summary_frame()
    pooled_cv = summary_df.loc['random effect', 'eff']
    pooled_sd = summary_df.loc['random effect', 'sd_eff']
    
    pooled_estimates.append(pooled_cv)  # mean after combine
    pooled_se.append(pooled_sd/np.sqrt(len(effects)-1))  # SE after combine
    

# calculate the pooled CV and SE after 500 iterations
mean_pooled_cv = np.mean(pooled_estimates)
mean_pooled_se = np.mean(pooled_se)

# calculate 95% CI
ci_low = mean_pooled_cv - 1.96 * mean_pooled_se
ci_high = mean_pooled_cv + 1.96 * mean_pooled_se

# output
print(f"Pooled CV = {mean_pooled_cv:.4f}")
print(f"95% CI = ({ci_low:.4f}, {ci_high:.4f})")