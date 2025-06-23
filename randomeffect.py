#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 17 22:33:54 2025

@author: Yifan Mai
"""

import numpy as np
import pandas as pd

parent_path = '/Users/ymai0110/Documents/medical_data/folder_June17/'
df = pd.read_excel(parent_path+'table2.xlsx', sheet_name='1cartilage')

cv = pd.to_numeric(df['CV_value_or_Lowerlimits'], errors='coerce')

n_iter = 500
mean_estimates = []

# Repeat 500 times
for _ in range(n_iter):
    # Randomly select one estimate per study
    sampled = df.groupby('Number').sample(n=1, replace=False)

    # Compute the average of selected estimates
    mean_val = sampled['cv_value_num'].mean()
    mean_estimates.append(mean_val)

# Convert to NumPy array for easier math
mean_estimates = np.array(mean_estimates)

# Final results
final_mean = mean_estimates.mean()
final_se = mean_estimates.std(ddof=1)  # standard error of the 500 means
z = 1.96
ci_lower = final_mean - z * final_se
ci_upper = final_mean + z * final_se

print(f"Final average of means: {final_mean:.3f}")
print(f"Standard error: {final_se:.3f}")
print(f"95% CI: ({ci_lower:.3f}, {ci_upper:.3f})")