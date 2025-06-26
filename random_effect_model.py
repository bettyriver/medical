#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 22 11:44:58 2025

@author: ymai0110
"""
import numpy as np
from statsmodels.stats.meta_analysis import combine_effects
import pandas as pd
import scipy.stats as stats
from scipy.stats import norm

def random_effect_model(df,value_cname,se_cname,study_cname='Number'):
    '''
    

    Parameters
    ----------
    df : pd.dataframe
        Ddataframe that include value and it's SE.
    value_cname : str
        the name of the column of value (CV, ICC, Kappa, SRM).
    se_cname : str
        the name of the column of SE.
    study_cname : str, optional
        the name of the column of study. The default is 'Number'.

    Returns
    -------
    None.

    '''
    n_iter = 500
    pooled_estimates = []
    pooled_se = []

    for i in range(n_iter):
        # Randomly select one estimate per study
        sampled = df.groupby(study_cname).sample(n=1, replace=False)

        # obtain CV
        effects = sampled[value_cname].values
        se = sampled[se_cname].values 
        variances = se**2  # no se available, simply use the same variances

        # note: variance = SE^2, SE = SD/sqrt(n-1)

        
        # use DerSimonian-Laird to perform random effect combine
        
        meta_result = combine_effects(effect=effects, variance=variances, method_re="dl")
        summary_df = meta_result.summary_frame()
        if len(effects)<=3:
            pooled_cv = summary_df.loc['fixed effect', 'eff']
            pooled_sd = summary_df.loc['fixed effect', 'sd_eff']
        else:
            pooled_cv = summary_df.loc['random effect', 'eff']
            pooled_sd = summary_df.loc['random effect', 'sd_eff']
        
        pooled_estimates.append(pooled_cv)  # mean after combine
        pooled_se.append(pooled_sd/np.sqrt(len(effects)))  # SE after combine
        

    # calculate the pooled CV and SE after 500 iterations
    mean_pooled_cv = np.nanmean(pooled_estimates)
    mean_pooled_se = np.nanmean(pooled_se)

    # calculate 95% CI
    ci_low = mean_pooled_cv - 1.96 * mean_pooled_se
    ci_high = mean_pooled_cv + 1.96 * mean_pooled_se

    # output
    #print(f"Pooled CV = {mean_pooled_cv:.4f}")
    #print(f"95% CI = ({ci_low:.4f}, {ci_high:.4f})")
    return mean_pooled_cv, ci_low, ci_high

def REM_ICC(df, effect_name='ICC'):
    '''
    

    Parameters
    ----------
    df : pandas.dataframe
        dataframe for ICC
    effect_name : str, default is 'ICC'
        effect name, should be 'ICC'

    Returns
    -------
    None.

    '''
    
    df = get_value_and_se(df=df, effect_name=effect_name)
    
    
    icc_cname= effect_name + '_value_final'
    icc_se_cname= effect_name + '_SE_final'
    
    icc_value = df[icc_cname]
    icc_se = df[icc_se_cname]
    
    # Fisher's Z transform of ICC
    #z = 0.5 * np.log((1 + icc_value) / (1 - icc_value))
    z = fisher_z(icc_value)
    
    # Transform SE_ICC to SE_Z using delta method
    #z_se_from_icc = icc_se/(1 - icc_value**2)
    z_se_from_icc = se_fisher_z(r=icc_value, se_r=icc_se)
    
    # note, a lot of study doesn't provide rater number
    ## If SE_ICC not available, estimate SE_Z from sample size (Fisher approximation)
    #z_se_from_n = 1 / np.sqrt(n - 3)
    
    ## Use SE_Z_from_SE_ICC when available, otherwise use SE_Z_from_n
    #z_se = z_se_from_icc.combine_first(z_se_from_n)
    
    df['z'] = z
    df['z_se'] = z_se_from_icc
    
    z, z_ci_low, z_ci_high = random_effect_model(df=df,
                                value_cname='z',
                                se_cname='z_se',
                                study_cname='Number')
    
    
    # back-transform
    pooled_icc = inverse_fisher_z(z)
    pooled_ci_low = inverse_fisher_z(z_ci_low)
    pooled_ci_high = inverse_fisher_z(z_ci_high)
    
    # output
    print(f"Pooled ICC = {pooled_icc:.4f}")
    print(f"95% CI = ({pooled_ci_low:.4f}, {pooled_ci_high:.4f})")
    return pooled_icc, pooled_ci_low, pooled_ci_high

def REM_Kappa(df, effect_name='Kappa'):
    '''
    

    Parameters
    ----------
    df : pandas.dataframe
        dataframe for Kappa
    effect_name : str, default is 'Kappa'
        effect name, should be 'Kappa'

    Returns
    -------
    None.

    '''
    
    df = get_value_and_se(df=df, effect_name=effect_name)
    
    
    kappa_cname= effect_name + '_value_final'
    kappa_se_cname= effect_name + '_SE_final'
    
    kappa_value = df[kappa_cname]
    kappa_se = df[kappa_se_cname]
    
    # Fisher's Z transform of kappa
    z = fisher_z(kappa_value)
    
    # Transform SE_ICC to SE_Z using delta method
    z_se_from_kappa = se_fisher_z(r=kappa_value, se_r=kappa_se)
    
    
    df['z'] = z
    df['z_se'] = z_se_from_kappa
    
    z, z_ci_low, z_ci_high = random_effect_model(df=df,
                                value_cname='z',
                                se_cname='z_se',
                                study_cname='Number')
    
    
    # back-transform
    pooled_kappa = inverse_fisher_z(z)
    pooled_ci_low = inverse_fisher_z(z_ci_low)
    pooled_ci_high = inverse_fisher_z(z_ci_high)
    
    # output
    print(f"Pooled Kappa = {pooled_kappa:.4f}")
    print(f"95% CI = ({pooled_ci_low:.4f}, {pooled_ci_high:.4f})")
    return pooled_kappa,pooled_ci_low, pooled_ci_high

def REM_CV(df, effect_name='CV'):
    '''
    

    Parameters
    ----------
    df : pandas.dataframe
        dataframe for CV
    effect_name : str, default is 'CV'
        effect name, should be 'CV'

    Returns
    -------
    None.

    '''
    
    df = get_value_and_se(df=df, effect_name=effect_name)
    
    
    cv_cname= effect_name + '_value_final'
    cv_se_cname= effect_name + '_SE_final'
    
    cv_value = df[cv_cname]
    cv_se = df[cv_se_cname]
    
    # log transform of cv
    log_cv_value = log_cv(cv_value)
    
    # log transform of cv SE
    se_log_cv_value = se_log_cv(cv=cv_value, se_cv=cv_se)
    
    
    df['z'] = log_cv_value
    df['z_se'] = se_log_cv_value
    
    z, z_ci_low, z_ci_high = random_effect_model(df=df,
                                value_cname='z',
                                se_cname='z_se',
                                study_cname='Number')
    
    
    # back-transform
    pooled_cv = inverse_log_cv(z)
    pooled_ci_low = inverse_log_cv(z_ci_low)
    pooled_ci_high = inverse_log_cv(z_ci_high)
    
    # output
    print(f"Pooled CV = {pooled_cv:.4f}")
    print(f"95% CI = ({pooled_ci_low:.4f}, {pooled_ci_high:.4f})")
    return pooled_cv, pooled_ci_low, pooled_ci_high

def REM_SRM(df, effect_name='SRM'):
    '''
    

    Parameters
    ----------
    df : pandas.dataframe
        dataframe for SRM
    effect_name : str, default is 'SRM'
        effect name, should be 'SRM'

    Returns
    -------
    None.

    '''
    
    df = get_value_and_se(df=df, effect_name=effect_name)
    
    
    srm_cname= effect_name + '_value_final'
    srm_se_cname= effect_name + '_SE_final'
    
    srm_value = df[srm_cname]
    srm_se = df[srm_se_cname]
    
    
    
    pooled_srm, srm_ci_low, srm_ci_high = random_effect_model(df=df,
                                value_cname=srm_cname,
                                se_cname=srm_se_cname,
                                study_cname='Number')
    
    # output
    print(f"Pooled SRM = {pooled_srm:.4f}")
    print(f"95% CI = ({srm_ci_low:.4f}, {srm_ci_high:.4f})")
    return pooled_srm, srm_ci_low, srm_ci_high
    

def fisher_z(r):
    """Apply Fisher Z-transformation"""
    return 0.5 * np.log((1 + r) / (1 - r))

def se_fisher_z(r, se_r):
    """Apply Fisher Z-transformation for SE"""
    return se_r / (1 - r**2)

def log_cv(cv_percent):
    """Apply log transformation to CV%"""
    return np.log(cv_percent / 100)

def se_log_cv(cv, se_cv):
    """
    Calculate SE of log(CV) given CV and SE of CV.
    Both cv and se_cv are in the same units (e.g., percentages).
    """
    return se_cv / cv
    
    
def inverse_fisher_z(z):
    """Back-transform Fisher's Z to r (ICC or Kappa)"""
    return (np.exp(2*z) - 1) / (np.exp(2*z) + 1)

def inverse_log_cv(log_cv):
    """Back-transform log(CV/100) to CV in percentage"""
    return 100 * np.exp(log_cv)    


def get_value_and_se(df,effect_name):
    
    
    value_cn = f"{effect_name}_value"
    value_lowlim_cn = f"{effect_name}_Lowerlimits"
    value_uplim_cn = f"{effect_name}_upperlimits"
    value_sd_cn = f"{effect_name}_SD"
    value_se_cn = f"{effect_name}_SE" 
    value_lowCI_cn = f"{effect_name}_LCL" 
    value_upCI_cn = f"{effect_name}_UCL"
    value_samplesize_cn = f"{effect_name}_sample_size"
    rater_cn = "Number_of_measurements"
    
    
    
    
    cols_to_convert = df.columns.difference(['Number'])
    
    #df = df.apply(pd.to_numeric, errors='coerce')
    df[cols_to_convert] = df[cols_to_convert].apply(pd.to_numeric, errors='coerce')
    
    # get the value or calculate the value from low/up limit
    value = df[value_cn]
    value_lowlim = df[value_lowlim_cn]
    value_uplim = df[value_uplim_cn]
    mask = df[value_cn].isna()
    #value[mask] = (value_lowlim[mask] + value_uplim[mask]) / 2
    df.loc[mask, value_cn] = (df.loc[mask, value_lowlim_cn] + df.loc[mask, value_uplim_cn]) / 2
    
    if effect_name=='ICC' or effect_name=='Kappa': # kappa=1 will cause error in fisher z transform
        query = df[value_cn]==1
        df.loc[query, value_cn] = 0.999 # icc=1 will cause estimate se=0, try to avoid it
    
    df[value_cn+'_final'] = df[value_cn]
    
    df = estimate_se_from_effect(df=df, effect_name=effect_name)
    
    return df
    
    
    

def estimate_se_from_samplesize(data_type, value, sample_size,rater=None):
    '''
    

    Parameters
    ----------
    data_type : str
        DESCRIPTION.
    value : float or array
        value (ICC, SRM, CV, Kappa) of estimate.
    sample_size : float or array
        sample size.
    rater : float or array, optional
        the number of rater or how many repeat within one reader. 
        This is needed for ICC. no need for others.
        The default is None.

    Returns
    -------
    se : float or array
        standard error estimation for given value.

    '''
    if data_type=='ICC':
        # Bonett+2002, acccurate when sample_size >= 30
        se = np.sqrt(2*((1-value)**2)*((1+(rater-1)*value)**2)/(rater*(sample_size-1)*(rater-1)))
    
    if data_type == 'SRM':
        se = np.sqrt(1/sample_size + value**2/(2*sample_size))
    
    if data_type == 'CV':
        se = np.sqrt(1/(2*(sample_size - 1))) *100 # percentage
    
    if data_type == 'Kappa':
        se = (1- value**2)/np.sqrt(sample_size)
    
    return se


def estimate_se_from_effect(df, effect_name):
    """
    Estimate the SE for the given effect (e.g., 'cv', 'icc', 'kappa') using the following logic:
    1. Use reported SE
    2. Use CI if available
    3. Use SD and n
    4. Estimate SE using assumed SD = 1.0 and n
    """
    value_col = f"{effect_name}_value"
    se_col = f"{effect_name}_SE"
    ci_low_col = f"{effect_name}_LCL"
    ci_high_col = f"{effect_name}_UCL"
    sd_col = f"{effect_name}_SE"
    samplesize_col = f"{effect_name}_sample_size"

    def get_se(row):
        # Step 1: Use reported SE
        if pd.notnull(row.get(se_col)):
            return row[se_col]
        
        # Step 2: Estimate SE from CI
        if pd.notnull(row.get(ci_low_col)) and pd.notnull(row.get(ci_high_col)):
            return (row[ci_high_col] - row[ci_low_col]) / (2 * 1.96)
        
        # Step 3: Estimate SE from SD and sample size
        if pd.notnull(row.get(sd_col)) and pd.notnull(row.get(samplesize_col))and row[samplesize_col] > 0:
            return row[sd_col] / np.sqrt(row[samplesize_col])
        
        # Step 4: Fallback estimate using sample size
        if pd.notnull(row.get(samplesize_col)) and row[samplesize_col] > 0:
            if effect_name == 'ICC':
                se = estimate_se_from_samplesize(data_type='ICC', 
                                                 value=row[value_col], 
                                                 sample_size=row[samplesize_col],
                                                 rater=row['Number_of_measurements'])
            else:
                se = estimate_se_from_samplesize(data_type=effect_name, 
                                                 value=row[value_col], 
                                                 sample_size=row[samplesize_col])
            return se
        
        return np.nan

    # Create new SE column
    df[f"{effect_name}_SE_final"] = df.apply(get_se, axis=1)
    return df


def one_study_result(df,effect_name):
    '''
    if one table only have one study, then use this to get mean and CI

    Parameters
    ----------
    df : pd.dataframe
        dataframe for the data.
    effect_name : str
        effect name

    Returns
    -------
    None.

    '''
    
    
    
    df = get_value_and_se(df=df,effect_name=effect_name)
    
    value_cn = f"{effect_name}_value"
    value = df[value_cn+'_final']
    se_cn = f"{effect_name}_SE_final"
    se = df[se_cn]
    
    if len(value)<=3: # if number of estimates too small, get_mean_and_CI() will be wrong
        value = np.mean(value)
        se = np.mean(se)
        if effect_name == 'ICC' or effect_name =='Kappa':
            z = fisher_z(value)
            z_se = se_fisher_z(r=value, se_r=se)
            
            # 95% CI in Z space
            #z_crit = norm.ppf(0.975)
            z_lower = z - 1.96 * z_se
            z_upper = z + 1.96 * z_se
            
            ci_lower = inverse_fisher_z(z_lower)
            ci_upper = inverse_fisher_z(z_upper)
            mean = value
        elif effect_name == 'CV':
            z = log_cv(value)
            z_se = se_log_cv(cv=value, se_cv=se)
            
            z_lower = z - 1.96 * z_se
            z_upper = z + 1.96 * z_se
            
            ci_lower = inverse_log_cv(z_lower)
            ci_upper = inverse_log_cv(z_upper)
            mean = value
            
            
        elif effect_name == 'SRM':
            ci_lower = value - 1.96 * se
            ci_upper = value + 1.96 * se
            mean = value
        
        ci_lower = ci_lower
        ci_upper = ci_upper
        print(f"Mean {effect_name}: {mean:.4f}")
        print(f"95% CI: ({ci_lower:.4f}, {ci_upper:.4f})")
        
        return mean, ci_lower, ci_upper
    
    
    if effect_name == 'ICC' or effect_name == 'Kappa':
        
        # transform
        z = fisher_z(r=value)
        
        mean_z, ci_lower_z, ci_upper_z = get_mean_and_CI(z)
        
        # back-transform
        mean = inverse_fisher_z(mean_z)
        ci_lower = inverse_fisher_z(ci_lower_z)
        ci_upper = inverse_fisher_z(ci_upper_z)
        
        
    elif effect_name == 'CV':
        # transform
        cv_trans = log_cv(value)
        mean_z, ci_lower_z, ci_upper_z = get_mean_and_CI(cv_trans)
        
        # back-transform
        mean = inverse_log_cv(mean_z)
        ci_lower = inverse_log_cv(ci_lower_z)
        ci_upper = inverse_log_cv(ci_upper_z)
    elif effect_name == 'SRM':
        mean, ci_lower, ci_upper = get_mean_and_CI(value)
        
    print(f"Mean {effect_name}: {mean:.4f}")
    print(f"95% CI: ({ci_lower:.4f}, {ci_upper:.4f})")
    return mean, ci_lower, ci_upper
        
def get_mean_and_CI(data):
    # Calculate mean
    mean = np.mean(data)
    
    # 95% confidence interval
    confidence = 0.95
    n = len(data)
    se = stats.sem(data)  # Standard error of the mean
    h = se * stats.t.ppf((1 + confidence) / 2., n-1)  # Margin of error using t-distribution
    
    ci_lower = mean - h
    ci_upper = mean + h
    
    return mean, ci_lower, ci_upper