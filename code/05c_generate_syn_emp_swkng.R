## Synthesizing employment with OPM and smooth quantile regression loss
# This code optimize the private loss function in OPM subject to linear 
# constraints using constrOptim. The linear constraints (upper bound on 
# employment based on the public data, lower bound on employment  = 0) 
# are set to improve the utility of synthetic data.

# The main difference between this version and OPM_constraints is that it
# does not use linear constraints to ensure the output satisfies >= 0
# requirement but instead use post-processing, hence the name nc (no constraints).

# Initial input: ./output/[sic3]/[sic3]_syn_data_year_mu.csv
# Final output: ./output/[sic3]/[sic3]_syndata_swkng_[kernel].csv
# Intermediate input and output will have similar syntax, but with an extra
# year parameter.

# Last updated: 4/27/24

data_jitter <- function(data){
    data = as.data.frame(data)
    for(i in 1:ncol(data)){
        names=names(data)
        tmp=data[,names[i]]
        tmp=tmp + runif(length(tmp),-1,1)
        data[,names[i]]=tmp
    }
    return(data)
}

de_id = function(x) {
    sep_loc = gregexpr('_', x)[[1]]
    lbdnum_id = as.numeric(substr(x, 1, sep_loc[1] -1))
    fye = as.numeric(substr(x, sep_loc[1] + 1, sep_loc[2] -1))
    lye = as.numeric(substr(x, sep_loc[2] + 1, nchar(x)))
    return(c(lbdnum_id, fye, lye))
}

generate_syn_emp = function(raw_data, syn_data, syn_var, taus, main_taus, h, epsilon, Cx, kernel, mod) {
    curr_n = nrow(syn_data)
    
    raw_data = data_jitter(raw_data)
    scale = ifelse(mod == 'continuer' & syn_var != 'emp_1977', 1e-2, 1e-2)
    year = as.numeric(substr(syn_var, 5, 8))
    if (mod == 'birth' | syn_var == 'emp_1977') {
      check_data = rbind(c(0, 0), c(0, 1), c(2000-year, 0), c(2000-year, 1))
    } else {
      check_data = rbind(c(0, 0), c(0, Cx), 
                         c(2000-year, 0), c(2000-year, Cx))
    }
    tmp = sandwichKNG(data = raw_data, total_eps = epsilon, tau = taus, main_tau = main_taus, scale = scale, #1e-1,
                      Cx = Cx, nbatch = 100000, method = 'varying', h = h, lb = 0, 
                      ub = Cx, check_data = check_data, lower_accept = 0.1, upper_accept = 0.7,
                      update_after = 25, adjust_scale_by = 2)[[1]]
    if (is.null(tmp)) {
        return(NULL)
    }
    
    # generate synthetic data
    samp = sample(1:length(taus), curr_n, replace = TRUE, prob = rep(0.1, length(taus)))
    y_hat = as.matrix(cbind(1, syn_data)) %*% tmp
    coord = cbind(c(1:curr_n), samp)
    ans = round(y_hat[coord]) #y_hat[coord] #
    ans[ans < 0] = 0
    ans = as.data.frame(ans)
    syn_data_id = rownames(syn_data)
    ans = cbind(ans, t(sapply(syn_data_id, de_id)))
    colnames(ans) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
    
    return(ans)
    
}

syn_data_all_years = function(seed) {
    start_year = ifelse(seed == 1, 1978, 1977)
    for (curr_year in start_year:2000) { # rmb to change back to 2021
        # set variable to synthesize and predictors
        syn_var = paste0('emp_', curr_year)
        print(curr_year)
        print(Sys.time())
        
        mu_var = paste0('mu_', mu_ref$mu_year[mu_ref$year == curr_year])
        prev = paste0('emp_', curr_year-1)
        if (curr_year == 1977) {
            pred_var_cont = c('year_to_death', mu_var)
        } else {
            pred_var_cont = c('year_to_death', mu_var, prev) #'age'
        }    
        pred_var_birth = c('year_to_death', mu_var)
        
        # read in previously generated synthetic data
        # if year = 1977, use the general synthetic year and multiunit file (used across different DP mechs)
        # else use the specific DP mech synthetic dataset
        if (curr_year == 1977) { # remember to change to _seed, seed
            syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_seed", seed, ".csv"))
        } else {
            syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_swkng_", kernel, "_", curr_year-1, "_seed", seed, ".csv")) # change method
        }
        # syn_data = left_join(syn_data, st_ref, by = 'st')
        # syn_data$region[is.na(syn_data$region)] = 99
        
        ## split raw and synthetic data into birth and continuer establishments for modeling
        
        # Within continuer establishments, there are two cases
        # 1. Have valid previous year employment (lag1)
        # => We model current year employment based on previous year employment
        # 2. Do not have valid previous year employment (lag1)
        # => For now treat it as birth model (i.e., model without predictors)
        # => Later, if have time, merge with data from lag2-employment.
        # If use lag2-employment if it exists, else treat like birth model.
        # Note: This solution can only starts from 1979 (first year with lag2 info)
        
        raw_data_cont = raw_data[raw_data$firstyear < curr_year & 
                                     raw_data$lastyear >= curr_year &
                                     !is.na(raw_data[, syn_var]), ]
        raw_data_birth = raw_data[raw_data$firstyear == curr_year & 
                                      raw_data$lastyear >= curr_year &
                                      !is.na(raw_data[, syn_var]), ]
        
        # split synthetic data into birth, continuer, and NA for modeling
        syn_data_cont = syn_data[syn_data$firstyear < curr_year & 
                                     syn_data$lastyear >= curr_year, ]
        syn_data_birth = syn_data[syn_data$firstyear == curr_year & 
                                      syn_data$lastyear >= curr_year, ]
        syn_data_na = syn_data[syn_data$firstyear > curr_year | 
                                   syn_data$lastyear < curr_year, ]
        
        if(curr_year > 1977) {
            if (curr_year == 1978 & sum(is.na(raw_data_cont[, prev])) > 0) {
                raw_data_cont_tmp = raw_data_cont[is.na(raw_data_cont[, prev]), ]
                raw_data_birth = rbind(raw_data_birth, raw_data_cont_tmp)
                raw_data_cont = raw_data_cont[!rownames(raw_data_cont) %in% rownames(raw_data_cont_tmp), ]
                
                # syn_data_birth = rbind(syn_data_birth, syn_data_cont[rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ])
                # syn_data_cont = syn_data_cont[!rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ]
                
            } else if (curr_year > 1978 & sum(is.na(raw_data_cont[, prev])) > 0) {
                prev_prev = paste0('emp_', curr_year-2)
                raw_data_cont_tmp = raw_data_cont[is.na(raw_data_cont[, prev]), ]
                raw_data_cont = raw_data_cont[!rownames(raw_data_cont) %in% rownames(raw_data_cont_tmp), ]
                # syn_data_cont_tmp = syn_data_cont[rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ]
                # syn_data_cont = syn_data_cont[!rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ]
                
                tmp = raw_data_cont_tmp[is.na(raw_data_cont_tmp[, prev_prev]), ]
                raw_data_birth = rbind(raw_data_birth, tmp)
                # syn_data_birth = rbind(syn_data_birth, syn_data_cont_tmp[rownames(syn_data_cont_tmp) %in% rownames(tmp), ])
                
                raw_data_cont_tmp = raw_data_cont_tmp[!rownames(raw_data_cont_tmp) %in% rownames(tmp), ]
                raw_data_cont_tmp[, prev] = raw_data_cont_tmp[, prev_prev]
                raw_data_cont = rbind(raw_data_cont, raw_data_cont_tmp)
                
                # syn_data_cont_tmp = syn_data_cont_tmp[!rownames(syn_data_cont_tmp) %in% rownames(tmp), ]
                # syn_data_cont_tmp[, prev] = syn_data_cont_tmp[, prev_prev]
                # syn_data_cont = rbind(syn_data_cont, syn_data_cont_tmp)        
                
            }
            
        }
        
        
        if(nrow(syn_data_na) > 0) {
            ## assign NA values for establishments that do not exist in the curr_year
            ans_na = cbind(NA, syn_data_na$lbdnum, syn_data_na$firstyear, syn_data_na$lastyear)
            ans_na = as.data.frame(ans_na)
            colnames(ans_na) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
            
        } else {
            ans_na = NULL
        }
        
        
        if (nrow(syn_data_birth) > 0) {
            print(paste('Year', curr_year, '- Birth model'))
            ## model and generate synthetic employment for establishments that were born 
            ## on the curr_year (birthers) by regions
            # emp_birthyear ~ most_current_mu + year_to_death
            
            # save as rownames to retrieve it easily later without impacting matrix multiplication
            # lbdnum will be used later on to merge to the main synthetic dataset
            row.names(syn_data_birth) = paste(syn_data_birth$lbdnum, syn_data_birth$firstyear, syn_data_birth$lastyear, sep = '_')
            
            # create year_to_death variable from lastyear and curr_year for birth model
            syn_data_birth$year_to_death = syn_data_birth$lastyear - curr_year
            raw_data_birth$year_to_death = raw_data_birth$lastyear - curr_year
            
            raw_data_birth = raw_data_birth[, c(pred_var_birth, syn_var)]
            syn_data_birth = syn_data_birth[, pred_var_birth]
            
            if (nrow (raw_data_birth) > 0 & nrow(syn_data_birth) > 0) {
                ans_birth = NULL
                while(is.null(ans_birth)) {
                    tic()
                    ans_birth = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_birth, 
                                                                       main_taus = c(0.1, 0.5, 0.9, 0.99),
                                                                       syn_data = syn_data_birth, syn_var = syn_var, 
                                                                       taus = taus, h = h, epsilon = epsilon, Cx = 25, 
                                                                       kernel = kernel, mod = 'birth')}, timeout= 5400), 
                                         error = function(cond) NULL)
                    toc()
                    
                }
                
            } else if(nrow(raw_data_birth) == 0 & nrow(syn_data_birth) > 0) {
                ## assign NA values for establishments that do not exist in the curr_year
                syn_data_id = rownames(syn_data_birth)
                ans_birth = cbind(0, t(sapply(syn_data_id, de_id)))
                colnames(ans_birth) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
            }
            
        } else {
            print(paste('Year', curr_year, '- Skip birth model'))
            ans_birth = NULL
        }
        
        
        
        if(nrow(syn_data_cont) > 0) {
            print(paste('Year', curr_year, '- Continuer model'))
            ## model and generate synthetic employment for establishments that were born 
            ## before the curr_year (continuers) by regions
            # emp_contyear ~ most_current_mu + age + year_to_death + prev
            
            # save as rownames to retrieve it easily later without impacting matrix multiplication
            # lbdnum will be used later on to merge to the main synthetic dataset
            row.names(syn_data_cont) = paste(syn_data_cont$lbdnum, syn_data_cont$firstyear, syn_data_cont$lastyear, sep = "_")
            # create year_to_death variable from lastyear and curr_year for birth model
            # syn_data_cont$age = curr_year - syn_data_cont$firstyear
            # raw_data_cont$age = curr_year - raw_data_cont$firstyear
            syn_data_cont$year_to_death = syn_data_cont$lastyear - curr_year
            raw_data_cont$year_to_death = raw_data_cont$lastyear - curr_year
            
            Cx_tmp = 2500
            raw_data_cont = raw_data_cont[, c(pred_var_cont, syn_var)]
            syn_data_cont = syn_data_cont[, pred_var_cont]
            
            if(syn_var == 'emp_1977') {
                ans_cont = NULL
                while(is.null(ans_cont)) {
                    tic()
                    ans_cont = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_cont, 
                                                                      syn_data = syn_data_cont, syn_var = syn_var, 
                                                                      main_taus = c(0.1, 0.5, 0.9, 0.99),
                                                                      taus = taus, h = h, epsilon = epsilon, Cx = max(Cx_tmp, 46), 
                                                                      kernel = kernel, mod = 'continuer')}, timeout= 3600), 
                                        error = function(cond) NULL)
                    toc()
                    
                }
                
            } else {
                raw_data_cont_tmp0 = raw_data_cont[raw_data_cont[,2] == 0, ]
                raw_data_cont_tmp1 = raw_data_cont[raw_data_cont[,2] == 1, ]
                syn_data_cont_tmp0 = syn_data_cont[syn_data_cont[,2] == 0, ]
                syn_data_cont_tmp1 = syn_data_cont[syn_data_cont[,2] == 1, ]
                
                raw_data_cont_tmp0 = raw_data_cont_tmp0[, -2]
                raw_data_cont_tmp1 = raw_data_cont_tmp1[, -2]
                syn_data_cont_tmp0 = syn_data_cont_tmp0[, -2]
                syn_data_cont_tmp1 = syn_data_cont_tmp1[, -2]
                
                
                
                
                if (nrow (raw_data_cont_tmp0) > 0 & nrow(syn_data_cont_tmp0) > 0) {
                    print(paste('Year', curr_year, '- Continuer model - Mu = 0'))
                    
                    ans_cont_tmp0 = NULL
                    while(is.null(ans_cont_tmp0)) {
                        tic()
                        ans_cont_tmp0 = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_cont_tmp0, 
                                                                               syn_data = syn_data_cont_tmp0, syn_var = syn_var, 
                                                                               main_taus = c(0.1, 0.5, 0.9, 0.99),
                                                                               taus = taus, h = h, epsilon = epsilon, Cx = max(Cx_tmp, 46), 
                                                                               kernel = kernel, mod = 'continuer')}, timeout= 3600), 
                                                 error = function(cond) NULL)
                        toc()
                        
                    }
                    
                } else if (nrow (raw_data_cont_tmp0) == 0 & nrow(syn_data_cont_tmp0) > 0) {
                    syn_data_id = rownames(syn_data_cont_tmp0)
                    ans_cont_tmp0 = cbind(0, t(sapply(syn_data_id, de_id)))
                    colnames(ans_cont_tmp0) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
                } else {
                    ans_cont_tmp0 = NULL
                }
                
                if (nrow (raw_data_cont_tmp1) > 0 & nrow(syn_data_cont_tmp1) > 0) {
                    print(paste('Year', curr_year, '- Continuer model - Mu = 1'))
                    ans_cont_tmp1 = NULL
                    while(is.null(ans_cont_tmp1)) {
                        tic()
                        ans_cont_tmp1 = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_cont_tmp1, 
                                                                               syn_data = syn_data_cont_tmp1, syn_var = syn_var, 
                                                                               main_taus = c(0.1, 0.5, 0.9, 0.99),
                                                                               taus = taus, h = h, epsilon = epsilon, Cx = max(Cx_tmp, 46), 
                                                                               kernel = kernel, mod = 'continuer')}, timeout= 3600), 
                                                 error = function(cond) NULL)
                        toc()
                        
                    }
                    
                } else if (nrow (raw_data_cont_tmp1) == 0 & nrow(syn_data_cont_tmp1) > 0) {
                    syn_data_id = rownames(syn_data_cont_tmp1)
                    ans_cont_tmp1 = cbind(0, t(sapply(syn_data_id, de_id)))
                    colnames(ans_cont_tmp1) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
                } else {
                    ans_cont_tmp1 = NULL
                }
                
                ans_cont = rbind(ans_cont_tmp0, ans_cont_tmp1)
                
            }
            
            
        } else {
            print(paste('Year', curr_year, '- Skip continuer model'))
            ans_cont = NULL
        }
        
        
        ans_all = rbind(ans_birth, ans_cont, ans_na)
        ans_all$lbdnum = as.numeric(ans_all$lbdnum)
        syn_data = left_join(syn_data, ans_all, by = c('lbdnum', 'firstyear', 'lastyear'))
        
        if (curr_year == 2000) {
            output_file = paste0("./output/", sic3, "/", sic3, "_syndata_swkng_", kernel, "_seed", seed, ".csv")
        } else {
            output_file = paste0("./output/", sic3, "/", sic3, "_syndata_swkng_", kernel, "_", curr_year, "_seed", seed, ".csv")
        }
        
        write.csv(syn_data, output_file, row.names = FALSE)
        
        old_output_file = paste0("./output/", sic3, "/", sic3, "_syndata_swkng_", kernel, "_", curr_year-1, "_seed", seed, ".csv")
        if (file.exists(old_output_file)) {
            file.remove(old_output_file)
        }
        
    }
    
}


library(dplyr)
library(caret)
library(quantreg)

sic3 = 239 #as.numeric(commandArgs(trailingOnly = TRUE)) #2371
sic2 = floor(sic3/10)

data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
data = data[data$sic3 == sic3, ]
data$lastyear[data$lastyear > 2000] = 2000
data$firstyear[data$firstyear < 1977] = 1977

n = nrow(data)
rownames(data) = c(1:n)

# st_ref = read.csv('st_ref.csv')
# st_ref = st_ref[, c('st', 'region')]
# maxemp_ref = read.csv('maxemp_ref_regionnaics2.csv')
# maxemp_ref = maxemp_ref[maxemp_ref$naics2 == floor(sic3/100), ]
mu_ref = data.frame(year = c(1977:2000), mu_year = c(1977, rep(seq(1982, 1997, 5), each = 5), rep(2002, 3)))

# rearrange columns
data = data[, c('lbdnum', 'firstyear', 'lastyear', paste('mu5', seq(1977, 2002, 5), sep = '_'), 
                paste('emp', c(1977:2000), sep = '_'))]
colnames(data) = c('lbdnum', 'firstyear', 'lastyear', paste('mu', seq(1977, 2002, 5), sep = '_'), 
                   paste('emp', c(1977:2000), sep = '_'))
taus = c(seq(0.1, 0.9, 0.1), 0.99)
h = 1
#h = 1.06*min(sse, iqr/1.34)*n^(-1/5)
kernel = 'logistic'

raw_data = data
# raw_data = left_join(raw_data, st_ref, by = 'st')
# raw_data$region[is.na(raw_data$region)] = 99

epsilon = 0.25

# 
# library(doParallel)
# library(tictoc)
# num_cores= 5 #detectCores()-1 #use all available core but 1
# print(num_cores)
# workers=makeCluster(num_cores, output_file = paste0('05b_generate_syn_emp_swkng_', sic3, '.Rout'))
# registerDoParallel(workers)
# 
# tic()
# Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
seed = as.integer(commandArgs(trailingOnly = TRUE)) #2371 2379 #
# load libraries
library(dplyr)
library(caret)
library(quantreg)
library(tictoc)
library(R.utils)
source('./code/source_functions_kng1.R')


set.seed(seed)
syn_data_all_years(seed)
# }
# toc()
# stopCluster(workers)

