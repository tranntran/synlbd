kng_qr_logistic = function(theta, x, y, tau, h, epsilon = 1, Cx = 1) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        exp_val = exp(-(x[i,] %*% theta - y[i])/h)
        tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = max(1-tau, tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M") - 1/2*norm(as.matrix(theta), "F")^2
    return(ans)   
}

kng_qr_logistic_constr = function(theta, x, y, tau, h, epsilon = 1, Cx = 1, check_theta, check_data) {
    n = nrow(x)
    tmp = 0
    for (i in 1:n) {
        exp_val = exp(-(x[i,] %*% theta - y[i])/h)
        tmp = tmp + (1/(1 + exp_val) - tau) %*% x[i, ]
    }
    delta = max(1-tau, tau)*Cx
    ans = -epsilon / (2*delta) * norm(tmp, "M") - 1/2*norm(as.matrix(theta), "F")^2
    
    check = all(check_data %*% as.matrix(theta) <= check_data %*% check_beta) & all(check_data %*% as.matrix(theta) >= 0)
    if (!check) {
        return(-Inf)
    } else {
        return(ans)   
}
}


stepwiseKNG = function(data, total_eps, median_eps = NULL, tau, scale = 1e-4, Cx,
                       nbatch = 10000, method = c("fixed", "varying"), h = 1,
                       ub = Inf, lb = -Inf, check_data = NULL, lower_accept = 0,
                       upper_accept = 1, update_after = 10, adjust_scale_by = 2,
                       formula = NULL){
    method = match.arg(method)
    if(is.null(median_eps)){
        median_eps = ifelse(method == "fixed", 0.3, 0.4)
    }
    ep = total_eps*(1-median_eps)/(length(tau)-1)
    i = ncol(data)
    Y = data[,i]
    X = as.matrix(cbind(rep(1, nrow(data)), data))
    X = as.matrix(X[, -ncol(X)])
    X[X > Cx] = Cx
    m = ncol(X) - 1
    
    if (is.null(formula)) {
        vars = colnames(data)
        formula = paste(vars[i], " ~ .")
    }
    
    print('0.5')
    nonpriv = quantreg::rq(formula, data = as.data.frame(data), tau = 0.5)
    start_points = matrix(rep(coef(nonpriv), 3), nrow= 3)
    tmp = MCMC(start_points, kng_qr_logistic,nsteps = nsteps, nchains = 3, x = X,  Cx = Cx,
               y = Y, tau = 0.4, epsilon = 0.1, h = 1, conv_checker = convergence_gelman(),
               kernel = kernel_ram())
    ans = t(tail(as.matrix(tmp[[1]]), 1))
    summary(check_data%*%ans)
    
    median_beta_kng = ans
    
    tau_lower = sort(tau[tau < 0.5], decreasing = TRUE)
    tau_upper = sort(tau[tau > 0.5])
    
    if (length(tau_lower) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_lower)){
            new_tau = tau_lower[i]
            print(new_tau)
            start_points = rbind(tail(as.matrix(tmp[[1]]), 1), tail(as.matrix(tmp[[2]]), 1), tail(as.matrix(tmp[[3]]), 1))
            tmp = MCMC(start_points, kng_qr_logistic,nsteps = nsteps, nchains = 3, x = X,  Cx = Cx,
                       y = Y, tau = new_tau, epsilon = ep, h = 1, conv_checker = convergence_gelman(),
                       kernel = kernel_ram())
            
            out = constrKNG(init = list(check_beta-1e-2, rep(0, m+1), check_beta/2), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 3,
                            check_data = check_data, ub = ub, lb = lb, method = method,
                            type = "lower", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            ans = cbind(t(proposed_beta), ans)
            scale_output = cbind(tail(out[[3]][[1]], 1), scale_output)
            accept_rate = cbind(out[[2]][[1]], accept_rate)
            
        }
    }
    
    if (length(tau_upper) > 0){
        check_beta = median_beta_kng
        for (i in 1:length(tau_upper)){
            new_tau = tau_upper[i]
            print(new_tau)
            out = constrKNG(init = list(check_beta+1e-2, rep(max(check_beta), m+1), check_beta*2), ep = ep, 
                            tau = new_tau, X = X, Y = Y, h = h, nbatch = nbatch, 
                            scale = scale, check_beta = check_beta, Cx = Cx, nchains = 3,
                            check_data = check_data, ub = ub , lb = lb, method = method,
                            type = "upper", lower_accept = lower_accept, upper_accept = upper_accept,
                            update_after = update_after, adjust_scale_by = adjust_scale_by)
            if (is.null(out)) {
                return(NULL)
            }
            proposed_beta = tail(out[[1]][[1]], 1)
            ans = cbind(ans, t(proposed_beta))
            scale_output = cbind(scale_output, tail(out[[3]][[1]], 1))
            accept_rate = cbind(accept_rate, out[[2]][[1]])
        }
    }
    
    tau = sort(tau)
    colnames(ans) = tau
    colnames(scale_output) = tau
    rownames(scale_output) = "Scale"
    colnames(accept_rate) = tau
    
    return(list(ans, scale_output, accept_rate))
}


## Synthesizing employment with KNG and smooth quantile regression loss
# The linear constraints (upper bound on employment based on the public data, 
# lower bound on employment  = 0) are set to improve the utility of synthetic data.
# (This version does not perform log transformation)

# Initial input: ./output/[naics4]/[naics4]_syn_data_year_mu.csv
# Final output: ./output/[naics4]/[naics4]_syndata_stepKNG_[kernel].csv
# Intermediate input and output will have similar syntax, but with an extra
# year parameter.

# Last updated: 4/23/24

# source code for KNG


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

generate_syn_emp = function(raw_data, syn_data, syn_var, taus, h, epsilon, Cx, kernel, mod) {
    curr_n = nrow(syn_data)
    
    # perform dpqr modeling
    # tmp = dpqr_KNG(data = raw_data, syn_var = syn_var, taus = taus, h = h, 
    #                epsilon = epsilon, Cx = Cx, kernel = kernel, mod = mod)
    raw_data = data_jitter(raw_data)
    scale = ifelse(mod == 'continuer' & syn_var != 'emp_1977', 1e-1, 1e-1)
    tmp = stepwiseKNG(data = raw_data, total_eps = epsilon, tau = taus, scale = scale, #1e-1,
                      Cx = Cx, nbatch = 100000, method = 'varying', h = h, lb = 0, 
                      ub = Cx + 1000, check_data = syn_data, lower_accept = 0.1, upper_accept = 0.7,
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
    de_id = function(x) {
        sep_loc = gregexpr('_', x)[[1]]
        lbdnum_id = as.numeric(substr(x, 1, sep_loc[1] -1))
        fye = as.numeric(substr(x, sep_loc[1] + 1, sep_loc[2] -1))
        lye = as.numeric(substr(x, sep_loc[2] + 1, nchar(x)))
        return(c(lbdnum_id, fye, lye))
    }
    ans = cbind(ans, t(sapply(syn_data_id, de_id)))
    colnames(ans) = c(syn_var, 'lbdnum', 'firstyear_emp', 'lastyear_emp')
    
    return(ans)
    
}

syn_data_all_years = function(seed) {
    start_year = ifelse(seed == 2 | seed == 4, 1980, 1981)
    for (curr_year in start_year:2021) { # rmb to change back to 1977
        # set variable to synthesize and predictors
        syn_var = paste0('emp_', curr_year)
        print(curr_year)
        print(Sys.time())
        
        mu_var = paste0('mu_', mu_ref$mu_year[mu_ref$year == curr_year])
        prev_emp = paste0('emp_', curr_year-1)
        if (curr_year == 1977) {
            pred_var_cont = c('year_to_death', mu_var)
        } else {
            pred_var_cont = c('year_to_death', mu_var, prev_emp) #age, ,
        }    
        pred_var_birth = c('year_to_death', mu_var)
        
        # read in previously generated synthetic data
        # if year = 1977, use the general synthetic year and multiunit file (used across different DP mechs)
        # else use the specific DP mech synthetic dataset
        if (curr_year == 1977) { # remember to change to _seed, seed
            syn_data = read.csv(paste0("./output/", naics4, "/", naics4, "_syn_data_year_mu_seed", seed, ".csv"))
        } else {
            syn_data = read.csv(paste0("./output/", naics4, "/", naics4, "_syndata_stepKNG_", kernel, "_", curr_year-1, "_seed", seed, ".csv")) # change method
        }
        syn_data = left_join(syn_data, st_ref, by = 'st')
        syn_data$region[is.na(syn_data$region)] = 99
        
        ## split raw and synthetic data into birth and continuer establishments for modeling
        
        # Within continuer establishments, there are two cases
        # 1. Have valid previous year employment (lag1)
        # => We model current year employment based on previous year employment
        # 2. Do not have valid previous year employment (lag1)
        # => For now treat it as birth model (i.e., model without predictors)
        # => Later, if have time, merge with data from lag2-employment.
        # If use lag2-employment if it exists, else treat like birth model.
        # Note: This solution can only starts from 1979 (first year with lag2 info)
        
        raw_data_cont = raw_data[raw_data$firstyear_emp < curr_year & 
                                     raw_data$lastyear_emp >= curr_year &
                                     !is.na(raw_data[, syn_var]), ]
        raw_data_birth = raw_data[raw_data$firstyear_emp == curr_year & 
                                      raw_data$lastyear_emp >= curr_year &
                                      !is.na(raw_data[, syn_var]), ]
        
        # split synthetic data into birth, continuer, and NA for modeling
        syn_data_cont = syn_data[syn_data$firstyear_emp < curr_year & 
                                     syn_data$lastyear_emp >= curr_year, ]
        syn_data_birth = syn_data[syn_data$firstyear_emp == curr_year & 
                                      syn_data$lastyear_emp >= curr_year, ]
        syn_data_na = syn_data[syn_data$firstyear_emp > curr_year | 
                                   syn_data$lastyear_emp < curr_year, ]
        
        if(curr_year > 1977) {
            if (curr_year == 1978 & sum(is.na(raw_data_cont[, prev_emp])) > 0) {
                raw_data_cont_tmp = raw_data_cont[is.na(raw_data_cont[, prev_emp]), ]
                raw_data_birth = rbind(raw_data_birth, raw_data_cont_tmp)
                raw_data_cont = raw_data_cont[!rownames(raw_data_cont) %in% rownames(raw_data_cont_tmp), ]
                
                # syn_data_birth = rbind(syn_data_birth, syn_data_cont[rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ])
                # syn_data_cont = syn_data_cont[!rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ]
                
            } else if (curr_year > 1978 & sum(is.na(raw_data_cont[, prev_emp])) > 0) {
                prev_prev_emp = paste0('emp_', curr_year-2)
                raw_data_cont_tmp = raw_data_cont[is.na(raw_data_cont[, prev_emp]), ]
                raw_data_cont = raw_data_cont[!rownames(raw_data_cont) %in% rownames(raw_data_cont_tmp), ]
                # syn_data_cont_tmp = syn_data_cont[rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ]
                # syn_data_cont = syn_data_cont[!rownames(syn_data_cont) %in% rownames(raw_data_cont_tmp), ]
                
                tmp = raw_data_cont_tmp[is.na(raw_data_cont_tmp[, prev_prev_emp]), ]
                raw_data_birth = rbind(raw_data_birth, tmp)
                # syn_data_birth = rbind(syn_data_birth, syn_data_cont_tmp[rownames(syn_data_cont_tmp) %in% rownames(tmp), ])
                
                raw_data_cont_tmp = raw_data_cont_tmp[!rownames(raw_data_cont_tmp) %in% rownames(tmp), ]
                raw_data_cont_tmp[, prev_emp] = raw_data_cont_tmp[, prev_prev_emp]
                raw_data_cont = rbind(raw_data_cont, raw_data_cont_tmp)
                
                # syn_data_cont_tmp = syn_data_cont_tmp[!rownames(syn_data_cont_tmp) %in% rownames(tmp), ]
                # syn_data_cont_tmp[, prev_emp] = syn_data_cont_tmp[, prev_prev_emp]
                # syn_data_cont = rbind(syn_data_cont, syn_data_cont_tmp)        
                
            }
            
        }
        
        
        ## assign NA values for establishments that do not exist in the curr_year
        ans_na = cbind(NA, syn_data_na$lbdnum, syn_data_na$firstyear_emp, syn_data_na$lastyear_emp)
        ans_na = as.data.frame(ans_na)
        colnames(ans_na) = c(syn_var, 'lbdnum', 'firstyear_emp', 'lastyear_emp')
        
        # save as rownames to retrieve it easily later without impacting matrix multiplication
        # lbdnum will be used later on to merge to the main synthetic dataset
        row.names(syn_data_cont) = paste(syn_data_cont$lbdnum, syn_data_cont$firstyear_emp, syn_data_cont$lastyear_emp, sep = "_")
        row.names(syn_data_birth) = paste(syn_data_birth$lbdnum, syn_data_birth$firstyear_emp, syn_data_birth$lastyear_emp, sep = '_')
        
        
        ## model and generate synthetic employment for establishments that were born 
        ## on the curr_year (birthers) by regions
        # emp_birthyear ~ most_current_mu + year_to_death
        
        # create year_to_death variable from lastyear_emp and curr_year for birth model
        syn_data_birth$year_to_death = syn_data_birth$lastyear_emp - curr_year
        raw_data_birth$year_to_death = raw_data_birth$lastyear_emp - curr_year
        
        # model by region
        ans_birth = NULL
        region_list = sort(unique(raw_data_birth$region))
        for(r in region_list) {
            print(paste('Year', curr_year, '- Birth model region', r))
            raw_data_birth_tmp = raw_data_birth[raw_data_birth$region == r, c(pred_var_birth, syn_var)]
            syn_data_birth_tmp = syn_data_birth[syn_data_birth$region == r, pred_var_birth]
            ans_birth_tmp = NULL
            if (nrow (raw_data_birth_tmp) > 0 & nrow(syn_data_birth_tmp) > 0) {
                while(is.null(ans_birth_tmp)) {
                    print('1')
                    #tic()
                    ans_birth_tmp = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_birth_tmp, 
                                                                           syn_data = syn_data_birth_tmp, syn_var = syn_var, 
                                                                           taus = taus, h = h, epsilon = epsilon, Cx = 46, kernel = kernel, 
                                                                           mod = 'birth')}, timeout=600), error = function(cond) NULL)
                    #toc()
                    
                }
            } else if(nrow  (raw_data_birth_tmp) == 0 & nrow(syn_data_birth_tmp) > 0) {
                ans_birth_tmp = cbind(0, syn_data_birth_tmp$lbdnum)
                conames(ans_birth_tmp) = c(syn_var, 'lbdnum')
            }
            ans_birth = rbind(ans_birth, ans_birth_tmp)
            
        }
        
        
        ## model and generate synthetic employment for establishments that were born 
        ## before the curr_year (continuers) by regions
        # emp_contyear ~ most_current_mu + age + year_to_death + prev_emp
        
        # create year_to_death variable from lastyear_emp and curr_year for birth model
        # syn_data_cont$age = curr_year - syn_data_cont$firstyear_emp
        # raw_data_cont$age = curr_year - raw_data_cont$firstyear_emp
        syn_data_cont$year_to_death = syn_data_cont$lastyear_emp - curr_year
        raw_data_cont$year_to_death = raw_data_cont$lastyear_emp - curr_year
        
        # model by region
        ans_cont = NULL
        region_list = sort(unique(raw_data_cont$region))
        Cx_tmp = ifelse(syn_var == 'emp_1977', 46, max(maxemp_ref$maxemp))
        for(r in region_list) {
            print(paste('Year', curr_year, '- Continuer model region', r))
            raw_data_cont_tmp = raw_data_cont[raw_data_cont$region == r, c(pred_var_cont, syn_var)]
            syn_data_cont_tmp = syn_data_cont[syn_data_cont$region == r, pred_var_cont]
            ans_cont_tmp = NULL
            if (nrow (raw_data_cont_tmp) > 0 & nrow(syn_data_cont_tmp) > 0) {
                if (syn_var == 'emp_1977') {
                    while(is.null(ans_cont_tmp)) {
                        tic()
                        ans_cont_tmp = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_cont_tmp, 
                                                                              syn_data = syn_data_cont_tmp, syn_var = syn_var, 
                                                                              taus = taus, h = h, epsilon = epsilon, Cx = max(Cx_tmp, 46), 
                                                                              kernel = kernel, mod = 'continuer')}, timeout= 1200), 
                                                error = function(cond) NULL)
                        toc()
                        
                    }
                    
                } else {
                    raw_data_cont_tmp0 = raw_data_cont_tmp[raw_data_cont_tmp[,2] == 0, ]
                    raw_data_cont_tmp1 = raw_data_cont_tmp[raw_data_cont_tmp[,2] == 1, ]
                    syn_data_cont_tmp0 = syn_data_cont_tmp[syn_data_cont_tmp[,2] == 0, ]
                    syn_data_cont_tmp1 = syn_data_cont_tmp[syn_data_cont_tmp[,2] == 1, ]
                    
                    raw_data_cont_tmp0 = raw_data_cont_tmp0[, -2]
                    raw_data_cont_tmp1 = raw_data_cont_tmp1[, -2]
                    syn_data_cont_tmp0 = syn_data_cont_tmp0[, -2]
                    syn_data_cont_tmp1 = syn_data_cont_tmp1[, -2]
                    print(paste('Year', curr_year, '- Continuer model region', r, '- Mu = 0'))
                    ans_cont_tmp0 = NULL
                    while(is.null(ans_cont_tmp0)) {
                        tic()
                        ans_cont_tmp0 = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_cont_tmp0, 
                                                                               syn_data = syn_data_cont_tmp0, syn_var = syn_var, 
                                                                               taus = taus, h = h, epsilon = epsilon, Cx = max(Cx_tmp, 46), 
                                                                               kernel = kernel, mod = 'continuer')}, timeout= 1200), 
                                                 error = function(cond) NULL)
                        toc()
                        
                    }
                    print(paste('Year', curr_year, '- Continuer model region', r, '- Mu = 1'))
                    ans_cont_tmp1 = NULL
                    while(is.null(ans_cont_tmp1)) {
                        tic()
                        ans_cont_tmp1 = tryCatch(withTimeout({generate_syn_emp(raw_data = raw_data_cont_tmp1, 
                                                                               syn_data = syn_data_cont_tmp1, syn_var = syn_var, 
                                                                               taus = taus, h = h, epsilon = epsilon, Cx = max(Cx_tmp, 46), 
                                                                               kernel = kernel, mod = 'continuer')}, timeout= 900), 
                                                 error = function(cond) NULL)
                        toc()
                        
                    }
                    ans_cont_tmp = rbind(ans_cont_tmp0, ans_cont_tmp1)
                    
                }
                
            } else if (nrow (raw_data_cont_tmp) == 0 & nrow(syn_data_cont_tmp) > 0) {
                ans_cont_tmp = cbind(0, syn_data_cont_tmp$lbd_num, syn_data_cont_tmp$firstyear_emp, syn_data_cont_tmp$lastyear_emp)
                colnames(ans_cont_tmp) = c(syn_var, 'lbdnum', 'firstyear_emp', 'lastyear_emp')
                
            }
            ans_cont = rbind(ans_cont, ans_cont_tmp)
            
        }
        
        ans_all = rbind(ans_birth, ans_cont, ans_na)
        ans_all$lbdnum = as.numeric(ans_all$lbdnum)
        syn_data = left_join(syn_data, ans_all, by = c('lbdnum', 'firstyear_emp', 'lastyear_emp'))
        syn_data = syn_data[, -which(colnames(syn_data) == 'region')]
        
        if (curr_year == 2021) {
            output_file = paste0("./output/", naics4, "/", naics4, "_syndata_stepKNG_", kernel, "_seed", seed, ".csv")
        } else {
            output_file = paste0("./output/", naics4, "/", naics4, "_syndata_stepKNG_", kernel, "_", curr_year, "_seed", seed, ".csv")
        }
        
        write.csv(syn_data, output_file, row.names = FALSE)
        
        old_output_file = paste0("./output/", naics4, "/", naics4, "_syndata_stepKNG_", kernel, "_", curr_year-1, "_seed", seed, ".csv")
        if (file.exists(old_output_file)) {
            file.remove(old_output_file)
        }
        
    }
    
}

library(dplyr)
library(caret)
library(quantreg)

naics4 = 2379 #as.numeric(commandArgs(trailingOnly = TRUE)) #2371 2379 #
naics3 = floor(naics4/10)

data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', naics3, '.csv'))
data = data[floor(data$bds_vcnaics/100) == naics4, ]
n = nrow(data)
rownames(data) = c(1:n)

st_ref = read.csv('st_ref.csv')
st_ref = st_ref[, c('st', 'region')]
maxemp_ref = read.csv('maxemp_ref_regionnaics2.csv')
maxemp_ref = maxemp_ref[maxemp_ref$naics2 == floor(naics4/100), ]
mu_ref = data.frame(year = c(1977:2021), mu_year = c(1977, rep(seq(1982, 2017, 5), each = 5), rep(2022, 4)))

# rearrange columns
data = data[, c('lbdnum', 'st', 'firstyear_emp', 'lastyear_emp', paste('mu5', seq(1977, 2022, 5), sep = '_'), 
                paste('emp', c(1977:2021), sep = '_'))]
colnames(data) = c('lbdnum', 'st', 'firstyear_emp', 'lastyear_emp', paste('mu', seq(1977, 2022, 5), sep = '_'), 
                   paste('emp', c(1977:2021), sep = '_'))
taus = c(seq(0.1, 0.9, 0.1), 0.99)
h = 1
#h = 1.06*min(sse, iqr/1.34)*n^(-1/5)
kernel = 'logistic'

raw_data = data
raw_data = left_join(raw_data, st_ref, by = 'st')
raw_data$region[is.na(raw_data$region)] = 99

epsilon = 0.25


# library(doParallel)
library(tictoc)
# num_cores= 5 #detectCores()-1 #use all available core but 1
# print(num_cores)
# workers=makeCluster(num_cores,type="SOCK")
# registerDoParallel(workers)

#tic()
# Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
# load libraries
source('./code/source_functions_kng.R')
library(dplyr)
library(caret)
library(quantreg)
library(tictoc)
library(R.utils)

seed = as.integer(commandArgs(trailingOnly = TRUE)) #2371 2379 #
#seed = 2
set.seed(seed)
syn_data_all_years(seed)
# }
#toc()
# stopCluster(workers)

