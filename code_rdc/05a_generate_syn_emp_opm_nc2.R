## Synthesizing employment with OPM and smooth quantile regression loss
# This code optimize the private loss function in OPM subject to linear 
# constraints using constrOptim. The linear constraints (upper bound on 
# employment based on the public data, lower bound on employment  = 0) 
# are set to improve the utility of synthetic data.

# The main difference between this version and OPM_constraints is that it
# does not use linear constraints to ensure the output satisfies >= 0
# requirement but instead use post-processing, hence the name nc (no constraints).

# Initial input: ./output/[sic3]/[sic3]_syn_data_year_mu.csv
# Final output: ./output/[sic3]/[sic3]_syndata_opmnc2_[kernel].csv
# Intermediate input and output will have similar syntax, but with an extra
# year parameter.

# Last updated: 4/27/24

# source code for opm
## function that output kernel specific variable and functions
select_kernel = function(kernel = c('gaussian', 'logistic', 'uniform')) {
    kernel = match.arg(kernel)
    if (kernel == 'gaussian') {
        k_bar = 1/sqrt(2*pi)
        integral_fn = function(u, h) {
            ans = sqrt(2/pi)*exp(-(u/h)^2/2) +(u/h)*(1-2*pnorm(-u/h))
            return(ans)
        }
        cdf_k_h = function(u, h) {
            return(pnorm(u/h))
        }
    } else if (kernel == 'logistic') {
        k_bar = 1/4
        integral_fn = function(u, h) {
            ans = u/h + 2*log(1+ exp(-u/h))
            return(ans)
        }
        cdf_k_h = function(u, h) {
            return(1/(1 + exp(-u/h)))
        }
    } else if (kernel == 'uniform') {
        k_bar = 1/2
        integral_fn = function(u, h) {
            ans = ifelse(abs(u/h) <= 1, u^2/(2*h^2) + 1/2, abs(u/h))
            return(ans)
        }
        cdf_k_h = function(u, h) {
            return(min((u/h + 1)/2, 1)*(u/h >= -1))
        }
    }
    return(list(k_bar, integral_fn, cdf_k_h))
}

## OPM (Kifer et al, 2012) with check loss smoothed by convolution methods
opm_convsm = function(theta, x, y, tau, h, epsilon = 1, Cx, b, 
                      kernel = c('gaussian', 'logistic', 'uniform')) {
    kernel = match.arg(kernel)
    kernel_info = select_kernel(kernel)
    k_bar = kernel_info[[1]]
    integral_fn = kernel_info[[2]]
    
    theta = matrix(theta, ncol = 1)
    n = nrow(x)
    p = ncol(x)
    xi = n*max(tau, 1-tau)*Cx*sqrt(p)
    lambda = n*k_bar*Cx^2/h*p
    ans = 0
    for (i in 1:n) {
        u = y[i] - x[i, ]%*%theta
        ans = ans + h/2*integral_fn(u, h) + (tau -1/2)*u
        # if(is.infinite(ans)) {
        #     print(i)
        #     break
        #     break
        # }
        
    }
    delta = 2*lambda/epsilon
    gamma = n*k_bar/h # parameter for strong convexity
    
    # sample from Gamma distribution to satisfy epsilon-DP
    # norm = rgamma(n = 1, shape = p, rate = epsilon/(2*xi))
    # vector = rnorm(n = p, m = 0, s = 1)
    # b = vector*norm/sqrt(sum(vector^2))
    
    ans = ans + max(delta - gamma, 0)*norm(theta, "F")^2/(2*n) + b%*%theta/n
    
    return(ans)
}

opm_convsm_grad = function(theta, x, y, tau, h, epsilon = 1, Cx, b,  
                           kernel = c('gaussian', 'logistic', 'uniform')) {
    kernel = match.arg(kernel)
    kernel_info = select_kernel(kernel)
    k_bar = kernel_info[[1]]
    cdf_k_h = kernel_info[[3]]
    
    theta = matrix(theta, ncol = 1)
    n = nrow(x)
    p = ncol(x)
    lambda = n*k_bar*Cx^2/h*p
    ans = 0
    for (i in 1:n) {
        u = x[i, ]%*%theta - y[i]
        ans = (cdf_k_h(u, h) - tau) %*% x[i, ]
    }
    delta = 2*lambda/epsilon
    gamma = n*k_bar/h # parameter for strong convexity
    
    ans = ans + max(delta - gamma, 0)*norm(theta, "F")/n + b/n
    
    return(ans)
}

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

# x = data[, -which(colnames(data) == syn_var)]
# x = cbind(1, x)
# colnames(x) = NULL    
# x = as.matrix(x)
# x[x > Cx] = Cx


# tmp = constrOptim(theta = as.matrix(rep(2, ncol(x))), f = opm_convsm, grad = opm_convsm_grad,
#                   ui = ui, ci = ci, x = x, y = y, Cx = Cx, tau = tau, h = h,
#                   b = b, epsilon = eps, kernel = kernel, method = 'BFGS')

dpqr_opm = function(data, syn_var, taus = c(seq(0.1, 0.9, 0.1), 0.99), 
                    h, epsilon = 1, Cx, mod, 
                    kernel = c('gaussian', 'logistic', 'uniform')) {
    data = data_jitter(data)
    y = as.matrix(data[, syn_var])
    x = data[, -which(colnames(data) == syn_var)]
    x = cbind(1, x)
    colnames(x) = NULL    
    x = as.matrix(x)
    x[x > Cx] = Cx
    p = ncol(x)
    n = nrow(x)
    eps = epsilon/length(taus)
    ct = 0
    ans = NULL
    trigger = FALSE
    curr_year = as.numeric(substr(syn_var, 5, 8))
    while ((is.null(ans) | trigger == TRUE) & ct < 25) {
        ans = NULL
        trigger = FALSE
        ct = ct + 1
        for (tau in taus) {
            print(tau)
            tmp = NULL
            attempt = 0
            xi = n*max(tau, 1-tau)*Cx*sqrt(p)
            # sample from Gamma distribution to satisfy epsilon-DP
            norm = rgamma(n = 1, shape = p, rate = eps/(2*xi))
            vector = rnorm(n = p, m = 0, s = 1)
            b = vector*norm/sqrt(sum(vector^2))
            
            while(is.null(tmp) & attempt <= 9) {
                attempt = attempt + 1
                if (attempt <= 2 | ct < 2) {
                    start_points = tryCatch(coef(quantreg::rq(y~ x[,-1], tau = tau)), 
                                            error = function(cond) {matrix(c(rep(2, ncol(x)-1), 1e-3), nrow = ncol(x), ncol = 1)})
                    start_points = data_jitter(start_points)
                } else if (attempt <= 2 | (ct >= 2 & ct <=15)) {
                    start_points = matrix(c(rep(2, ncol(x)-1), 1e-3), nrow = ncol(x), ncol = 1)
                } else {
                    start_points = matrix(c(rep(1, ncol(x)-1), 0), nrow = ncol(x), ncol = 1)
                }
                try(
                    # double check optimal h
                    {tmp = optim(par = as.matrix(start_points), f = opm_convsm, gr = opm_convsm_grad,
                                 x = x, y = y, Cx = Cx, tau = tau, h = h, b = b, epsilon = eps, 
                                 kernel = kernel, method = 'BFGS')$par}, #
                    silent = TRUE
                )
            }
            if (is.null(tmp) & attempt == 10) {
                print(paste(tau, '- Optimization error, trying again.'))
                trigger = TRUE
                break
                break
            }
            
            ans = cbind(ans, tmp)
        }
    }
    colnames(ans) = taus
    
    return(ans)
}

generate_syn_emp = function(raw_data, syn_data, syn_var, taus, h, epsilon, Cx, kernel, mod) {
    curr_n = nrow(syn_data)
    
    # perform dpqr modeling
    tmp = dpqr_opm(data = raw_data, syn_var = syn_var, taus = taus, h = h, 
                   epsilon = epsilon, Cx = Cx, kernel = kernel, mod = mod)
    
    # generate synthetic data
    samp = sample(1:length(taus), curr_n, replace = TRUE, prob = rep(0.1, length(taus)))
    y_hat = as.matrix(cbind(1, syn_data)) %*% tmp
    coord = cbind(c(1:curr_n), samp)
    ans = round(y_hat[coord]) #y_hat[coord] #
    n_neg = sum(ans < 0)
    if (!is.na(n_neg) & n_neg > 0) {
        ans[ans < 0] = rbinom(n_neg, 1, 0.5) # post processing
    }
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
    colnames(ans) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
    
    return(ans)
    
}

syn_data_all_years = function(seed) {
    for (curr_year in 1977:2000) { # rmb to change back to 2021
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
            syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_opmnc2_", kernel, "_", curr_year-1, "_seed", seed, ".csv")) # change method
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
                ans_birth = generate_syn_emp(raw_data = raw_data_birth, syn_data = syn_data_birth, 
                                         syn_var = syn_var, taus = taus, h = h, epsilon = epsilon, 
                                         Cx = 25, kernel = kernel, mod = 'birth') #allow max year to death as 25
            } else if(nrow(raw_data_birth) == 0 & nrow(syn_data_birth) > 0) {
                ## assign NA values for establishments that do not exist in the curr_year
                ans_birth = cbind(NA, syn_data_birth$lbdnum, syn_data_birth$firstyear, syn_data_birth$lastyear)
                ans_birth = as.data.frame(ans_birth)
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
            
            Cx_tmp = 200
            raw_data_cont = raw_data_cont[, c(pred_var_cont, syn_var)]
            syn_data_cont = syn_data_cont[, pred_var_cont]
            if (nrow (raw_data_cont) > 0 & nrow(syn_data_cont) > 0) {
                ans_cont = generate_syn_emp(raw_data = raw_data_cont, syn_data = syn_data_cont, 
                                                syn_var = syn_var, taus = taus, h = h, epsilon = epsilon, 
                                                Cx = max(Cx_tmp, 25), kernel = kernel, mod = 'continuer')
            } else if (nrow (raw_data_cont) == 0 & nrow(syn_data_cont) > 0) {
                ans_cont = cbind(0, syn_data_cont$lbd_num, syn_data_cont$firstyear, syn_data_cont$lastyear)
                colnames(ans_cont) = c(syn_var, 'lbdnum', 'firstyear', 'lastyear')
                
            }

            
        } else {
            print(paste('Year', curr_year, '- Skip continuer model'))
            ans_cont = NULL
        }
        
        
        ans_all = rbind(ans_birth, ans_cont, ans_na)
        ans_all$lbdnum = as.numeric(ans_all$lbdnum)
        syn_data = left_join(syn_data, ans_all, by = c('lbdnum', 'firstyear', 'lastyear'))

        if (curr_year == 2000) {
            output_file = paste0("./output/", sic3, "/", sic3, "_syndata_opmnc2_", kernel, "_seed", seed, ".csv")
        } else {
            output_file = paste0("./output/", sic3, "/", sic3, "_syndata_opmnc2_", kernel, "_", curr_year, "_seed", seed, ".csv")
        }
        
        write.csv(syn_data, output_file, row.names = FALSE)
        
        old_output_file = paste0("./output/", sic3, "/", sic3, "_syndata_opmnc2_", kernel, "_", curr_year-1, "_seed", seed, ".csv")
        if (file.exists(old_output_file)) {
            file.remove(old_output_file)
        }
        
    }
    
}


library(dplyr)
library(caret)
library(quantreg)

sic3 = as.numeric(commandArgs(trailingOnly = TRUE)) #2371
sic2 = floor(sic3/10)

data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
data = data[data$sic3 == sic3, ]
data$lastyear[data$lastyear > 2000] = 2000
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


library(doParallel)
library(tictoc)
num_cores= 5 #detectCores()-1 #use all available core but 1
print(num_cores)
workers=makeCluster(num_cores,type="SOCK")
registerDoParallel(workers)

tic()
Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
    # load libraries
    library(dplyr)
    library(caret)
    library(quantreg)
    
    set.seed(seed)
    syn_data_all_years(seed)
}
toc()
stopCluster(workers)

