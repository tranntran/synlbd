# Last updated: 3/12/24

# load library
library(dplyr)

# first mu status is recorded in 1977, hense the name of the variable
generate_syn_first_mu = function(sic3, year, syn_data, noise_dist, epsilon) {
    ## generate noisy probability of multiunit status using multinomial distribution
    ## and the specified noise dist
    set.seed(sic3)
    
    sic2 = floor(sic3/10)
    data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv')) #remember to change back file name
    data_sic3 = data[data$sic3 == sic3, ]
    data_sic3$lastyear[data_sic3$lastyear > 2000] = 2000
    data_sic3 = data_sic3[, c('firstyear', 'lastyear', paste0('mu5_', year))]
    colnames(data_sic3) = c('firstyear', 'lastyear', 'curr_mu')
    data_sic3 = data_sic3[!is.na(data_sic3$curr_mu), ]
    data_sic3$numyear = data_sic3$lastyear - data_sic3$firstyear
    tmp = data_sic3 %>% group_by(numyear, curr_mu) %>% summarize(freq = n())
    
    all_st_numyear = expand.grid(c(0:25), c(0, 1))
    colnames(all_st_numyear) = c('numyear', 'curr_mu')
    
    all_st_numyear = left_join(all_st_numyear, tmp, by = c('numyear', 'curr_mu'))
    all_st_numyear$freq[is.na(all_st_numyear$freq)] = 0
    
    # add laplace or geom noise to each count
    if (noise_dist == 'laplace') {
        all_st_numyear$freq_noisy = all_st_numyear$freq + rmutil::rlaplace(nrow(all_st_numyear), 0, 1/epsilon)
    } else if (noise_dist == 'geom') {
        p = 1 - exp(-epsilon)
        all_st_numyear$freq_noisy = all_st_numyear$freq + (rgeom(nrow(all_st_numyear), p) - rgeom(nrow(all_st_numyear), p))
    } else if (noise_dist == 'none') {
        all_st_numyear$freq_noisy = all_st_numyear$freq
        
    }
    
    all_st_numyear$freq_noisy[all_st_numyear$freq_noisy < 0] = 0
    
    # calculate the probability of mu given state and number of year with reported emp
    all_st_numyear = group_by(all_st_numyear, numyear) %>% mutate(prob_noisy = freq_noisy/sum(freq_noisy))
    all_st_numyear$prob_noisy[is.na(all_st_numyear$prob_noisy)] = 0
    
    # probability of mu given old mu only (without numyear_emp information)
    # this is to account for the fact that data is too sparse (i.e., no positive probability within numyear)
    all_st_numyear_tmp = group_by(all_st_numyear, curr_mu) %>% summarize(freq = sum(freq), freq_noisy = sum(freq_noisy))
    all_st_numyear_tmp = group_by(all_st_numyear_tmp) %>% mutate(prob_noisy = freq_noisy/sum(freq_noisy))
    all_st_numyear_tmp = cbind(99, all_st_numyear_tmp)
    colnames(all_st_numyear_tmp)[1] = 'numyear'
    all_st_numyear = rbind(all_st_numyear, all_st_numyear_tmp)
    
    # rearrange noisy probability data
    all_st_numyear = all_st_numyear[, c('numyear', 'curr_mu', 'prob_noisy')]
    all_st_numyear = arrange(all_st_numyear, numyear, curr_mu)
    
    
    ## synthesize mu based on synthetic first and last year
    syn_data$numyear = syn_data$lastyear - syn_data$firstyear
    syn_data$curr_mu = NA
    for (i in 1:nrow(syn_data)) {
        curr_st = syn_data$st[i]
        curr_numyear = syn_data$numyear[i]
        tmp = all_st_numyear[all_st_numyear$numyear == curr_numyear,]
        
        # if there is no positive probability in a combination of st and numyear, then use the
        # probability given numyear only (coarsen up the geo unit)
        if(sum(tmp$prob_noisy) == 0) {
            tmp = all_st_numyear[all_st_numyear$numyear == 99,]
        }
        
        syn_data$curr_mu[i] = c(0, 1)[which(t(rmultinom(1, 1, tmp$prob_noisy)) == 1)]

    }
    
    syn_data = syn_data[, c('lbdnum', 'firstyear', 'lastyear', 'curr_mu')]
    colnames(syn_data) = c('lbdnum', 'firstyear', 'lastyear', paste0('mu_', year))
    return(syn_data)
    
}


generate_syn_other_mu = function(sic3, year, syn_data, noise_dist, epsilon) {
    ## generate noisy probability of multiunit status using multinomial distribution
    ## and the specified noise dist
    set.seed(sic3)
    sic2 = floor(sic3/10)
    data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv')) #remember to change back file name
    data_sic3 = data[data$sic3 == sic3, ]
    data_sic3$lastyear[data_sic3$lastyear > 2000] = 2000
    
    data_sic3 = data_sic3[, c('firstyear', 'lastyear', paste0('mu5_', year-5), paste0('mu5_', year))]
    colnames(data_sic3) = c('firstyear', 'lastyear', 'old_mu', 'new_mu')
    data_sic3$numyear = data_sic3$lastyear - data_sic3$firstyear
    tmp = data_sic3 %>% group_by(numyear, old_mu, new_mu) %>% summarize(freq = n())
    
    all_st_numyear = expand.grid(c(0:25), c(0, 1), c(0:1))
    colnames(all_st_numyear) = c('numyear', 'old_mu', 'new_mu')
    
    all_st_numyear = left_join(all_st_numyear, tmp, by = c('numyear', 'old_mu', 'new_mu'))
    all_st_numyear$freq[is.na(all_st_numyear$freq)] = 0
    
    # add laplace or geom noise to each count
    if (noise_dist == 'laplace') {
        all_st_numyear$freq_noisy = all_st_numyear$freq + rmutil::rlaplace(nrow(all_st_numyear), 0, 1/epsilon)
    } else if (noise_dist == 'geom') {
        p = 1 - exp(-epsilon)
        all_st_numyear$freq_noisy = all_st_numyear$freq + (rgeom(nrow(all_st_numyear), p) - rgeom(nrow(all_st_numyear), p))
    } else if (noise_dist == 'none') {
        all_st_numyear$freq_noisy = all_st_numyear$freq
        
    }
    
    all_st_numyear$freq_noisy[all_st_numyear$freq_noisy < 0] = 0
    
    # calculate the probability of mu given state and number of year with reported emp and old mu
    all_st_numyear = group_by(all_st_numyear, numyear, old_mu) %>% mutate(prob_noisy = freq_noisy/sum(freq_noisy))
    all_st_numyear$prob_noisy[is.na(all_st_numyear$prob_noisy)] = 0
    
    # probability of mu given numyear_emp only (without state information)
    # this is to account for the fact that data is too sparse (i.e., no positive probability within and numyear)
    all_st_numyear_tmp = group_by(all_st_numyear, old_mu, new_mu) %>% summarize(freq = sum(freq), freq_noisy = sum(freq_noisy))
    all_st_numyear_tmp = group_by(all_st_numyear_tmp, old_mu) %>% mutate(prob_noisy = freq_noisy/sum(freq_noisy))
    all_st_numyear_tmp = cbind(99, all_st_numyear_tmp)
    colnames(all_st_numyear_tmp)[1] = 'numyear'
    all_st_numyear = rbind(all_st_numyear, all_st_numyear_tmp)
    
    # rearrange noisy probability data
    all_st_numyear = all_st_numyear[, c('numyear', 'old_mu', 'new_mu', 'prob_noisy')]
    all_st_numyear = arrange(all_st_numyear, numyear, old_mu, new_mu)
    
    
    ## synthesize mu based on synthetic first and last year
    syn_data = syn_data[, c('lbdnum', 'firstyear', 'lastyear', paste0('mu_', year-5))]
    colnames(syn_data)[ncol(syn_data)] = 'old_mu' 
    syn_data$numyear = syn_data$lastyear - syn_data$firstyear
    syn_data$new_mu = NA
    for (i in 1:nrow(syn_data)) {
        curr_st = syn_data$st[i]
        curr_numyear = syn_data$numyear[i]
        curr_oldmu = syn_data$old_mu[i]
        tmp = all_st_numyear[all_st_numyear$numyear == curr_numyear &
                                 all_st_numyear$old_mu == curr_oldmu,]

        # if there is no positive probability in a combination of st and numyear, then use the
        # probability given numyear only (coarsen up the geo unit)
        if(sum(tmp$prob_noisy) == 0) {
            tmp = all_st_numyear[all_st_numyear$numyear == 99 &
                                     all_st_numyear$old_mu == curr_oldmu,]
        }
        
        syn_data$new_mu[i] = c(0, 1)[which(t(rmultinom(1, 1, tmp$prob_noisy)) == 1)]
        
    }
    
    syn_data = syn_data[, c('lbdnum', 'firstyear', 'lastyear', 'new_mu')]
    colnames(syn_data) = c('lbdnum', 'firstyear', 'lastyear', paste0('mu_', year))
    return(syn_data)
}


# sic3 = 237110 # change parameter
noise_dist = 'laplace' # change parameter
epsilon = 0.1  # change parameter

# naics2 = 23 # parameter to change
# sic3_data = read.csv('../../data/sic3_inds/sic3_inds.csv', header = FALSE)
# colnames(sic3_data) = 'sic3'
# sic3_list = unique(sic3_data$sic3[floor(sic3_data$sic3/100) == naics2])
# sic3_list = sic3_list[c(1:2, 7:10)] #c(sic3_list[c(3:6)])#, sic3_list[c(1:2, 7:10)])

year_list = seq(1977, 2002, 5)
sic3_list = 239 #c(178, 239, 542, 829)

library(doParallel)
library(tictoc)
num_cores= 5 #detectCores()-1 #use all available core but 1
print(num_cores)
workers=makeCluster(num_cores,type="SOCK")
registerDoParallel(workers)

tic()

#Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
    for (seed in 1:5) {
    print(seed)
    # load library
    library(dplyr)
    set.seed(seed)
    for (sic3 in sic3_list) {
        input_file = paste0("./output/", sic3, "/", sic3, "_syn_data_year_seed", seed, ".csv")
        if (!file.exists(input_file)) {
            print(paste(sic3, 'syn year file does not exist, run step 3 first and retry.'))
        } else {
            for (year in year_list) {
                print(paste(sic3, '-', year))
                if (year == 1977) {
                    syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syn_data_year_seed", seed, ".csv"))
                    syn_data_na = syn_data[syn_data$firstyear > 1977, ]
                    syn_data_na$mu_1977 = NA
                    syn_data_na = syn_data_na[, c('lbdnum', 'firstyear', 'lastyear', 'mu_1977')]
                    
                    syn_data_mod = syn_data[syn_data$firstyear <= 1977, ]
                    syn_data_mod = generate_syn_first_mu(sic3, year, syn_data_mod, noise_dist, epsilon)
                    
                    syn_data_new = rbind(syn_data_na, syn_data_mod)
                    syn_data = left_join(syn_data, syn_data_new, by = c('lbdnum', 'firstyear', 'lastyear'))
                    
                    write.csv(syn_data, paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_1977_seed", seed, ".csv"), row.names = FALSE)
                } else {
                    mu_var = paste0('mu_', year)
                    syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_", year-5, "_seed", seed, ".csv"))
                    if (sum(syn_data$firstyear > year | syn_data$lastyear <= year - 5) > 0) {
                        syn_data_na = syn_data[syn_data$firstyear > year | syn_data$lastyear <= year - 5, ]
                        syn_data_na = as.data.frame(cbind(syn_data_na$lbdnum, syn_data_na$firstyear, syn_data_na$lastyear, NA))
                        colnames(syn_data_na) = c('lbdnum', 'firstyear', 'lastyear', mu_var)
                    } else {
                        syn_data_na = NULL
                    }

                    
                    # subset data for establishments that were born between last mu update and current mu update
                    # fit this subset data to first mu model
                    syn_data_first = syn_data[syn_data$firstyear <= year & syn_data$firstyear > year-5, ]
                    if (nrow(syn_data_first) > 0) {
                        syn_data_first = generate_syn_first_mu(sic3, year, syn_data_first, noise_dist, epsilon)
                        
                    } else {
                        syn_data_first = NULL
                    }
                    
                    
                    # subset data for establishments that already went through at least one mu update
                    syn_data_other = syn_data[syn_data$firstyear <= year-5 & syn_data$lastyear > year-5, ]
                    if (nrow(syn_data_other) > 0) {
                        syn_data_other = generate_syn_other_mu(sic3, year, syn_data_other, noise_dist, epsilon)
                        
                    } else {
                        syn_data_other = NULL
                    }
                    
                    
                    
                    syn_data_new = rbind(syn_data_na, syn_data_first, syn_data_other)
                    syn_data = left_join(syn_data, syn_data_new, by = c('lbdnum', 'firstyear', 'lastyear'))
                    
                    if (year == 2002) {
                        output_file = paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_seed", seed, ".csv")
                    } else {
                        output_file = paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_", year, "_seed", seed, ".csv")
                    }
                    write.csv(syn_data, output_file, row.names = FALSE)
                    file.remove(paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_", year-5, "_seed", seed, ".csv"))
                    
                }
            }
            
            
        }
        
    }

}

toc()
stopCluster(workers)

