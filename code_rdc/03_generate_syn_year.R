# Last updated: 4/22/24
# note added parallel code

# load libraries
library(dplyr)



generate_syn_year = function(sic3, seed) {
    set.seed(sic3)
    sic2 = floor(sic3/10)
    
    # read in raw data and noisy year probability
    data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
    year_probs_data = read.csv(paste0('./output/year_noisy_probs/year_noisy_probs_laplace_', sic2, "_seed", seed, '.csv'))
    year_probs_data$category = c(1:nrow(year_probs_data))
    
    # subset data
    data_sic3 = data[data$sic3 == sic3, ]
    n = nrow(data_sic3)
    
    # generate synthetic year using multinomial distribution
    ans = t(rmultinom(n, 1, year_probs_data$prob_noisy))
    colnames(ans) = year_probs_data$category
    ans = colnames(ans)[max.col(ans)]
    ans = data.frame(category = as.numeric(ans))
    ans = dplyr::left_join(ans, year_probs_data, by = "category")
    ans = ans[, c('firstyear', 'lastyear')]
    
    # prepare output
    # copy over state data as we do not synthesize state values
    ans = cbind(data_sic3$lbdnum, data_sic3$sic3, ans)
    colnames(ans) = c('lbdnum', 'sic3', 'firstyear', 'lastyear')
    
    #return(ans)
    write.csv(ans, paste0("./output/", sic3, "/", sic3, "_syn_data_year_seed", seed, ".csv"), row.names = FALSE)
}

sic3_list = c(239, 542, 829)

# naics2 = 23 # parameter to change
# sic3_data = read.csv('../../data/sic3_inds/sic3_inds.csv', header = FALSE)
# colnames(sic3_data) = 'sic3'
# sic3_list = unique(sic3_data$sic3[floor(sic3_data$sic3/100) == naics2])
# sic3_list = c(sic3_list[c(3:6)], sic3_list[c(1:2, 7:10)])
# 
# library(doParallel)
# library(tictoc)
# num_cores= 5 #detectCores()-1 #use all available core but 1
# print(num_cores)
# workers=makeCluster(num_cores,type="SOCK")
# registerDoParallel(workers)
# 
# tic()
#Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
for (seed in 1:5) {
    library(dplyr)
    print(seed)
    set.seed(seed)
    for (sic3 in sic3_list) {
        print(sic3)
        output_folder = paste0('./output/', sic3)
        if (!dir.exists(output_folder)) {
            dir.create(output_folder)
        }
        
        generate_syn_year(sic3, seed)
        
    }
}


# toc()
# stopCluster(workers)























