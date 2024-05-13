# last updated: 4/22/24
# note: currently not parallelize the utility code as we don't want each seed version
# to open the output file simultaneously

# load packages
library(dplyr)
library(synthpop)
library(tictoc)

sic3 = 542 #as.numeric(commandArgs(trailingOnly = TRUE))

var_list = c('firstyear', 'lastyear', paste('mu5', seq(1977, 2002, 5), sep = '_'), 
             paste('emp', c(1977:2000), sep = '_'))

tic()
for(seed in 1:5) {
    # load libraries
    library(dplyr)
    library(synthpop)
    print(seed)
    set.seed(seed)
    tab = data.frame(sic3 = rep(NA, 1),
                     seed = rep(seed, 1),
                     pMSE_logit = rep(NA, 1),
                     pMSE_CART = rep(NA, 1))
    tab$sic3[1] = sic3
    sic2 = floor(sic3/10)
    syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_opmnc2_logistic_seed", seed, ".csv"))
    syn_data = syn_data[,-which(colnames(syn_data) %in% c('lbdnum', 'sic3'))]
    
    raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
    raw_data = raw_data[raw_data$sic3 == sic3, ]
    raw_data = raw_data[, var_list]
    raw_data$firstyear[raw_data$firstyear < 1977] = 1977
    raw_data$lastyear[raw_data$lastyear > 2000] = 2000
    colnames(raw_data) = colnames(syn_data)
    tab$pMSE_logit[1] = utility.gen(syn_data, raw_data, method = 'logit', maxorder = 0)$pMSE
    # currently not running logit with interaction due to the very high number of possible interactions
    #tab$pMSE_logit_inter[i] = utility.gen(syn_data, raw_data, not.synthesised = c('st'), method = 'logit')$pMSE
    tab$pMSE_CART[1] = utility.gen(syn_data, raw_data, method = "cart", nperms = 20)$pMSE
    output = paste0("./output/general_utility/pMSE_utility_syndata_opmnc2_updated.csv")
    if (file.exists(output)) {
        tab_og = read.csv(output)
        tab = rbind(tab_og, tab)
    }
    write.csv(tab, output, row.names = FALSE)
}
toc()



