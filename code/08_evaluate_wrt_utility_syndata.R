# load packages
library(dplyr)
library(synthpop)
library(transport)


wasserstein_randomization_test <- function(a, b, n_rep = 1000){
  require(transport)
  rand_dist <- NULL
  if(class(a) == "factor"){
    sample_dist <- transport::wasserstein1d(table(a), table(b))
  } else {
    sample_dist <- transport::wasserstein1d(a, b)
  }
  for(i in 1:n_rep){  
    tmp <- c(a, b)
    sel <- sample(1:length(tmp), length(a))
    tmp_a <- tmp[sel]
    tmp_b <- tmp[-sel]
    if(class(a) == "factor"){
      w_dist <- transport::wasserstein1d(table(tmp_a), table(tmp_b))
    } else {
      w_dist <- transport::wasserstein1d(tmp_a, tmp_b)
    }
    
    rand_dist <- c(rand_dist, w_dist)
  }
  
  p <- sum(sample_dist < rand_dist)/n_rep
  
  return(list(dist = rand_dist, sample_distance = sample_dist, p = p))
}


sic3 = 542 #as.numeric(commandArgs(trailingOnly = TRUE))

var_list = c(paste('emp', c(1977:2000), sep = '_'))
methods = c('opmnc2', 'stepkng1', 'swkng')

raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', floor(sic3/10), '.csv'))
raw_data = raw_data[raw_data$sic3 == sic3, ]
raw_data = raw_data[, var_list]

for (m in 1:length(methods)) { 
  method = methods[m]
  print(method)
  for(seed in 1:5) {
    print(seed)
    set.seed(seed)
    tab = data.frame(sic3 = rep(sic3, 1),
                     seed = rep(seed, 1),
                     wrt = rep(NA, 1))

    syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_", method, "_logistic_seed", seed, ".csv"))
    syn_data = syn_data[,-which(colnames(syn_data) %in% c('lbdnum', 'sic3', c('firstyear', 'lastyear', paste('mu', seq(1977, 2002, 5), sep = '_'))))]
    tmp = NULL
    for (var in colnames(syn_data)) {
      print(var)
      syn_var_tmp = na.omit(syn_data[, var])
      raw_var_tmp = na.omit(raw_data[, var])
      wrt = wasserstein_randomization_test(raw_var_tmp, syn_var_tmp, n_rep = 1000)
      tmp = c(tmp, wrt$sample_distance/median(wrt$dist))
    }
    tab$wrt = mean(tmp)
    output = paste0("./output/general_utility/wrt_utility_syndata_", method, "_emp.csv")
    if (file.exists(output)) {
      tab_og = read.csv(output)
      tab = rbind(tab_og, tab)
    }
    write.csv(tab, output, row.names = FALSE)
  }
}
