library(ggplot2)
library(data.table)
library(dplyr)


sic3 = 542 #as.numeric(commandArgs(trailingOnly = TRUE)) #2373
set.seed(sic3)


k = 3
combn = 100
all_samp = t(combn(24, 3))
all_samp = all_samp[sample(1:nrow(all_samp), combn), ]


methods = c('opmnc2', 'stepkng1', 'swkng')
for (m in 1:3){#length(methods)) {
  method = methods[m]
  nist_score = matrix(NA, nrow = combn, ncol = 6)
  
  for (l in 1:combn) {
    print(l)
    samp = paste0('emp_', 1977:2000)[all_samp[l, ]]
    
    raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', floor(sic3/10), '.csv'))
    raw_data = raw_data[raw_data$sic3 == sic3, ]
    raw_data = raw_data[, samp]
    n = nrow(raw_data)
    bins = apply(raw_data, 2, quantile, na.rm = TRUE)
    vars = colnames(raw_data)
    tab = CJ(x1=c(0:5, NA), x2=c(0:5, NA), x3=c(0:5, NA))
    colnames(tab) = vars
    
    for (j in 1:6) {
      if (j == 6) {
        dt = raw_data 
      } else {
        dt = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_", method, "_logistic_seed", j, ".csv"))
        dt = dt[, samp]
      }
      tmp = NULL
      for (i in 1:k) {
        tmp = cbind(tmp, 
                    findInterval(dt[, i], bins[, i], left.open = T))
      }
      colnames(tmp) = vars
      tmp = as.data.table(tmp)
      tmp_tab = tmp[, .N, by=vars]
      
      new_name = ifelse(j==6, 'raw', paste0('syn', j))
      setnames(tmp_tab, "N", new_name)
      tab = merge(tab, tmp_tab, by = vars, all.x = T)
    }
    tab[is.na(tab)] = 0
    colSums(tab)
    tab = as.data.frame(tab)
    for (j in 1:6) {
      tmp = abs(tab[, j+3] - tab[, 9])/n
      raw_score = sum(tmp)
      nist_score[l, j] = (2-raw_score)/2*1000
    }
    
  }
  
  if(all(nist_score[, 6] == 1000)) {
    nist_score = nist_score[, -6]
    nist_score_tmp = apply(nist_score, 2, mean)
  } else {
    print("Something wrong, double check")
  }
  
  tab = data.frame(sic3 = rep(sic3, 1),
                   km_mean = rep(mean(nist_score_tmp), 1),
                   km_mean = rep(sd(nist_score_tmp)/sqrt(5), 1))
  
  output = paste0("./output/general_utility/km_utility_syndata_", method, "_emp.csv")
  if (file.exists(output)) {
    tab_og = read.csv(output)
    tab = rbind(tab_og, tab)
  }
  write.csv(tab, output, row.names = FALSE)
  
}



# if(all(nist_score[, 6] == 1000)) {
#   nist_score = nist_score[, -6]
#   png(filename = paste0('./output/kmarginals/', sic3, '_kmarginals_swkng.png'), width = 650, height = 450)
#   hist(apply(nist_score, 1, mean), breaks = 30, 
#        main = paste0('Sandwich KNG'),
#        xlab = 'Score')
#   dev.off()
#   
# } else {
#   print("Something wrong, double check")
# }
