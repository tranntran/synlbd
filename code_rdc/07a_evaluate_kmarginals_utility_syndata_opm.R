#setwd("/projects/programs/lbd_naics")
library(ggplot2)
library(data.table)
library(dplyr)


sic3 = 178 #as.numeric(commandArgs(trailingOnly = TRUE)) #2373
set.seed(sic3)

k = 3
combn = 100
all_samp = t(combn(24, 3))
all_samp = all_samp[sample(1:nrow(all_samp), combn), ]
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
            dt = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_swkng_logistic_seed", j, ".csv"))
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
    png(filename = paste0('./output/kmarginals/', sic3, '_kmarginals_swkng.png'), width = 650, height = 450)
    hist(apply(nist_score, 1, mean), breaks = 30, 
         main = paste0('Sandwich KNG'),
         xlab = 'Score')
    dev.off()
    #hist(apply(nist_score, 1, sd)/sqrt(5))
    
} else {
    print("Something wrong, double check")
}
