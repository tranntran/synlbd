# This code can be used to evaluate the synthetic data with first and last year, as well
# as multiunit status.

# Last updated: 3/13/24

# load packages
library(dplyr)
library(synthpop)

# functions to calculate pMSE scsore (Snoke et al, 2018)

# utility_glmnet = function(data, alpha){
#     x = model.matrix(inds ~ ., data = data)[,-1]
#     y = as.factor(inds)
#     lambda = exp(seq(log(0.001), log(1), length.out = 100))
#     cv_glmnet = glmnet::cv.glmnet(x, y, alpha = alpha, family = "binomial", nfolds=10, lambda= lambda,
#                                   type.measure = "auc")
#     preds = predict(cv_glmnet, newx = x, type = "response", 
#                     s= cv_glmnet$lambda.1se, exact=FALSE)
#     score = sum((preds-0.5)^2)/nrow(data)
#     return(score)
#     
# }

# naics2 = 23 # parameter to change
# sic3_data = read.csv('../../data/sic3_inds/sic3_inds.csv', header = FALSE)
# colnames(sic3_data) = 'sic3'
# sic3_list = unique(sic3_data$sic3[floor(sic3_data$sic3/100) == naics2])
# sic3_list = sic3_list[-1]

sic3_list = c(178, 239, 542, 829)
# tab = read.csv(paste0("./output/general_utility/utility_syndata_yearmu.csv"))
var_list = c('firstyear', 'lastyear', paste('mu5', seq(1977, 2002, 5), sep = '_'))


for (i in 2:length(sic3_list)) {
    sic3 = sic3_list[i]
    print(sic3)
    for (seed in 1:5) {
        print(seed)
        tab = data.frame(sic3 = rep(NA, 1),
                         pMSE_logit = rep(NA, 1),
                         pMSE_logit_inter = rep(NA, 1),
                         pMSE_CART = rep(NA, 1))
        tab$sic3[1] = sic3
        sic2 = floor(sic3/10)
        syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_seed", seed, ".csv"))
        syn_data = syn_data[,-which(colnames(syn_data) %in% c('lbdnum', 'sic3'))]
        
        raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
        raw_data = raw_data[raw_data$sic3 == sic3, ]
        raw_data = raw_data[, var_list]
        colnames(raw_data) = colnames(syn_data)
        
        
        tab$pMSE_logit[1] = utility.gen(syn_data, raw_data,  method = 'logit', maxorder = 0)$pMSE
        tab$pMSE_CART[1] = utility.gen(syn_data, raw_data, method = "cart", nperms = 20)$pMSE
        
        output = "./output/general_utility/utility_syndata_yearmu.csv"
        if (file.exists(output)) {
            tab_og = read.csv(output)
            tab = rbind(tab_og, tab)
        }
        write.csv(tab, output, row.names = FALSE)
        
    }

}


# dat = read.csv(paste0("./output/general_utility/utility_syndata_yearmu.csv"))
# dat$naics3 = floor(dat$sic3/10)
# dat %>% group_by(naics3) %>% summarise(mean_pmse_logit = mean(pMSE_logit),
#                                        mean_pmse_logit_int = mean(pMSE_logit_inter),
#                                        mean_pmse_cart = mean(pMSE_CART))
# 
# dat %>% group_by(naics3) %>% summarise(se_pmse_logit = sd(pMSE_logit)/sqrt(n()),
#                                        se_pmse_logit_int = sd(pMSE_logit_inter)/sqrt(n()),
#                                        se_pmse_cart = sd(pMSE_CART)/sqrt(n()))
# 
# dat %>% group_by(naics3) %>% summarise(se_pmse_logit = sd(pMSE_logit),
#                                        se_pmse_logit_int = sd(pMSE_logit_inter),
#                                        se_pmse_cart = sd(pMSE_CART))
# 
# dat %>% group_by(naics3) %>% summarise(se_pmse_logit = n(),
#                                        se_pmse_logit_int = n(),
#                                        se_pmse_cart = n())
