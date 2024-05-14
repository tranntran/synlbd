# This code can be used to generate the noisy joint distribution of first year and 
# last year in the LBD at the 3-digit NAICS. The output is the noisy joint probability
# saved as a csv file, so that we can use this file to generate synthetic year data 
# at 6-digit NAICS level later. 

# The strategy to synthesize year is similar to Kinney et al (2011), where we count 
# the number of observations in each category and use the Dirichlet-Multinomial 
# distribution to sample the synthetic data. The main difference is that this code
# satisfy DP (by adding Laplace noise or Exponential noise).

# We are calculatiing at the 3-digit NAICS level because the data can get sparse 
# when we look at the joint distribution for some industries. This can also help
# save runtime because we only need to do this step once for every 3-digit NAICS.

# Last updated: 4/22/24
# added parallel code

# libraries needed
library(dplyr)

generate_noisy_probs = function(data, epsilon, noise_dist = c('laplace', 'geom'), seed) {
  # input variable 
  # data: merged_by_lbdnum_[sic2].csv
  # output: ./output/year_noisy_probs/year_noisy_probs_[sic2].csv
  sic2 = floor(data$sic3[1]/10)
  set.seed(seed)
  
  # create all possible combination between first and last year
  all_year = expand.grid(c(1977:2000), c(1977:2000))
  colnames(all_year) = c('firstyear', 'lastyear')
  all_year = all_year[! all_year$firstyear > all_year$lastyear, ]
  all_year = all_year[order(all_year$firstyear, decreasing = FALSE),]
  
  # tabulate the frequency of each combination of first and last year from data
  data = data[, c('firstyear', 'lastyear')]
  year_freq = as.data.frame(table(data))
  year_freq[c('firstyear', 'lastyear')] = sapply(year_freq[c('firstyear', 'lastyear')], 
                                                 function(x) as.numeric(as.character(x)))
  year_freq = year_freq[! year_freq$firstyear > year_freq$lastyear, ]
  
  # left join to create master table and set the NA values to be 0 
  all_year = dplyr::left_join(x = all_year, y = year_freq, by = c('firstyear', 'lastyear'))
  all_year$Freq[is.na(all_year$Freq)] = 0
  
  # add laplace or geom noise to each count
  if (noise_dist == 'laplace') {
    all_year$Freq_noisy = all_year$Freq + rmutil::rlaplace(nrow(all_year), 0, 1/epsilon)
  } else if (noise_dist == 'geom') {
    p = 1 - exp(-epsilon)
    all_year$Freq_noisy = all_year$Freq + (rgeom(nrow(all_year), p) - rgeom(nrow(all_year), p))
  } else if (noise_dist == 'none') {
    all_year$Freq_noisy = all_year$Freq
    
  }
  
  all_year$Freq_noisy[all_year$Freq_noisy < 0] = 0
  
  # calculate the probability of last year given first year
  # all_year = group_by(all_year, firstyear_emp) %>% mutate(Prob_noisy = Freq_noisy/sum(Freq_noisy))
  all_year$prob_noisy = all_year$Freq_noisy/sum(all_year$Freq_noisy)
  all_year$prob_noisy[is.na(all_year$prob_noisy)] = 0
  
  all_year = all_year[, c('firstyear', 'lastyear', 'prob_noisy')]
  
  write.csv(all_year, paste0("./output/year_noisy_probs/year_noisy_probs_", noise_dist, "_", sic2, "_seed", seed, ".csv"), row.names = FALSE)
  
}

sic2_list = c(23, 54, 82)



#naics2 = 23 # parameter to change
#naics4_data = read.csv('../../data/naics4_inds/naics4_inds.csv', header = FALSE)
#naics4_data$naics3 = floor(naics4_data$V1/10)
#naics3_list = unique(naics4_data$naics3[floor(naics4_data$V1/100) == naics2])
# 
# library(doParallel)
# library(tictoc)
# num_cores= 5 #detectCores()-1 #use all available core but 1
# print(num_cores)
# workers=makeCluster(num_cores,type="SOCK")
# registerDoParallel(workers)
# 
# tic()

# Out = foreach(seed=c(1:5), .errorhandling='stop') %dopar% {
for (seed in 1:5) {
  library(dplyr)
  for (sic2 in sic2_list) {
    print(sic2)
    data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
    data$lastyear[data$lastyear > 2000] = 2000
    data$firstyear[data$firstyear < 1977] = 1977
    generate_noisy_probs(data = data, epsilon = 1, noise_dist = 'laplace', seed = seed)
  }
}

# toc()
# stopCluster(workers)