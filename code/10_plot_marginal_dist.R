# This code can be used to generate overlaying distribution plots
# between synthetic and raw data. It is used as a way to visually assess
# the utility of the generated synthetic data.

# Last updated: 3/13/24

# load libraries
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)

sic3 = 829
sic2 = floor(sic3/10)


### year and mu variables 

var_list = c('firstyear', 'lastyear', paste('mu5', seq(1977, 2002, 5), sep = '_'))

syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syn_data_year_mu_seed5.csv"))
syn_data = syn_data[,-c(1:2)]
raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', sic2, '.csv'))
raw_data = raw_data[raw_data$sic3 == sic3, ]
raw_data = raw_data[, var_list]
raw_data$lastyear[raw_data$lastyear > 2000] = 2000
raw_data$firstyear[raw_data$firstyear < 1977] = 1977
#raw_data= raw_data[, 1:2]

colnames(raw_data) = colnames(syn_data)
n = nrow(raw_data)

all_data = rbind(raw_data, syn_data)
all_data$Type = c(rep('Raw', n), rep('Synthetic', n))

# overall
p1 = ggplot(data = all_data, aes(x = firstyear, fill = Type)) +
  geom_histogram(position = 'identity', alpha = 0.5, binwidth = 1) + 
  labs(x = 'First Year', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A")) #+ 

p2 = ggplot(data = all_data, aes(x = lastyear, fill = Type)) +
  geom_histogram(position = 'identity', alpha = 0.5, binwidth = 1) + 
  labs(x = 'Last Year', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A")) #+ 

png(paste0('./output/plots/', sic3, '_year.png'), width = 625, height = 325)
ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
dev.off()

p1 = ggplot(data = all_data, aes(x = as.factor(mu_1977), fill = Type)) +
  geom_bar(position = 'identity', alpha = 0.5) + 
  labs(subtitle = '1977', x = 'Multiunit Status', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A"))
p2 = ggplot(data = all_data, aes(x = as.factor(mu_1982), fill = Type)) +
  geom_bar(position = 'identity', alpha = 0.5) + 
  labs(subtitle = '1982', x = 'Multiunit Status', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A"))
p3 = ggplot(data = all_data, aes(x = as.factor(mu_1987), fill = Type)) +
  geom_bar(position = 'identity', alpha = 0.5) + 
  labs(subtitle = '1987', x = 'Multiunit Status', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A"))
p4 = ggplot(data = all_data, aes(x = as.factor(mu_1992), fill = Type)) +
  geom_bar(position = 'identity', alpha = 0.5) + 
  labs(subtitle = '1992', x = 'Multiunit Status', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A"))
p5 = ggplot(data = all_data, aes(x = as.factor(mu_1997), fill = Type)) +
  geom_bar(position = 'identity', alpha = 0.5) + 
  labs(subtitle = '1997', x = 'Multiunit Status', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A"))
p6 = ggplot(data = all_data, aes(x = as.factor(mu_2002), fill = Type)) +
  geom_bar(position = 'identity', alpha = 0.5) + 
  labs(subtitle = '2002', x = 'Multiunit Status', y = 'Count') +
  theme_minimal() + 
  scale_fill_manual(values = c("#78B7C5", "#EBCC2A"))
png(paste0('./output/plots/', sic3, '_mu.png'), width = 725, height = 425)
ggarrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
dev.off()


# by state
ggplot(data = all_data, aes(x = firstyear_emp, fill = type)) +
  geom_histogram(position = 'identity', alpha = 0.4, binwidth = 1) +
  facet_wrap(~ st)


ggplot(data = all_data, aes(x = lastyear_emp, fill = type)) +
  geom_histogram(position = 'identity', alpha = 0.3, binwidth = 1) +
  facet_wrap(~ st)


ggplot(data = all_data, aes(x = mu, fill = type)) +
  geom_bar(position = 'identity', alpha = 0.3) +
  facet_wrap(~ st)


### employment variables between 1 synthetic data version and the confidential data
library(ggplot2)
sic3 = 178
var_list = c('firstyear', 'lastyear', paste('mu5', seq(1977, 2002, 5), sep = '_'), 
             paste('emp', c(1977:2000), sep = '_'))

syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_swkng_logistic_seed5.csv"))
syn_data = syn_data[,-which(colnames(syn_data) %in% c('lbdnum', 'sic3'))]


raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', floor(sic3/10), '.csv'))
raw_data = raw_data[raw_data$sic3 == sic3, ]
raw_data = raw_data[, var_list]
colnames(raw_data) = colnames(syn_data)

# plot histogram of emp
all_data = rbind(raw_data, syn_data)
all_data$type = c(rep('raw data', nrow(raw_data)), rep('syn data', nrow(syn_data)))

for (i in 14:32) {
  i = 8
  i = i+1
  all_data_tmp = all_data[, c(i, 33)]
  colnames(all_data_tmp) = c('value', 'type')
  print(ggplot(data = all_data_tmp, aes(x = value, fill = type)) + 
          geom_histogram(position = 'identity', alpha = 0.4, binwidth = 1) +
          labs(subtitle = colnames(syn_data)[i]) + 
          scale_x_continuous(limits = c(0, 100))) 
  print(summary(all_data_tmp$value[all_data_tmp$type  == 'syn data']))
  print(summary(all_data_tmp$value[all_data_tmp$type  == 'raw data']))
  
  print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'syn data']), c(seq(0.1, 0.9, 0.1), seq(0.91, 0.99, 0.01), 0.995)))
  print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'raw data']), c(seq(0.1, 0.9, 0.1), seq(0.91, 0.99, 0.01), 0.995)))
}


# summary statistics visual check
for (i in 14: ncol(syn_data)) {
  print(colnames(syn_data)[i])
  print('Syn data')
  print(summary(syn_data[, i]))
  print('Raw data')
  print(summary(raw_data[, i]))
}



### employment variables between multiple synthetic versions and confidential data
library(ggplot2)
library(wesanderson)
sic3 = 178

var_list = c('firstyear', 'lastyear', paste('mu5', seq(1977, 2002, 5), sep = '_'), 
             paste('emp', c(1977:2000), sep = '_'))

syn_data = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_opmnc2_logistic_seed5.csv"))
syn_data = syn_data[,-which(colnames(syn_data) %in% c('lbdnum', 'sic3'))]
syn_data1 = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_stepkng1_logistic_seed5.csv"))
syn_data1 = syn_data1[,-which(colnames(syn_data1) %in% c('lbdnum', 'sic3'))]
syn_data2 = read.csv(paste0("./output/", sic3, "/", sic3, "_syndata_swkng_logistic_seed5.csv"))
syn_data2 = syn_data2[,-which(colnames(syn_data2) %in% c('lbdnum', 'sic3'))]

raw_data = read.csv(paste0('./output/merged_by_lbdnum/merged_by_lbdnum_', floor(sic3/10), '.csv'))
raw_data = raw_data[raw_data$sic3 == sic3, ]
raw_data = raw_data[, var_list]
colnames(raw_data) = colnames(syn_data)

# plot histogram of emp
all_data = rbind(raw_data, syn_data, syn_data1, syn_data2)
all_data$type = c(rep('Raw Data', nrow(raw_data)), 
                  rep('OPM', nrow(syn_data)),
                  rep('Stepwise KNG', nrow(syn_data1)),
                  rep('Sandwich KNG', nrow(syn_data2)))
all_data$type = factor(all_data$type, c('Raw Data', 'OPM', 'Stepwise KNG', 'Stepwise KNG'))

for (i in 14:33) {
  i = i + 1
  all_data_tmp = all_data[, c(i, 33)]
  colnames(all_data_tmp) = c('value', 'type')
  #all_data_tmp = all_data_tmp[all_data_tmp$type %in% c('opmv2b', 'opmv2d', 'raw data'), ]
  print(ggplot(data = all_data_tmp, aes(x = value, fill = type)) + 
          geom_histogram(position = 'identity', alpha = 1, binwidth = 1) + #, color = wes_palette('Zissou1')) +
          labs(subtitle = 'Distribution of Establishment Employment for SIC 178 in 2000', x = 'Employee Numbers', y = 'Count') +
          #labs(subtitle = colnames(syn_data)[i]) + 
          scale_x_continuous(limits = c(0, 40)) + 
          facet_wrap(~type) + 
          theme_minimal() + theme(legend.position = 'none')+
          scale_fill_manual(values = c('#999999', "#F11B00", "#78B7C5", "#EBCC2A")) #+ 
  )
  print(summary(all_data_tmp$value[all_data_tmp$type  == 'OPM']))
  print(summary(all_data_tmp$value[all_data_tmp$type  == 'Stepwise KNG']))
  print(summary(all_data_tmp$value[all_data_tmp$type  == 'Sandwich KNG']))
  print(summary(all_data_tmp$value[all_data_tmp$type  == 'Raw Data']))
  print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'OPM']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)))
  print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Stepwise KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)))
  print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Sandwich KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)))
  print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Raw Data']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98,0.99, 0.995)))
}

tab = rbind(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'OPM']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)),
            quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Stepwise KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)),
            quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Sandwich KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)),
            quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Raw Data']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98,0.99, 0.995)))
rownames(tab) = c('OPM', 'Stepwise KNG', 'Sandwich KNG', 'Raw Data')
tab = as.data.frame(tab)
write.csv(tab, './disclosure/178_emp2000_summary_statistics.csv')

# summary statistics visual check
for (i in 14: ncol(syn_data)) {
  print(colnames(syn_data)[i])
  print('Syn data')
  print(summary(syn_data[, i]))
  print('Raw data')
  print(summary(raw_data[, i]))
}