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
all_data$type = factor(all_data$type, c('Raw Data', 'OPM', 'Stepwise KNG', 'Sandwich KNG'))

i = 32
all_data_tmp = all_data[, c(i, 33)]
colnames(all_data_tmp) = c('value', 'type')
#all_data_tmp = all_data_tmp[all_data_tmp$type %in% c('opmv2b', 'opmv2d', 'raw data'), ]
png(paste0('./output/plots/', sic3, "_", colnames(syn_data)[i], '.png'), width = 500, height = 350)
print(ggplot(data = all_data_tmp, aes(x = value, fill = type)) + 
          geom_histogram(position = 'identity', alpha = 1, binwidth = 1) + #, color = wes_palette('Zissou1')) +
          labs(subtitle = 'Distribution of Establishment Employment for SIC 178 in 2000', x = 'Employee Numbers', y = 'Count') +
          #labs(subtitle = colnames(syn_data)[i]) + 
          scale_x_continuous(limits = c(0, 40)) + 
          facet_wrap(~type) + 
          theme_minimal() + theme(legend.position = 'none')+
          scale_fill_manual(values = c('#999999', "#F11B00", "#78B7C5", "#EBCC2A")) #+ 
)
dev.off()
print(summary(all_data_tmp$value[all_data_tmp$type  == 'OPM']))
print(summary(all_data_tmp$value[all_data_tmp$type  == 'Stepwise KNG']))
print(summary(all_data_tmp$value[all_data_tmp$type  == 'Sandwich KNG']))
print(summary(all_data_tmp$value[all_data_tmp$type  == 'Raw Data']))
print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'OPM']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)))
print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Stepwise KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)))
print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Sandwich KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)))
print(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Raw Data']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98,0.99, 0.995)))


tab = rbind(quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'OPM']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)),
            quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Stepwise KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)),
            quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Sandwich KNG']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98, 0.99, 0.995)),
            quantile(na.omit(all_data_tmp$value[all_data_tmp$type  == 'Raw Data']), c(seq(0.1, 0.9, 0.1), 0.93, 0.95, 0.97, 0.98,0.99, 0.995)))
rownames(tab) = c('OPM', 'Stepwise KNG', 'Sandwich KNG', 'Raw Data')
tab = as.data.frame(tab)
write.csv(tab, './disclosure/178_emp2000_summary_statistics.csv')

