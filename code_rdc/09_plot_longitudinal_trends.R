# code to plot job creation rate and net job creation rate
rm(list = ls())

library(ggplot2)
library(dplyr)
library(data.table)
drops = c("firstyear", "lastyear", paste0("mu_", seq(1977, 2002, 5)))
sic =  178 #c(178, 239, 542, 829)


jc = list()
njc = list()

methods = c("opmnc2", 'stepkng1', 'swkng')
for (seed in 1:5){
    print(seed)
    tab_jc = matrix(NA, nrow = 4, ncol = 23)
    tab_jd = matrix(NA, nrow = 4, ncol = 23)
    tab_n = matrix(NA, nrow = 4, ncol = 23)
    df = NULL
    for (k in 1:length(sic)){
        tmp = read.csv(paste0("./output/merged_by_lbdnum/merged_by_lbdnum_", floor(sic/10), ".csv"))
        tmp = tmp[, paste0('emp_', c(1977:2000))]
        df = rbind(df, tmp)
    }
    # df[df == ""] = NA
    # df = df[, !drops, with = FALSE]
    
    for (t in 2:24){
        keep = c((t-1):t)
        df_t = na.omit(df[, keep])
        z_et = rowSums(df_t)/2
        z = sum(z_et)
        diff = df_t[,2] - df_t[,1]
        diff = cbind(diff, rep(0, length(diff)))
        tab_jc[1, (t-1)] = sum(abs(apply(diff, 1, max))/z)
        tab_jd[1, (t-1)] = sum(abs(apply(diff, 1, min))/z)
        tab_n[1, (t-1)] = tab_jc[1, (t-1)] - tab_jd[1, (t-1)]
        
    }
    
    for (i in 1:length(methods)){
        df = NULL
        for (k in 1:length(sic)){
            tab = read.csv(paste0('./output/', sic, '/', sic, '_syndata_', methods[i], '_logistic_seed', seed, '.csv'))
            tab = tab[, paste0('emp_', c(1977:2000))]
            df = rbind(df, tab)
        }
        
        for (t in 2:24){
            keep = c((t-1):t)
            df_t = na.omit(df[, keep])
            z_et = rowSums(df_t)/2
            z = sum(z_et)
            diff = df_t[,2] - df_t[,1]
            diff = cbind(diff, rep(0, length(diff)))
            tab_jc[(i+1), (t-1)] = sum(abs(apply(diff, 1, max))/z)
            tab_jd[(i+1), (t-1)] = sum(abs(apply(diff, 1, min))/z)
            tab_n[(i+1), (t-1)] = tab_jc[(i+1), (t-1)] - tab_jd[(i+1), (t-1)]
            
        }
        
    }
    
    jc[[seed]] = tab_jc*100
    njc[[seed]] = tab_n*100
}

######################## 
# plot job creation rate by year

plot_jc = matrix(NA, ncol = 23, nrow = 4)
plot_jc_sd = matrix(NA, ncol = 23, nrow = 4)
for (i in 1:4){
    tmp = rbind(jc[[1]][i,], jc[[2]][i,], jc[[3]][i,], jc[[4]][i,], jc[[5]][i,])
    plot_jc[i,] = apply(tmp, 2, mean)
    plot_jc_sd[i,] = apply(tmp, 2, sd)
}

vars = paste0("emp_", c(1978:2000))
plot_jc = as.data.table(plot_jc)
plot_jc_sd = as.data.table(plot_jc_sd)
colnames(plot_jc) = vars
colnames(plot_jc_sd) = vars
plot_jc$methods= c("Raw Data", "OPM", 
                   "Stepwise Varying", "Sandwich Varying") 
plot_jc_sd$methods = c("Raw Data", "OPM", 
                       "Stepwise Varying", "Sandwich Varying")  

plot_jc= melt(plot_jc, id.vars = c("methods"), measure.vars = vars)
plot_jc_sd = melt(plot_jc_sd, id.vars = c("methods"), measure.vars = vars)
plot_final = merge(plot_jc, plot_jc_sd, by = c('methods', 'variable'), all.x = T)
colnames(plot_final) = c('methods', 'variable', 'mean', 'sd')
plot_final$min_val = plot_final$mean - plot_final$sd
plot_final$max_val = plot_final$mean + plot_final$sd
plot_final$year = as.numeric(substr(plot_final$variable, 5,8))


# Dynamically generate default color values, but have Raw Data = "black"
adj_names = sort(setdiff(unique(plot_final$methods), "Raw Data"))
thickness = c(1.5, 1.5, 1.5, 1.5)
names(thickness) = c(adj_names, 'Raw Data')

palette = c("#000000", "#F11B00", "#3B9AB2", "#E1AF00")
names(palette) = c('Raw Data', adj_names)

png(paste0('./output/plots/', sic3, '_njcr.png'), width = 750, height = 500)
ggplot(data = plot_final, aes(x = year, group = methods)) + 
    geom_line(aes(y = mean, color = methods, size = methods)) + 
    geom_ribbon(aes(y = mean, ymin = min_val, ymax = max_val, fill = methods), alpha = .3) +
    labs(title="Job Creation Rate by Year", 
         x ="Year", y = "Job Creation Rate (Percent)", 
         color = 'Method', fill = 'Method', size = 'Method') +
    scale_x_continuous(breaks=seq(1980, 2000, by = 5)) + 
    scale_colour_manual(values=palette) +
    scale_fill_manual(values=palette) +
    scale_size_manual(values = thickness) +
    theme(axis.text = element_text(size=13),
          legend.position = 'bottom',
          plot.title = element_text(size = 14, face = "bold"),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=13),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'),
          axis.title=element_text(size=13))
dev.off()
######################## 
# plot net job creation rate by year

plot_njc = matrix(NA, ncol = 23, nrow = 4)
plot_njc_sd = matrix(NA, ncol = 23, nrow = 4)
for (i in 1:4){
    tmp = rbind(njc[[1]][i,], njc[[2]][i,], njc[[3]][i,], njc[[4]][i,], njc[[5]][i,])
    plot_njc[i,] = apply(tmp, 2, mean)
    plot_njc_sd[i,] = apply(tmp, 2, sd)
}

vars = paste0("emp_", c(1978:2000))
plot_njc = as.data.table(plot_njc)
plot_njc_sd = as.data.table(plot_njc_sd)
colnames(plot_njc) = vars
colnames(plot_njc_sd) = vars
plot_njc$methods= c("Raw Data", "OPM", 
                    "Stepwise Varying", "Sandwich Varying") 
plot_njc_sd$methods = c("Raw Data", "OPM", 
                        "Stepwise Varying", "Sandwich Varying") 

plot_njc= melt(plot_njc, id.vars = c("methods"), measure.vars = vars)
plot_njc_sd = melt(plot_njc_sd, id.vars = c("methods"), measure.vars = vars)
plot_final = merge(plot_njc, plot_njc_sd, by = c('methods', 'variable'), all.x = T)
colnames(plot_final) = c('methods', 'variable', 'mean', 'sd')
plot_final$min_val = plot_final$mean - plot_final$sd
plot_final$max_val = plot_final$mean + plot_final$sd
plot_final$year = as.numeric(substr(plot_final$variable, 5,8))


palette = c("#000000", "#F11B00", "#3B9AB2", "#E1AF00")
names(palette) = c('Raw Data', adj_names)

png(paste0('./output/plots/', sic3, '_njcr.png'), width = 750, height = 500)
ggplot(data = plot_final, aes(x = year, group = methods)) + 
    geom_line(aes(y = mean, color = methods, size = methods)) + 
    geom_ribbon(aes(y = mean, ymin = min_val, ymax = max_val, fill = methods), alpha = .3) +
    labs(title="Net Job Creation Rate by Year", 
         x ="Year", y = "Net Job Creation Rate (Percent)", 
         color = 'Method', fill = 'Method', size = 'Method') +
    scale_x_continuous(breaks=seq(1980, 2000, by = 5)) + 
    scale_colour_manual(values=palette) +
    scale_fill_manual(values=palette) +
    scale_size_manual(values = thickness) +
    theme(axis.text = element_text(size=13),
          legend.position = 'bottom',
          plot.title = element_text(size = 14, face = "bold"),
          legend.title=element_text(size=14), 
          legend.text=element_text(size=13),
          legend.key.size = unit(1, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'),
          axis.title=element_text(size=13))
dev.off()