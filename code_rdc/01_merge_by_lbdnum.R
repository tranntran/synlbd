# # This code is used to merge establishments longitudinally using lbdnum as the key.
# # After merging, it will output files merged_by_lbdnum_[naics3].csv, where [naics3] 
# # is all the unique naics3 within the naics2 parameter, which needs to be specified 
# # bewlow by users. Output file merged_by_lbdnum_[naics3].csv will have the following 
# # columns lbdnum, firstyear_emp, lastyear_emp, bds_vcnaics, bds_vcnaics3, st, mu, 
# # emp_[year] (year range from 1977 to 2021).
# 
# # By doing this, we can help save time from reading in data and merging establishments
# # when modeling at the NAICS6 level. However, we work with each unique NAICS3 within 
# # the specified NAICS2 (instead of just saving every possible NAICS3), because saving 
# # everything will drain the memory space. Additionally, BDS tabulation code is at NAICS2
# # level. To automate this code, we can parse the naics2 parameter from the command line.
# 
# # Last updated: 04/22/24
# Sys.time()
# # parameter to change
#naics2 = 23 # as.numeric(commandArgs(trailingOnly = TRUE))
# 
## Libraries needed
library(dplyr)

## Now, we want to read in the relevant data and merge the data frames together
vars.list=c("lbdnum", "firstyear", "lastyear", "sic3", "mu", "emp", "pay")

## Now, read in for every year and merge together
for(i in 1977:2000){ #2021){
    read.file = paste("../../data/synlbd", i, ".csv", sep="")
    df.tmp = read.csv(read.file)
    df.tmp = df.tmp[, vars.list] ##Restrict on variables
    names(df.tmp)[2:length(vars.list)] = paste(vars.list, "_", i, sep="")[2:length(vars.list)]##
    if(i == 1977){df.main = df.tmp} else {df.main = dplyr::full_join(df.main, df.tmp, by="lbdnum")}
    print(i)
}#end of for
print('complete merge')
Sys.time()
# 
# ## Now, we only need one first year var, last year var, st var, and mu var
# names(df.main)
# ## First, we want to reorganize the variables
# ## We will do this with a loop
fyemp = c()
lyemp = c()
sic3 = c()
mu = c()
emp = c()
pay = c()
for(i in 1977:2000){
    fyemp = c(fyemp, paste("firstyear_", i, sep = ""))
    lyemp = c(lyemp, paste("lastyear_", i, sep = ""))
    sic3 = c(sic3, paste("sic3_", i, sep = ""))
    mu = c(mu, paste("mu_", i, sep = ""))
    emp = c(emp, paste("emp_", i,sep = ""))
    pay = c(pay, paste("pay_", i,sep = ""))
    
}

names = c("lbdnum", fyemp, lyemp, sic3, mu, emp, pay)

df.main = df.main[, names] ##Reordering all the variables

##Now, we need to check that the data is behaving like we think it should be
##First, just checking the first few rows
head(df.main)

##Now, we need to only construct one firstyear, lastyear, and mu variables
##Checked for validity. They are all the same across years
nyears = 2000 - 1977 + 1 #2021

#df.main = head(df.main, 100000)
# tmp = head(df.main, 100)
# test = as.vector(unlist(apply(tmp[, (2+2*nyears):(2+3*nyears-1)], 1, FUN = function(x) names(which.max(table(x))))))

df.main$firstyear=apply(df.main[, 2:(2+nyears-1)],1,FUN=min,na.rm=TRUE)
df.main$lastyear=apply(df.main[, (2+nyears):(2+2*nyears-1)], 1,FUN=min,na.rm=TRUE)

Sys.time()
tmp = apply(df.main[, (2+2*nyears):(2+3*nyears-1)], 1, FUN = function(x) names(which.max(table(x))))
tmp[sapply(tmp, is.null)] = NA
df.main$sic3= unlist(tmp)
df.main$sic2 = as.numeric(substr(df.main$sic3, 1, 2))
Sys.time()
# subset to naics2 to reduce running time and memory usage
# df.main = df.main[floor(df.main$bds_vcnaics3/10) == naics2, ]


# tmp = apply(df.main[, (2+3*nyears):(2+4*nyears-1)], 1, FUN = function(x) names(which.max(table(x))))
# tmp[sapply(tmp, is.null)] = NA
# df.main$st= as.numeric(unlist(tmp))


#write.csv(df.main, paste0("./output/merged_by_lbdnum/merged_by_lbdnum_", naics2, ".csv"), row.names = FALSE)
# nyears = 2000 - 1977 + 1 #2021
# df.main = read.csv("./output/merged_by_lbdnum/merged_by_lbdnum_23.csv")

df_tmp = NULL
tmp = df.main[, (2+3*nyears):(2+4*nyears-1)]
#tmp = head(tmp, 1000)

for (i in 0:5) {
    print(i)
    if (i == 0) {
        df_tmp = cbind(df_tmp, tmp[, 1])
    } else {
        tmp_tmp = apply(tmp[, c((5*(i-1)+2): min(24, 5*i+1))], 1, FUN = function(x) names(which.max(table(x))))
        tmp_tmp[sapply(tmp_tmp, is.null)] = NA
        df_tmp= cbind(df_tmp, as.numeric(unlist(tmp_tmp)))
    }
    
}
df_tmp = as.data.frame(df_tmp)
mu5 = paste0('mu5_', seq(1977, 2002, 5))
colnames(df_tmp) = mu5
df.main = cbind(df.main, df_tmp)

# test = matrix(c(1, 0, 0, 0, NA), nrow = 1)

# tmp = apply(, 1, FUN = function(x) names(which.max(table(x))))
# tmp[sapply(tmp, is.null)] = NA
# df.main$mu= as.numeric(unlist(tmp))

##Now, subsetting the variables into the smaller data frame
names=c("lbdnum", "sic3", "firstyear","lastyear", mu5, emp, pay)
df.main=df.main[, names]
write.csv(df.main, paste0("./output/merged_by_lbdnum/merged_by_lbdnum.csv"), row.names = FALSE)

df.main = read.csv("./output/merged_by_lbdnum/merged_by_lbdnum.csv")
sic2_list = sort(unique(floor(df.main$sic3/10)))
sic2_list = sic2_list[sic2_list>9]
# naics4_data = read.csv('../../data/naics4_inds/naics4_inds.csv', header = FALSE)
# naics4_data$naics3 = floor(naics4_data$V1/10)
# naics3_list = unique(naics4_data$naics3[floor(naics4_data$V1/100) == naics2])

for (sic2 in sic2_list) {
    print(sic2)
    tmp = df.main[!is.na(df.main$sic3) & floor(df.main$sic3/10) == sic2, ]
    write.csv(tmp, paste0("./output/merged_by_lbdnum/merged_by_lbdnum_", sic2, ".csv"), row.names = FALSE)
    
}
Sys.time()



## debugging 
# test = read.csv("./output/merged_by_lbdnum/merged_by_lbdnum_237.csv")
# test = head(test, 100)
# test_lbdnum =test$lbdnum
# 
# df.main = read.csv("./output/merged_by_lbdnum/merged_by_lbdnum_23.csv")
# df.main = df.main[df.main$lbdnum %in% test_lbdnum, ]
# tmp = df.main[, (2+4*nyears):(2+5*nyears-1)]
# #tmp = head(tmp, 1000)
# df_tmp = NULL
# for (i in 0:9) {
#     print(i)
#     if (i == 0) {
#         df_tmp = cbind(df_tmp, tmp[, 1])
#     } else {
#         tmp_tmp = apply(tmp[, c((5*(i-1)+2): min(45, 5*i+1))], 1, FUN = function(x) names(which.max(table(x))))
#         tmp_tmp[sapply(tmp_tmp, is.null)] = NA
#         df_tmp= cbind(df_tmp, as.numeric(unlist(tmp_tmp)))
#     }
#     
# }
# df_tmp = as.data.frame(df_tmp)
# mu5 = paste0('mu_', seq(1977, 2022, 5))
# colnames(df_tmp) = mu5
# df.main = cbind(df.main, df_tmp)