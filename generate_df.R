##############################################################################
#
# generate_df.R
#
# Thao P. Le, Thomas K. Waring, Chris M. Baker
#
# 2023
# 
# This file provides the code to calculate the recommended minimum sampling
# for different prior data and leakage/risk tolerances.
#
# After running, the code saves the results into "recommended_samples.csv"
# 
##############################################################################
library(tidyverse)
library(tidyr)
source("core_functions.R")


# run and save

prior_N_list <- seq(500,10E3,500)
t1_prior_level_list <- c(.94,.95,.96)
t2_list <- seq(0.004,0.01,0.001)
rel_tol = 0.001

largest_measured_rate <- 3E-3
max_num_y <- floor(prior_N_list[length(prior_N_list)]*largest_measured_rate)+1
n_rows <- max_num_y*length(t1_prior_level_list)*length(t2_list)

DF <- data.frame(prior_samples=rep(NA, n_rows),
                 past_leakage=rep(NA, n_rows),
                 T1_level_percent=rep(NA, n_rows),
                 T1_percent = rep(NA,n_rows),
                 T2_percent=rep(NA, n_rows),
                 optimal_samples=rep(NA, n_rows),
                 stringsAsFactors=FALSE)

data_name <- 'outputs/recommended_samples.csv'
prev_csv <- tryCatch(read.csv(data_name), error = function(e) {
  DF
})

overwrite = FALSE # set to true to ignore cached values
i = 1
for (prior_N in prior_N_list) {
  prior_y_list <- seq(0,floor(prior_N*largest_measured_rate),1)
  for (t1_prior_level in t1_prior_level_list) {
    for (t2 in t2_list) {
      for (prior_y in prior_y_list) {
        alpha = 0.5 + prior_y
        beta =  0.5 + prior_N - prior_y
        t1 <- calc_threshold_1(t1_prior_level,alpha,beta)
        if (t1<=t2 && !is_red(t2,alpha,beta) ){
          prev_result <- prev_csv  %>% filter(prior_samples == prior_N, past_leakage == prior_y, abs(T1_level_percent - t1_prior_level*100)<0.00001, abs(T2_percent- t2*100)< 0.0000001)
          found_result = TRUE
          if (nrow(prev_result) == 0 || overwrite) {
            # print(paste0("prior samples = ", prior_N, ", past leakage = ",prior_y,", T1_level_percent = ",t1_prior_level*100,", T2_percent = ",t2*100))
            print(paste0("Calculating new result: prior_y = ", prior_y, ", N = ", prior_N, ", t1 level = ",t1_prior_level, ", t2 = ", t2))
            opt_samples <- tryCatch(
              recommended_sample_size (t1, t2, prior_N, prior_y, rel_tol, dist_type="normal")$minimum, 
              error = function(e){
                found_result = FALSE
              })
          } else {
            t1 <- prev_result$T1_percent / 100
            opt_samples <- prev_result$optimal_samples
          }
          
          if (found_result && opt_samples>0) {
            DF[i,] <- list(prior_N,prior_y,t1_prior_level*100,t1*100,t2*100,opt_samples)
            i <- i+1
          }
          
        }
        
      }
    }
  }
}

# print(DF)
write.csv(DF,data_name, row.names = FALSE)