##############################################################################
#
# test.R
#
# Thao P. Le
#
# 2023
# 
# Testing script
# 
##############################################################################
library(tidyverse)
library(tidyr)
source("core_functions.R")

rel_tol = 0.001
t1_prior_level = 0.95
threshold2 <- 0.005
prior_N <- 10000
prior_y <- 6
alpha = 0.5 + prior_y
beta =  0.5 + prior_N - prior_y

threshold1 <- calc_threshold_1(t1_prior_level,alpha,beta)

expected_t1<-0.00111778405

abs_error <- abs(threshold1 - expected_t1)
if (abs_error > 0.00001) {
  stop('Threshold 1 has error', abs_error)
}

expected_N<-756.9404707

calculated_N <- recommended_sample_size(threshold1, threshold2,prior_N, prior_y, rel_tol, dist_type="normal")$minimum


abs_error <- abs(expected_N - calculated_N)
if (abs_error > 0.00001) {
  stop('Recommended sample has error', abs_error)
}


################
rel_tol<- 0.00001
t1_prior_level = 0.95
threshold2 <- 0.01
prior_N <- 1000
prior_y <- 5
alpha = 0.5 + prior_y
beta =  0.5 + prior_N - prior_y

threshold1 <- calc_threshold_1(t1_prior_level,alpha,beta)
threshold1

# note that this shouldn't work because threshold 1 (T_change) is too close to threshold 2 (T_risk)
# which makes it difficult to calculate certain normalised probablities 
calculated_N <- tryCatch(
  recommended_sample_size(threshold1, threshold2,prior_N, prior_y, rel_tol, dist_type="normal")$minimum, 
  error = function(e){
    print(e)
    print("solution cannot be calculated")
  })


if (calculated_N  != "solution cannot be calculated"){
  stop("What has changed?")
}