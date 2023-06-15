##############################################################################
#
# "600 samples" scenarios
#
# 2023
# 
# Here we run scenarios with fixed sample number
# 
##############################################################################
library(tidyverse)
library(tidyr)
source("allocate.r")

##############################################################################

# it will be red if there is any detected leakage 
colour_status_fixed <- function(leakage){
  if(leakage==0){
    return("green")
  }
  else{
    return("red")
  }
  
}

# doesn't use prior data
simulate_quarters_fixed <- function(num_quarters,rates, sample_size=600,verbose=FALSE){
  
  
  N_t <- rep(c(NA),num_quarters)
  y_t <- rep(c(NA),num_quarters)
  
  prior_N_t <- rep(c(NA),num_quarters)
  prior_y_t <- rep(c(NA),num_quarters)
  colour_t <- rep(c(NA),num_quarters)
  a_post_t <- rep(c(NA),num_quarters)
  b_post_t <- rep(c(NA),num_quarters)
  
  
  for (reporting in 1:(num_quarters)){
    
    print(paste0("reporting period: ", reporting))
    
    # get the next sample size
    N <- sample_size
    
    print("N:")
    print(N)
    
    # now onto sampling!
    N_t[reporting] <- N
    
    rate <- rates[reporting]
    samples <- rbernoulli(N, rate)
    
    y <- sum(samples)
    y_t[reporting] <- y
    
    
    # there is probably no clear break condition here
    if (y>0){
      print("Detected contamination, should this be a break condition?")
    }
    
    
    # we don't use this
    # a <- 0.5 + y_prior + y
    # a_post_t[reporting] <- a
    # b <- 0.5 + N_prior + N - y_prior - y
    # b_post_t[reporting] <- b
    
    colour_t[reporting] <-colour_status_fixed(y)
    
  }
  
  print(paste0("num quarters: ",num_quarters))
  
  simulation_quarters = seq(1,num_quarters,1)
  
  quarter_nums = seq(1,num_quarters)
  
  df <- data.frame(quarter_nums = quarter_nums, N_t = as.numeric(N_t[simulation_quarters]),
                   y_t = as.numeric(y_t[simulation_quarters]),
                   rates = as.numeric(rates[quarter_nums]),
                   colour_t = colour_t[quarter_nums],
                   a_post_t = as.numeric(a_post_t[quarter_nums]),
                   b_post_t = as.numeric(b_post_t[quarter_nums]))
  
  df$quarter_nums <- factor(df$quarter_nums, levels=df$quarter_nums)
  
  df
  
}



################################################################################
# plotting
################################################################################

num_quarters<- 12
low_rate <- 0.0012
high_rate <- 0.02
sample_size<-600

pathway_name<-"Fixed_Scenario_10k_6"

take_previous_method <- list("quarters", 2)

routine_save_name<- paste0(pathway_name,"-",sample_size ,"_", low_rate, "_")
highrisk_save_name<- paste0(pathway_name,"-",sample_size,"_", high_rate, "_")

routine_rates <- rep(low_rate, num_quarters)

num_routine <- 4
highrisk_rates <-c(rep(low_rate,num_routine),rep(high_rate,num_quarters-num_routine)) 


############################### fluid sampling

sampling_method = list("fixed")
additional_save_name <- paste0("prior_",take_previous_method[[1]],"_",sampling_method[[1]])



simulation_routine <-simulate_quarters_fixed(num_quarters,routine_rates,sample_size)
plot_sim(simulation_routine,"routine",routine_save_name,additional_save_name,sampling_method[[1]])

simulation_red <- simulate_quarters_fixed(num_quarters,highrisk_rates,sample_size)
plot_sim(simulation_red,"risky",highrisk_save_name,additional_save_name,sampling_method[[1]])




























