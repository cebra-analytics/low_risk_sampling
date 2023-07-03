##############################################################################
#
# plots_probability_status.R
#
# Thao P. Le, Thomas K. Waring, Chris M. Baker
#
# 2023
# 
#  
# 
##############################################################################
library(tidyverse)
library(tidyr)
source("core_functions.R")

##############################################################################
# preparatory code
##############################################################################
rel_tol <- 0.001

data_name <- 'outputs/recommended_samples.csv'
prev_csv <- tryCatch(read.csv(data_name), error = function(e){
  stop("Didn't find past values.")
})

colour_status <- function(threshold1, threshold2, a,b){
  if(is_green(threshold1,a,b)){
    return("green")
  }
  else if(is_orange(threshold2,a,b)){
    return("orange")
  }
  else{
    return("red")
  }
  
}





get_fluid_N <- function(N_prior, y_prior,threshold1_at_prior_level,threshold2,rel_tol = 0.0001){
  prev_result <- prev_csv  %>% filter(prior_samples == N_prior, past_leakage == y_prior, abs(T1_level_percent - threshold1_at_prior_level*100)<0.00001, abs(T2_percent - threshold2*100)<0.000001)
  if (nrow(prev_result) == 0) {
    N <- tryCatch( recommended_sample_size(threshold1, threshold2,N_prior, y_prior, rel_tol, dist_type="normal")$minimum, 
                   error = function(e){
                     if (verbose) {
                       print(paste0("Exiting during reporting period ",reporting,", reason was:"))
                       print(paste0("-- ",geterrmessage()))
                     }
                     print(paste0("Calculating new result: prior_y = ", y_prior, ", N = ", N_prior, ", t1 level = ",threshold1_at_prior_level, ", t2 = ", threshold2))
                     -1
                   })
  } else {
    N <- prev_result$optimal_samples
  }
  return(as.integer(N)) 
}

##############################################################################
# TK look here !!!!!!!!!


# find y value so that posterior rate is 95% above threshold
find_y_cutoff <- function (threshold, a0, b0, N1) {
  sq_residual <- function (y1) {
    (pbeta (threshold, a0 + y1, b0 + N1 - y1) - 0.95)^2
  }
  result <- optimize (sq_residual, interval=c(0,100))
  return (result$minimum)
}

# find probabilities of colour
colour_probs <- function (t1, t2, a0, b0, N1, rate) {
  y_orange <- find_y_cutoff (t1, a0, b0, N1)
  y_red <- find_y_cutoff (t2, a0, b0, N1)

  pr_green <- pbinom(y_orange, size=N1, prob=rate)
  pr_orange <- pbinom(y_red, size=N1, prob=rate) - pr_green
  pr_red = 1 - pr_green - pr_orange

  return(c('green'=pr_green, 'orange'=pr_orange, 'red'=pr_red))
}

prior_N <- 10000
prior_y <- 6
threshold1_at_prior_level<-0.95
threshold2 <- 0.005

alpha = 0.5 + prior_y
beta =  0.5 + prior_N - prior_y
threshold1 <- calc_threshold_1(threshold1_at_prior_level,alpha,beta)

N <- get_fluid_N(prior_N, prior_y,threshold1_at_prior_level,threshold2 )
print(N)


true_rate <- seq(0.0055, 0.05,0.0005)
length(true_rate)

M_run_times = 100

true_rate_list <- rep(true_rate, each=3)
colour <- rep(c("green","orange","red"), length(true_rate))
value <- rep(0,length(true_rate_list))


colour_status_results <- data.frame(true_rate_list,colour,value)

for( rate in true_rate){

  prs <- colour_probs (threshold1, threshold2, alpha, beta, N, rate)
  colour_status_results$value[colour_status_results$true_rate_list==rate] <- prs

  # for(i in seq(1,M_run_times,1)){
  #   samples <- rbernoulli(N, rate)
    
  #   y <- sum(samples)
    
  #   print(y)
    
    
  #   a <- 0.5 + prior_y + y
  #   b <- 0.5 + prior_N + N - prior_y - y
    
  #   colour <- colour_status(threshold1, threshold2, a,b)
    
  #   print(colour)
  #   colour_status_results$value[colour_status_results$true_rate_list==rate & colour_status_results$colour==colour] <- colour_status_results$value[colour_status_results$true_rate_list==rate & colour_status_results$colour==colour]+1
    
  # }
}

colour_status_results

# gather results and plot a stacked bar graph.
colour_key = c("green"="chartreuse3", "orange"="orange", "red"="brown1")

# ggplot(colour_status_results, aes(fill=colour, y=value, x=true_rate_list)) + 
#   geom_bar(position="fill", stat="identity")+
#   scale_fill_manual(values = colour_key, name='Colour status') +
#   xlab('True rate') + ylab('Probability of colour status')

ggplot(colour_status_results, aes(fill=colour, y=value, x=true_rate_list)) + 
  geom_area()+
  scale_fill_manual(values = colour_key, name='Colour status') +
  xlab('True rate') + ylab('Probability of colour status')


ggsave(paste0("outputs/plot_probability_status_new.png"),width = 12, height = 4)


################################################################################
# starting from green instead

prior_N <- 10000
prior_y <- 6
threshold1_at_prior_level<-0.95
threshold2 <- 0.005

alpha = 0.5 + prior_y
beta =  0.5 + prior_N - prior_y
threshold1 <- calc_threshold_1(threshold1_at_prior_level,alpha,beta)

N <- get_fluid_N(prior_N, prior_y,threshold1_at_prior_level,threshold2 )
print(N)


true_rate <- seq(0.001, 0.02,0.0005)
length(true_rate)

M_run_times = 100

true_rate_list <- rep(true_rate, each=3)
colour <- rep(c("green","orange","red"), length(true_rate))
value <- rep(0,length(true_rate_list))


colour_status_results <- data.frame(true_rate_list,colour,value)

for( rate in true_rate){
  for(i in seq(1,M_run_times,1)){
    samples <- rbernoulli(N, rate)
    
    y <- sum(samples)
    
    print(y)
    
    
    a <- 0.5 + prior_y + y
    b <- 0.5 + prior_N + N - prior_y - y
    
    colour <- colour_status(threshold1, threshold2, a,b)
    
    print(colour)
    colour_status_results$value[colour_status_results$true_rate_list==rate & colour_status_results$colour==colour] <- colour_status_results$value[colour_status_results$true_rate_list==rate & colour_status_results$colour==colour]+1
    
  }
}

colour_status_results

# gather results and plot a stacked bar graph.
colour_key = c("green"="chartreuse3", "orange"="orange", "red"="brown1")

# ggplot(colour_status_results, aes(fill=colour, y=value, x=true_rate_list)) + 
#   geom_bar(position="fill", stat="identity")+
#   scale_fill_manual(values = colour_key, name='Colour status') +
#   xlab('True rate') + ylab('Probability of colour status')

ggplot(colour_status_results, aes(fill=colour, y=value, x=true_rate_list)) + 
  geom_area()+
  scale_fill_manual(values = colour_key, name='Colour status') +
  xlab('True rate') + ylab('Probability of colour status')


ggsave(paste0("outputs/plot_probability_status_100_0.05.png"),width = 12, height = 4)



################################################################################
# fixed sampling:

colour_status_fixed <- function(leakage){
  if(leakage==0){
    return("green")
  }
  else if(leakage==1){
    return("orange")
  }
  else{
    return("red")
  }
  
}

true_rate <- seq(0.0055, 0.05,0.0005)

M_run_times = 100

N <- 600 # fixed at 600

true_rate_list <- rep(true_rate, each=3)
colour <- rep(c("green","orange","red"), length(true_rate))
value <- rep(0,length(true_rate_list))

colour_status_results_fixed <- data.frame(true_rate_list,colour,value)

for( rate in true_rate){
  for(i in seq(1,M_run_times,1)){
    samples <- rbernoulli(N, rate)
    
    y <- sum(samples)
    
    print(y)
    
    colour <- colour_status_fixed(y)
    
    print(colour)
    colour_status_results_fixed$value[colour_status_results_fixed$true_rate_list==rate & colour_status_results_fixed$colour==colour] <- colour_status_results_fixed$value[colour_status_results_fixed$true_rate_list==rate & colour_status_results_fixed$colour==colour]+1
    
  }
}

colour_status_results_fixed

colour_key = c("green"="chartreuse3", "orange"="orange", "red"="brown1")

ggplot(colour_status_results_fixed, aes(fill=colour, y=value, x=true_rate_list)) + 
  geom_area()+
  scale_fill_manual(values = colour_key, name='Colour status') +
  xlab('True rate') + ylab('Probability of colour status')
ggsave(paste0("outputs/plot_probability_status_100_0.05_fixed.png"),width = 12, height = 4)

