##############################################################################
#
# core_functions.R 
#
# Thao P. Le, Thomas K. Waring, Chris M. Baker
#
# 06/06/2023
#
# This file provides various functions to calculate the recommended sampling
# size, given various parameters and prior data sizes.
#
# The key function is :
# recommended_sample_size(threshold1, threshold2, prior_N, prior_y, rel_tol,dist_type)
# 
##############################################################################

library(tidyverse)
library(memoise)
library( Rmpfr)

# calculates the threshold 1 at a certain level relative to the beta prior B(a,b)
calc_threshold_1 <- function(threshold1_at_prior_level,a,b){
  qbeta(threshold1_at_prior_level,a,b) 
}



# checks whether a beta distribution B(a,b) satisfies the colours green/orange/red
is_green <- function(threshold1, a,b){
  pbeta(threshold1, a,b) >= 0.95
}
is_orange <-function(threshold2,a,b){
  pbeta(threshold2, a,b) >= 0.95
}
is_red <- function(threshold2,a,b){
  pbeta(threshold2, a,b) < 0.95
}


# Calculates the density function of the rate (using the beta
# distribution given by prior alpha and beta) conditional on r >= T2 
# i.e., gives the normalised probabilities at rates above threshold 2
norm_probabilites_above_T2_internal <- function(rate,threshold2,prior_N, prior_y,rel_tol=0.00001){
  if(rate<threshold2){
    stop("rate is lower than threshold 2")
  }
  
  beta_prior <- 0.5 + prior_N - prior_y
  alpha_prior <- 0.5 + prior_y
  
  bottom_prob <- 1- pbeta(threshold2, alpha_prior, beta_prior)
  if(bottom_prob>0){
    probability <- dbeta(rate, alpha_prior,beta_prior) / bottom_prob
  }
  else{
    func_to_int <- function(t){
      exp(-alpha_prior*t) * (1-exp(-t))^(beta_prior-1)
    }
    top_probability_without_Bab <- (alpha_prior-1)*log( rate) + (beta_prior-1)*log((1-rate))
    
    bottom_probability_without_Bab <- log(integrate(func_to_int, lower =0, upper =-log(threshold2),rel.tol = rel_tol)$value)
    
    probability <- exp(top_probability_without_Bab - bottom_probability_without_Bab)
    
    if (probability==Inf){
      func_to_int_mpfr <- function(t){
        exp(-alpha_prior*t) * (1-exp(-mpfr(as.character(t),128)))^(beta_prior-1)
      }
      
      top_probability_without_Bab <- (alpha_prior-1)*log( rate) + (beta_prior-1)*log((1-rate))
      bottom_probability_without_Bab <- log(integrateR(func_to_int_mpfr, lower =0, upper =-log(threshold2),rel.tol = rel_tol)$value)
      
      probability_mpfr <- exp(top_probability_without_Bab - bottom_probability_without_Bab)
      
      return(asNumeric(probability_mpfr))
    }
  }
  
  return(probability)
}

# Calculates the probability of correctly coloured posterior betas. 
# Note that rate > threshold2
prob_correct_colour <- function(new_N, rate, threshold1, threshold2,prior_N,prior_y,dist_type="normal"){
  alpha_prior <- 0.5 + prior_y
  beta_prior <- 0.5 + prior_N - prior_y
  
  binomial_version <-function(){
    new_N_int<- floor(new_N)
    
    total_prob_of_wanted_colour = 0
    
    for (new_y in seq(0, new_N_int, 1)) {
      alpha_post <- alpha_prior + new_y
      beta_post <- beta_prior + new_N_int - new_y
      
      if (rate<threshold1){
        if (!is_red(threshold2,alpha_post, beta_post)){ # want green or orange
          total_prob_of_wanted_colour <- total_prob_of_wanted_colour +dbinom(new_y,new_N_int,rate)
        }
      } else if (rate>threshold2){ # want orange or red
        if(!is_green(threshold1,alpha_post,beta_post)){
          total_prob_of_wanted_colour <- total_prob_of_wanted_colour +dbinom(new_y,new_N_int,rate)
        }
        
      } else { # threshold1 <= rate <= threshold2
        total_prob_of_wanted_colour <- total_prob_of_wanted_colour +dbinom(new_y,new_N_int,rate)
        # any colour is okay
        
      }
      
    }
    total_prob_of_wanted_colour
  }
  
  func_to_int <- function(new_y) {
    
    total_prob_of_wanted_colour = 0
    sd_likelihood <-  sqrt(new_N * rate * (1 - rate))
    
    alpha_post <- alpha_prior + new_y
    beta_post <- beta_prior + new_N - new_y
    
    if (rate<threshold1){
      total_prob_of_wanted_colour <- ifelse(!is_red(threshold2,alpha_post, beta_post),dnorm(new_y, new_N * rate, sd_likelihood),0)
      
    } else if (rate>threshold2){ # want orange or red
      total_prob_of_wanted_colour <- ifelse(!is_green(threshold1,alpha_post,beta_post),dnorm(new_y, new_N * rate, sd_likelihood),0)
      
      
    } else { # threshold1 <= rate <= threshold2
      total_prob_of_wanted_colour <- dnorm(new_y, new_N * rate, sd_likelihood)
      # any colour is okay
      
    }
    
    return(total_prob_of_wanted_colour)
  }
  
  
  
  if(new_N!=0 & dist_type=="normal"){
    
    full_total_prob_of_wanted_colour <- tryCatch(
      {
        integrate(func_to_int, lower = 0, upper = new_N)$value
      },
      error = function(e){
        print(e)
        print("trying binomial version instead")
        binomial_version()
      }
    )
    
  }  else if(new_N ==0 | dist_type=="binom" )    {
    
    full_total_prob_of_wanted_colour <- binomial_version()
    
  } else{
    stop("either N<0 or dist type is neither normal nor binom")
  }
  
  return(full_total_prob_of_wanted_colour)
  
}



# gives the recommended sample size N when the rate is above threshold 2 (where the 
# prior is green or orange); make sure you have a small enough rel_tol value
# we can chose either dist_type = "normal" or dist_type = "binom". The "binom"  
# method is slower and produces un-smooth plots, due to its discrete nature
recommended_sample_size <- function(threshold1, threshold2, prior_N, prior_y, rel_tol,dist_type="normal"){
  
  if (threshold1 > threshold2) {
    stop("T1>T2, hence the prior is already 'red'")
  }
  else{
    m_norm_probabilites_above_T2_internal <- memoise(norm_probabilites_above_T2_internal)
    
    second_normalisation <- function(rate_list){
      sapply(rate_list, function(x) m_norm_probabilites_above_T2_internal(x,threshold2,prior_N, prior_y))
    }
    normalisation_ans <-integrate(second_normalisation,threshold2,1,rel.tol = rel_tol)$value
    
    expected_probability_of_right_colour <- function (N){
      funct_to_integrate <- function(rate_list){
        sapply(rate_list, function(x) m_norm_probabilites_above_T2_internal(x,threshold2,prior_N, prior_y,rel_tol)/normalisation_ans*prob_correct_colour(N,x, threshold1, threshold2,prior_N,prior_y,dist_type=dist_type))
      }
      ans <-integrate(funct_to_integrate,threshold2,1,rel.tol = rel_tol)
      ans_value <- ans$value
      return(ans_value)
    }
    
    m_expected_probability_of_right_colour <- memoise(expected_probability_of_right_colour)
  }
  
  # finding a good start point for the optimiser, assuming monotonicly increasing function
  start_N = 1000
  start_N_increment = 1000
  prob_over_95 = FALSE
  previous_prob_of_colour<- 0.0
  while(!prob_over_95 ){
    
    prob_of_colour <- m_expected_probability_of_right_colour(start_N)
    if(prob_of_colour>0.95){
      prob_over_95 <- TRUE
    }
    else{
      if (previous_prob_of_colour > prob_of_colour & dist_type=="normal"){
        # something has gone wrong, try a smaller increment
        start_N = 0 
        start_N_increment <- start_N_increment/2
        
      }else{
        start_N <- start_N + start_N_increment
        
      }
      previous_prob_of_colour <- prob_of_colour
      
    }
    if(start_N> 1000000){
      stop("start_N is way too large")
    }
    if(start_N_increment<100 & dist_type=="normal"){
      # the non-monotone part has been reached
      stop("T1 is too close to T2, recommend reevaluating pathway risk")
      
    }
    else if(start_N_increment<100 & dist_type=="binom"){
      stop("start_N_increment is getting too small, and dist_type = binom :(")
    }
  }
  
  N_lb = max(0,start_N - start_N_increment)
  N_ub = N_lb + start_N_increment
  
  result<-optimise(function (x) (m_expected_probability_of_right_colour(x)-.95)^2,interval =c(N_lb, N_ub))
  if(abs(result$minimum - N_ub)<0.0001){
    stop("T1 and T2 too close to calculate properly or some other nonlinearity")
  }
  
  
  return(result)
}





