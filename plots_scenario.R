##############################################################################
#
# plots_scenario.R
#
# Thao P. Le, Thomas K. Waring, Chris M. Baker
#
# 2023
# 
# Plots scenarios using our sampling method (Figures 4, 5, 6)
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


simulate_quarters <- function(N_prior_list,y_prior_list,num_quarters,rates,threshold1_at_prior_level, threshold2, take_previous_method = list("number",10000), sampling_method=list("fixed",1875), verbose=FALSE ){
  
  pre_quarters = length(N_prior_list)
  
  N_t <- rep(c(NA),pre_quarters+num_quarters)
  y_t <- rep(c(NA),pre_quarters+num_quarters)
  
  for(i in seq(1,length(N_prior_list))){
    # first add in the prior samples
    N_t[i]<-N_prior_list[i]
    y_t[i]<-y_prior_list[i]
  }
  
  all_samples <- c()
  # generate fake samples that could correspond to the prior data
  for(i in seq(1,length(N_prior_list))){
    correct_number_of_simulated_prior_leakages = FALSE
    while(!correct_number_of_simulated_prior_leakages){
      simulated_samples <- rbinom(N_prior_list[i],1,y_prior_list[i]/N_prior_list[i])
      simulated_prior_leakages <- sum(simulated_samples)
      if (simulated_prior_leakages==y_prior_list[i]){
        correct_number_of_simulated_prior_leakages = TRUE
        all_samples <- c(all_samples,simulated_samples)
      }
    }
  }
  
  
  
  prior_N_t <- rep(c(NA),num_quarters)
  prior_y_t <- rep(c(NA),num_quarters)
  colour_t <- rep(c(NA),num_quarters)
  a_post_t <- rep(c(NA),num_quarters)
  b_post_t <- rep(c(NA),num_quarters)
  
  
  for (reporting in 1:(num_quarters)){
    
    print(paste0("reporting period: ", reporting))
    
    # first calculate the prior
    if(take_previous_method[[1]]=="number"){
      N_prior <- take_previous_method[[2]]
      if (length(all_samples)<N_prior){
        stop("not enough prior samples given the take_previous_method!")
      }
      
      y_prior <- sum(tail(all_samples,N_prior))
      
      prior_N_t[reporting] <- N_prior
      prior_y_t[reporting] <- y_prior
      
    } else if (take_previous_method[[1]]=="quarters"){
      index_shift = pre_quarters
      quarters_to_take <- seq(reporting-take_previous_method[[2]]+index_shift, reporting-1+index_shift,1)
      
      
      N_prior <- sum(N_t[quarters_to_take])
      y_prior <- sum(y_t[quarters_to_take])
      
      prior_N_t[reporting] <- N_prior
      prior_y_t[reporting] <- y_prior
      
    } else {
      
      stop("wrong take_previous_method!")
      
    }
    
    prior_a <- 0.5 + y_prior
    prior_b <- 0.5 + N_prior - y_prior
    
    threshold1 <- calc_threshold_1(threshold1_at_prior_level, prior_a, prior_b)
    
    prior_colour_status<-colour_status(threshold1, threshold2, prior_a, prior_b)
    
    
    if(threshold1>threshold2){
      print("prior is red, stop code")
      # prior is now red! stop code!
      num_quarters = reporting - 1
      break
    }
    
    rel_tol = 0.0001
    get_fluid_N <- function(){
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
    
    
    # get the next sample size
    
    if(sampling_method[[1]]=="fixed"){
      
      N<-sampling_method[[2]]
      
    } else if (sampling_method[[1]]=="fluid"){
      
      N<-get_fluid_N()
      
    } else if (sampling_method[[1]]=="fixed-with-targeted"){
      
      fixed_capacity<-sampling_method[[2]]
      target_capacity<-sampling_method[[4]]
      
      if (sampling_method[[3]]=="orange"){
        
        if (prior_colour_status=="orange"){
          #when the posterior goes orange and the fluid sampling tells us to sample more, and so more is sampled, up the recommended minimum of the fluid sampling 
          
          fluid_N<-get_fluid_N()
          if (fluid_N>(fixed_capacity)){
            N <- min(fluid_N,fixed_capacity+target_capacity)
          } else {
            N <- fixed_capacity
          }
        } else {
          N <- fixed_capacity
        }
        
      } else if(sampling_method[[3]]=="date"){
        
        quarter_to_start_targeting<-sampling_method[[5]]
        # when secondary sources tell us that we might want to increase sampling in the next time period
        
        if (reporting>= quarter_to_start_targeting){
          N <- fixed_capacity + target_capacity
        } else {
          N <- fixed_capacity
        }
        
      } else if(sampling_method[[3]]=="matching"){
        
        # when the fixed sampling is say 1/2 of the fluid sampling, apply the targeted capacity (regardless of colour) 
        when_this_smaller_than_fluid_fraction <- sampling_method[[5]]
        
        fluid_N<-get_fluid_N()
        
        if (fixed_capacity < when_this_smaller_than_fluid_fraction*fluid_N){
          N <- min(fluid_N,fixed_capacity+target_capacity)
        } else {
          N <- fixed_capacity
        }
        
      } else {
        stop("wrong kind of targeted sampling method!")
      }
      
    } else{
      
      stop("wrong kind of sampling method!")
      
    }
    
    print("N:")
    print(N)
    
    if (N < 0) {
      num_quarters = reporting - 1
      break
    }
    
    
    # now onto sampling!
    N_t[pre_quarters+reporting] <- N
    
    rate <- rates[reporting]
    samples <- rbernoulli(N, rate)
    all_samples <- c(all_samples,samples)
    
    y <- sum(samples)
    y_t[pre_quarters+reporting] <- y
    
    
    a <- 0.5 + y_prior + y
    a_post_t[reporting] <- a
    b <- 0.5 + N_prior + N - y_prior - y
    b_post_t[reporting] <- b
    
    colour_t[reporting] <- colour_status(threshold1, threshold2, a,b)
    
    
  }
  print(paste0("num quarters: ",num_quarters))
  
  simulation_quarters = seq(pre_quarters+1,pre_quarters+num_quarters,1)
  
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


plot_sim <- function (sim_df,name,pathway_name="pathway",additional_save_name="prior_quarters_fixed_sampling",fluid_or_fixed="fixed") {
  sim_df['leakage_proportion'] <- 100* sim_df['y_t'] / sim_df['N_t']
  
  print(name)
  print(sim_df)
  
  colour_key = c("green"="#cff983", "orange"="#ff9448", "red"="#C70039")
  
  colour_key_border = c("green"="green3", "orange"="darkorange", "red"="red2")
  
  
  
  max_sampling = max(sim_df['N_t'])
  x_lim_max = max(c(3000,max_sampling))
  
  factor <- x_lim_max/4
  
  sim_df['leakage_proportion'] <- 100*factor* sim_df['y_t'] / sim_df['N_t']
  
  ggplot(data =sim_df) +
    geom_bar(aes(x=quarter_nums, y=N_t, fill=colour_t,color = colour_t),stat='identity') +
    scale_fill_manual(values = colour_key, name='Status') +
    scale_color_manual(values = colour_key_border, name='Status') +
    geom_line(aes(x=as.numeric(quarter_nums), y=as.numeric(leakage_proportion)),stat="identity",color="black",size=1.5)+
    geom_point(aes(x=as.numeric(quarter_nums), y=as.numeric(leakage_proportion)),shape=21, color="black", fill="#69b3a2", size=3)+
    xlab('Quarters') + ylab('Number of samples') + 
    scale_y_continuous(sec.axis=sec_axis(~./factor,name="Proportion of leakages (%)",labels = function(x) paste0(x, "%")),limits=c(0,x_lim_max))+
    theme_light()+
    theme( legend.position=c(0.24,0.83),legend.background = element_rect(fill = "white", color = "grey"),text = element_text(size = 17))
  
  print(paste0("outputs/",pathway_name,"_",name,"_",additional_save_name, "_combined.png"))
  
  ggsave(paste0("outputs/",pathway_name,"_",name,"_",additional_save_name, "_combined.png"),width = 4, height = 4)
  
  
}


################################################################################
# Plotting
################################################################################

########## parameters #######################################


prior_N <- 10000
prior_y <- 6

N_prior_list <- c(5000,5000)
y_prior_list <- c(3,3)


num_quarters<- 5
low_rate <- 0.0012
high_rate <- 0.02
threshold1_at_prior_level <- 0.95
threshold2<-0.005

pathway_name<-"Scenario_10k_6"

take_previous_method <- list("quarters", 2)

routine_save_name<- paste0(pathway_name,"-",threshold1_at_prior_level,"-",threshold2 ,"_", low_rate, "_")
highrisk_save_name<- paste0(pathway_name,"-",threshold1_at_prior_level,"-",threshold2 ,"_", high_rate, "_")


routine_rates <- rep(low_rate, num_quarters)

num_routine <- 4
highrisk_rates <-c(rep(low_rate,num_routine),rep(high_rate,num_quarters-num_routine)) git
git
very_low_rate <- 0.0001
very_low_rates <- rep(very_low_rate,num_quarters)
very_low_save_name<- paste0(pathway_name,"-",threshold1_at_prior_level,"-",threshold2 ,"_", very_low_rate, "_")

semi_high_rate <- 0.005
slowly_rising_rates <- seq(low_rate,semi_high_rate,(semi_high_rate-low_rate)/(num_quarters-1))
slowly_rising_save_name<- paste0(pathway_name,"-",threshold1_at_prior_level,"-",threshold2 ,"_slow", low_rate,'-',semi_high_rate , "_")


########## sampling #######################################


sampling_method = list("fluid")
additional_save_name <- paste0("prior_",take_previous_method[[1]],"_",sampling_method[[1]])

simulation_routine <- simulate_quarters(N_prior_list,y_prior_list,num_quarters,routine_rates,threshold1_at_prior_level, threshold2, take_previous_method , sampling_method)

simulation_red <- simulate_quarters(N_prior_list,y_prior_list,num_quarters,highrisk_rates,threshold1_at_prior_level, threshold2, take_previous_method , sampling_method)

simulation_routine_low <- simulate_quarters(N_prior_list,y_prior_list,num_quarters,very_low_rates,threshold1_at_prior_level, threshold2, take_previous_method , sampling_method)

simulation_routine_rising <- simulate_quarters(N_prior_list,y_prior_list,num_quarters,slowly_rising_rates,threshold1_at_prior_level, threshold2, take_previous_method , sampling_method)

save(simulation_routine, simulation_red, simulation_routine_low,simulation_routine_rising, file = "outputs/plots_scenario_saved_outputs.rdata")


load("outputs/plots_scenario_saved_outputs.rdata")

plot_sim(simulation_routine,"routine",routine_save_name,additional_save_name,sampling_method[[1]])
plot_sim(simulation_red,"risky",highrisk_save_name,additional_save_name,sampling_method[[1]])
plot_sim(simulation_routine_low,"very_low_risk",very_low_save_name,additional_save_name,sampling_method[[1]])
plot_sim(simulation_routine_rising,"slowly_rising",slowly_rising_save_name,additional_save_name,sampling_method[[1]])
