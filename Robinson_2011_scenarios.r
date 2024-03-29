##############################################################################
#
# Robinson_2011_scenarios
#
# 2023
# 
# Here we run scenarios using the sampling picking method of Robinson 2011
# 
##############################################################################
library(tidyverse)
library(tidyr)
source("allocate.r")

##############################################################################

# it will be red if the *future* sample size is 100% 
colour_status_robinson <- function(sample_size_upper_limit,sample_size_result){
  if(sample_size_result<sample_size_upper_limit){
    return("green")
  }
  else{
    return("red")
  }
  
}

simulate_quarters_robinson <- function(N_prior_list, y_prior_list,num_quarters,rates, sample_size_upper_limit,risk_cutoff, take_previous_method = list("number",10000),verbose=FALSE){
  
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
    
    initial_samplesize <- N_prior
    initial_leakage <- y_prior 
    
    # get the next sample size
    N <-  inspect(initial_leakage, initial_samplesize, sample_size_upper_limit, risk.cutoff = risk_cutoff,verbose = TRUE)$sample.size
    
    print("N:")
    print(N)
    
    # there is probably no clear break condition here
    if (N==sample_size_upper_limit){
      print("Sampling 100%, should this be a break condition?")
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
    
    if(reporting>1){
      colour_t[reporting-1] <- colour_status_robinson(sample_size_upper_limit,N)
    }
    
    
    if(reporting==num_quarters){
      # need to calculate the next time period too
      if(take_previous_method[[1]]=="number"){
        N_prior <- take_previous_method[[2]]
        if (length(all_samples)<N_prior){
          stop("not enough prior samples given the take_previous_method!")
        }
        
        y_prior <- sum(tail(all_samples,N_prior))
        
      } else if (take_previous_method[[1]]=="quarters"){
        index_shift = pre_quarters
        quarters_to_take <- seq((reporting+1)-take_previous_method[[2]]+index_shift, reporting+index_shift,1)
        
        
        N_prior <- sum(N_t[quarters_to_take])
        y_prior <- sum(y_t[quarters_to_take])
        
        
      } 
      
      initial_samplesize <- N_prior
      initial_leakage <- y_prior 
      
      # get the next sample size
      N <-  inspect(initial_leakage, initial_samplesize, sample_size_upper_limit, risk.cutoff = risk_cutoff,verbose = TRUE)$sample.size
      
      colour_t[reporting] <- colour_status_robinson(sample_size_upper_limit,N)
      
    }
    
    
    
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
  
  
# this is slightly different fromthe plot_scenarios.R one
plot_sim_robinson <- function (sim_df,name,pathway_name="pathway",additional_save_name="prior_quarters_fixed_sampling",fluid_or_fixed="fixed") {
  sim_df['leakage_proportion'] <- 100* sim_df['y_t'] / sim_df['N_t']

  print(name)
  print(sim_df)
  
  colour_key = c("green"="#7d7d7d", "orange"="#7d7d7d", "red"="#7d7d7d")
  colour_key = c("green"="#88b9ec", "orange"="#88b9ec", "red"="#88b9ec")
  
  max_sampling = max(sim_df['N_t'])
  x_lim_max = max(c(3000,max_sampling))
  
  factor <- x_lim_max/4
  
  sim_df['leakage_proportion'] <- 100*factor* sim_df['y_t'] / sim_df['N_t']

  ggplot(data =sim_df) +
    geom_bar(aes(x=quarter_nums, y=N_t, fill=colour_t),stat='identity',show.legend = FALSE) +
    scale_fill_manual(values = colour_key, name='Colour status') +
    geom_line(aes(x=as.numeric(quarter_nums), y=as.numeric(leakage_proportion)),stat="identity",color="black",size=1.5)+
    geom_point(aes(x=as.numeric(quarter_nums), y=as.numeric(leakage_proportion)),shape=21, color="black", fill="#69b3a2", size=3)+
    xlab('Quarters') + ylab('Number of samples') +
    scale_y_continuous(sec.axis=sec_axis(~./factor,name="Proportion of leakages (%)",labels = function(x) paste0(x, "%")),limits=c(0,x_lim_max))+
    #ggtitle(paste0("Simulation: ", name, " with ", fluid_or_fixed , " sampling"))+
     theme_light()+
    theme(text = element_text(size = 17))

  ggsave(paste0("outputs/",pathway_name,"_",name,"_",additional_save_name, "_combined.png"),width = 4, height = 4)

}



##############################################################################
# plotting
##############################################################################

prior_N <- 10000
prior_y <- 6

N_prior_list <- c(5000,5000)
y_prior_list <- c(3,3)


num_quarters<- 5
low_rate <- 0.0012
high_rate <- 0.02
risk_cutoff<-0.005
sample_size_upper_limit <- 50000

pathway_name<-"Robinson_2011_Scenario_10k_6"

take_previous_method <- list("quarters", 2)

routine_save_name<- paste0(pathway_name,"-",risk_cutoff ,"_",sample_size_upper_limit,"_", low_rate, "_")
highrisk_save_name<- paste0(pathway_name,"-",risk_cutoff ,"_",sample_size_upper_limit,"_",  high_rate, "_")

routine_rates <- rep(low_rate, num_quarters)

num_routine <- 4
highrisk_rates <-c(rep(low_rate,num_routine),rep(high_rate,num_quarters-num_routine)) 

very_low_rate <- 0.0001
very_low_rates <- rep(very_low_rate,num_quarters)
very_low_save_name<- paste0(pathway_name,"-",threshold1_at_prior_level,"-",threshold2 ,"_", very_low_rate, "_")

semi_high_rate <- 0.005
slowly_rising_rates <- seq(low_rate,semi_high_rate,(semi_high_rate-low_rate)/(num_quarters-1))
slowly_rising_save_name<- paste0(pathway_name,"-",threshold1_at_prior_level,"-",threshold2 ,"_slow", low_rate,'-',semi_high_rate , "_")

############################### fluid sampling
print("=========== fluid sampling ============")

sampling_method = list("fluid")
additional_save_name <- paste0("prior_",take_previous_method[[1]],"_",sampling_method[[1]])

simulation_routine_robinson <- simulate_quarters_robinson(N_prior_list,y_prior_list,num_quarters,routine_rates,sample_size_upper_limit,risk_cutoff,  take_previous_method )
simulation_red_robinson <- simulate_quarters_robinson (N_prior_list,y_prior_list,num_quarters,highrisk_rates,sample_size_upper_limit,risk_cutoff, take_previous_method)
simulation_routine_low_robinson <- simulate_quarters_robinson(N_prior_list,y_prior_list,num_quarters,very_low_rates,sample_size_upper_limit,risk_cutoff,  take_previous_method )
simulation_routine_rising_robinson <- simulate_quarters_robinson(N_prior_list,y_prior_list,num_quarters,slowly_rising_rates,sample_size_upper_limit,risk_cutoff,  take_previous_method )


save(simulation_routine_robinson, simulation_red_robinson, simulation_routine_low_robinson,simulation_routine_rising_robinson, file = "outputs/Robinson_2011_scenarios_saved_outputs.rdata")

load("outputs/Robinson_2011_scenarios_saved_outputs.rdata")

plot_sim_robinson(simulation_routine_robinson,"routine",routine_save_name,additional_save_name,sampling_method[[1]])
plot_sim_robinson(simulation_red_robinson,"risky",highrisk_save_name,additional_save_name,sampling_method[[1]])
plot_sim_robinson(simulation_routine_low_robinson,"very_low_risk",very_low_save_name,additional_save_name,sampling_method[[1]])
plot_sim_robinson(simulation_routine_rising_robinson,"slowly_rising",slowly_rising_save_name,additional_save_name,sampling_method[[1]])


