# Anne S. Warlaumont

library(lme4)
library(lmerTest)

two_agent_vocal_sim <- function(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh) {

  a1_voc_record = c() # initialize a variable that will record when agent 1 (the infant) vocalizes
  a2_voc_record = c() # initialize a variable that will record when agent 2 (the adult) vocalizes
  a2toa1_r_record = c() # initialize the record of when agent 2 responded to agent 1

  for (t in 1:sim_length){
    
    a1_p_voc_multiplier = rlnorm(1,meanlog=a1_meanlog,sdlog=a1_sdlog)
    a1_p_voc = max(a1_minp,min(a1_maxp,a1_p_voc_multiplier*a1_p_voc))
  
    a2_p_voc_multiplier = rlnorm(1,meanlog=a2_meanlog,sdlog=a2_sdlog)
    a2_p_voc = max(a2_minp,min(a2_maxp,a2_p_voc_multiplier*a2_p_voc))
    
    # Determine whether a1 vocalizes at time t
    a1_voc_event = rbinom(1,1,a1_p_voc)
    a1_voc_record = c(a1_voc_record,a1_voc_event)
    
    # Determine whether a2 vocalizes at time t, or if a2 cannot vocalize because a1 is vocalizing, increase a2's probability of vocalizing next time
    if(a1_voc_event==1){
      a2_p_voc = max(a2_minp,min(a2_maxp,a2_p_voc*a2_othersensitivity))
      a2_voc_event = 0
    } else{
      a2_voc_event = rbinom(1,1,a2_p_voc)
    }
    a2_voc_record = c(a2_voc_record,a2_voc_event)
    
    # Determine if a1 received a response from a2, and also update a1's future probability of vocalizing if so
    if (t>rthresh){
      if (a1_voc_record[t-rthresh] == 1){ # If a1 vocalized at time t-rthresh
        if (sum(a2_voc_record[(t-rthresh+1):t])>0){ # If a2 vocalized within rthresh seconds afterward
          a2toa1_r_record = c(a2toa1_r_record,1) # Record that a2 responded to a1
          a1_p_voc = max(a1_minp,min(a1_maxp,a1_p_voc*a1_respsensitivity)) # Increase a1's probability of vocalizing (unless a1_respsensitivity = 1)
        }else{ # If a2 did not vocalize within rthresh seconds afterward
          a2toa1_r_record = c(a2toa1_r_record,0) # Record that a2 did not respond to a1
        }
      }else{ # If a1 did not vocalize at time t-rthresh
        a2toa1_r_record = c(a2toa1_r_record,NA) # Record that whether a2 responded to a1 is not applicable
      }
    }

  }

  return(list(a1_voc_record,a2_voc_record,a2toa1_r_record))

}

extreme_sim <- function(sim_length,rthresh){
  
  a1_voc_record = c() # initialize a variable that will record when agent 1 (the infant) vocalizes
  a2_voc_record = c() # initialize a variable that will record when agent 2 (the adult) vocalizes
  a2toa1_r_record = c() # initialize the record of when agent 2 responded to agent 1
  
  for (t in 1:sim_length){
    
    # Determine whether a1 vocalizes at time t
    if ((t %% 20 <= 10) & (t %% 2 < 1)){
      a1_voc_event = 1
    } else if ((t %% 20 >= 10) & t%%4 == 1){
      a1_voc_event = 1
    } else {
      a1_voc_event = 0
    }
    a1_voc_record = c(a1_voc_record,a1_voc_event)
    
    # Determine whether a2 vocalizes at time t
    if ((t %% 20 < 10) & (t %% 4 == 1)){
      a2_voc_event = 1
    } else {
      a2_voc_event = 0
    }
    a2_voc_record = c(a2_voc_record,a2_voc_event)
    
    
    # Determine if a1 received a response from a2
    if (t>rthresh){
      if (a1_voc_record[t-rthresh] == 1){ # If a1 vocalized at time t-rthresh
        if (sum(a2_voc_record[(t-rthresh+1):t])>0){ # If a2 vocalized within rthresh seconds afterward
          a2toa1_r_record = c(a2toa1_r_record,1) # Record that a2 responded to a1
        }else{ # If a2 did not vocalize within rthresh seconds afterward
          a2toa1_r_record = c(a2toa1_r_record,0) # Record that a2 did not respond to a1
        }
      }else{ # If a1 did not vocalize at time t-rthresh
        a2toa1_r_record = c(a2toa1_r_record,NA) # Record that whether a2 responded to a1 is not applicable
      }
    }
    
  }
  
  return(list(a1_voc_record,a2_voc_record,a2toa1_r_record))
  
}

get_ivis_and_responses <- function(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record){
  
  sim_length = length(a1_voc_record)
  
  a1_ivi_record = c() # initialize the record of a1 inter-vocalization intervals (IVIs)
  a1_ivi_response_record = c() # initialize the record of whether those IVIs are with response or without response
  
  a1_voc1_t = NA # initialize the time of the first vocalization in the current IVI pair
  a1_voc2_t = NA # initialize the time of the second vocalization in the current IVI pair
  
  for (t in 1:sim_length){
    
    if (is.na(a1_voc1_t) & a1_voc_record[t]==1){
      
      a1_voc1_t = t
      
    } else if (a1_voc_record[t]==1){
      
      a1_voc2_t = t
      ivi = a1_voc2_t-a1_voc1_t
      a1_ivi_record = c(a1_ivi_record,ivi)
      
      if (ivi<=rthresh){
        # If the IVI is less than the response time window, then
        # we want to exclude it from analysis comparing IVIs with vs. without responses:
        a1_ivi_response_record = c(a1_ivi_response_record,NA)
      } else if (a2toa1_r_record[a1_voc1_t]==1){
        a1_ivi_response_record = c(a1_ivi_response_record,1)
      } else{
        a1_ivi_response_record = c(a1_ivi_response_record,0)
      }
      
      a1_voc1_t = a1_voc2_t
      
    }
  }
  
  return(list(a1_ivi_record,a1_ivi_response_record))
  
}

analyze_ivis <- function(ivi_records,ivi_response_records,previvi_resids,simIDs){
  
  # Compare IVI with vs. without response, without controlling for previous IVI
  #uncontrolled_response_model = lmer(scale(log(ivi_records))~ivi_response_records+(ivi_response_records|simIDs))
  uncontrolled_response_model = lm(scale(log(ivi_records))~ivi_response_records)
  
  # Compare with vs. without response residuals of the correlation between current IVI and previous IVI
  #residual_response_model = lmer(scale(previvi_resids)~ivi_response_records+(ivi_response_records|simIDs))
  residual_response_model = lm(scale(previvi_resids)~ivi_response_records)
  
  return(list(uncontrolled_response_model,residual_response_model,previvi_model))
  
}

# # Run simulations and analyze data:
# 
# sim_length = 10*60*60
# rthresh = 5 #1
# 
# for (a2_othersensitivity in 100){ #c(1,100)){
#   for (a2_respsensitivity in 1){ #c(1,100)){
#     for (a1_othersensitivity in 1){ #c(1,100)){
#       for (a1_respsensitivity in c(1,100)){ #1){
#         
#         print("Agent parameters:")
#         print(paste("agent 1 (i.e. infant) other sensitivity:",
#                     a1_othersensitivity))
#         print(paste("agent 1 response sensitivity:",
#                     a1_respsensitivity))
#         print(paste("agent 2 (i.e. adult) other sensitivity:",
#                     a2_othersensitivity))
#         print(paste("agent 2 response sensitivity:",
#                     a2_respsensitivity))
#         print("* The values indicate the factor that probability of vocalizing is multiplied by after the other agent vocalizes or responds. 1 is thus completely insensitive and higher values indicate higher degrees of sensitivity.")
#           
#         simID = 0
#         a1_ivi_records = c()
#         a1_ivi_response_records = c()
#         simIDs=c()
#         previvi_resids = c()
#         
#         for (i in 1:30){ #1){
#           
#           simID = simID+1
#           
#           a1_meanlog = 0
#           a2_meanlog = 0
#           a1_sdlog = .1 #runif(1,.01,.5)
#           a2_sdlog = .1 #runif(1,.01,.5)
#           a1_minp = .1 #runif(1,.01,.1)
#           a2_minp = .1 #runif(1,.01,.1)
#           a1_maxp = .5 #runif(1,.25,1)
#           a2_maxp = .5 #runif(1,.25,1)
#           a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
#           a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)
#           
#           voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
#           a1_voc_record = voc_records[[1]]
#           a2_voc_record = voc_records[[2]]
#           a2toa1_r_record = voc_records[[3]]
#           
#           ivi_records = get_ivis_and_responses(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record)
#           
#           a1_ivi_record = ivi_records[[1]]
#           n_ivi = length(a1_ivi_record)
#           a1_ivi_records = c(a1_ivi_records,a1_ivi_record[1:(n_ivi-1)])
#           
#           # Correlate current IVI with previous IVI and get the residuals of that correlation
#           previvi_model = lm(scale(log(a1_ivi_record[2:n_ivi]))~scale(log(a1_ivi_record[1:(n_ivi-1)])))
#           previvi_resid = resid(previvi_model)
#           previvi_resids = c(previvi_resids,previvi_resid[1:(n_ivi-1)])
#           
#           a1_ivi_response_record = ivi_records[[2]]
#           a1_ivi_response_records = c(a1_ivi_response_records,a1_ivi_response_record[2:length(a1_ivi_response_record)])
#           
#           simIDs = c(simIDs,rep(simID,n_ivi-1))
#           
#         }
#         
#         ivi_models = analyze_ivis(a1_ivi_records,a1_ivi_response_records,previvi_resids,simIDs)
#         a1_uncontrolled_response_model = ivi_models[[1]]
#         a1_residual_response_model = ivi_models[[2]]
#         a1_previvi_model = ivi_models[[3]]
#         print(summary(a1_uncontrolled_response_model))
#         print(summary(a1_residual_response_model))
#         
#       }
#     }
#   }
# }
# 
# View(data.frame(a1_voc_record,a2_voc_record,c(a2toa1_r_record,NA)))
# View(data.frame(a1_ivi_record,a1_ivi_response_record,c(NA,previvi_resid)))
# hist(a1_ivi_record)

# Run an extreme simulation:
sim_length = 10*60*60
rthresh = 1
simID = 0
a1_ivi_records = c()
a1_ivi_response_records = c()
simIDs=c()
previvi_resids = c()
for (i in 1:30){ #1){
  
  simID = simID+1
  
  voc_records = extreme_sim(sim_length,rthresh)
  a1_voc_record = voc_records[[1]]
  a2_voc_record = voc_records[[2]]
  a2toa1_r_record = voc_records[[3]]
  
  ivi_records = get_ivis_and_responses(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record)
  
  a1_ivi_record = ivi_records[[1]]
  n_ivi = length(a1_ivi_record)
  a1_ivi_records = c(a1_ivi_records,a1_ivi_record[1:(n_ivi-1)])
  
  # Correlate current IVI with previous IVI and get the residuals of that correlation
  previvi_model = lm(scale(log(a1_ivi_record[2:n_ivi]))~scale(log(a1_ivi_record[1:(n_ivi-1)])))
  previvi_resid = resid(previvi_model)
  previvi_resids = c(previvi_resids,previvi_resid[1:(n_ivi-1)])
  
  a1_ivi_response_record = ivi_records[[2]]
  a1_ivi_response_records = c(a1_ivi_response_records,a1_ivi_response_record[2:length(a1_ivi_response_record)])
  
  simIDs = c(simIDs,rep(simID,n_ivi-1))
  
}
ivi_models = analyze_ivis(a1_ivi_records,a1_ivi_response_records,previvi_resids,simIDs)
a1_uncontrolled_response_model = ivi_models[[1]]
a1_residual_response_model = ivi_models[[2]]
a1_previvi_model = ivi_models[[3]]
print(summary(a1_uncontrolled_response_model))
print(summary(a1_residual_response_model))
print(a1_voc_record)
print(a2_voc_record)
View(data.frame(a1_voc_record,a2_voc_record,c(a2toa1_r_record,NA)))
View(data.frame(a1_ivi_record,a1_ivi_response_record,c(NA,previvi_resid)))

# Run a simulation designed to make infant very bursty and not at all affected by adult and to make adult much more likely to respond when an infant makes multiple vocalizations within a burst
# I think it is possible the reason we're not getting a spurious result of significantly negative response effect when adult is sensitive to infant but infant is not sensitive to adult could be because our adult sensitivity drops off too quickly over time, i.e. is not adding up over multiple successive infant vocalizations to a sufficient extent that adult vocalizations are actually more likely within vs. outside of bursts.


# To-do:
# * Consider implementing the random lognormal drift in p_voc separately from the decaying influence of the other agent's vocalization or response
# * Figure out why there are not spurious significant effects where with-response IVIs are smaller than without-response IVIs when not controlling for previous IVI. Perhaps the response window needs to be longer or other parameters need to be different for this to show up?
# * Write the simulation data and records, or a subset of it, to a data frame to help to check that this code is computing everything as expected. I haven't done any substantive checking yet for accuracy.
# * Decouple computing response for increase in a1's voc probability from computing response for analysis purposes
# * Then test the effects of analyzing the data with different time windows, assuming different infant response sensitivities
# * Develop an analysis method that can identify sensitivity to other vocalization regardless of contingency
# * Try using a more flexible prediction of next vocalization onset, e.g. using reservoir computing for time/event series prediction

