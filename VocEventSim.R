# Anne S. Warlaumont

library(lme4)
library(lmerTest)
library(ggplot2)

two_agent_vocal_sim <- function(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh) {
  
  a1_voc_record = c() # initialize a variable that will record when agent 1 (the infant) vocalizes
  a2_voc_record = c() # initialize a variable that will record when agent 2 (the adult) vocalizes
  a2toa1_r_record = c() # initialize the record of when agent 2 responded to agent 1
  a1toa2_r_record = c() # initialize the record of when agent 1 responded to agent 2
  
  for (t in 1:sim_length){
    
    a1_p_voc_multiplier = rlnorm(1,meanlog=a1_meanlog,sdlog=a1_sdlog)
    a1_p_voc = max(a1_minp,min(a1_maxp,a1_p_voc_multiplier*a1_p_voc))
    
    a2_p_voc_multiplier = rlnorm(1,meanlog=a2_meanlog,sdlog=a2_sdlog)
    a2_p_voc = max(a2_minp,min(a2_maxp,a2_p_voc_multiplier*a2_p_voc))
    
    # Determine whether a1 vocalizes at time t
    a1_voc_event = rbinom(1,1,a1_p_voc)
    a1_voc_record = c(a1_voc_record,a1_voc_event)
    
    # Determine whether a2 vocalizes at time t
    # If so, change a1's probability of vocalizing according to a1_othersensitivity
    # Or if a2 cannot vocalize because a1 is vocalizing, change a2's probability of vocalizing according to a2_othersensitivity
    if(a1_voc_event==1){
      a2_p_voc = max(a2_minp,min(a2_maxp,a2_p_voc*a2_othersensitivity))
      a2_voc_event = 0
    } else{
      a2_voc_event = rbinom(1,1,a2_p_voc)
      if (a2_voc_event==1){
        a1_p_voc = max(a1_minp,min(a1_maxp,a1_p_voc*a1_othersensitivity))
      }
    }
    a2_voc_record = c(a2_voc_record,a2_voc_event)
    
    # Determine if a1 received a response from a2, and if so update a1's probability of vocalizing according to a1_respsensitivity
    if (t>rthresh){
      if (a1_voc_record[t-rthresh] == 1){ # If a1 vocalized at time t-rthresh
        if (sum(a2_voc_record[(t-rthresh+1):t])>0){ # If a2 vocalized within rthresh seconds afterward
          a2toa1_r_record = c(a2toa1_r_record,1) # Record that a2 responded to a1
          a1_p_voc = max(a1_minp,min(a1_maxp,a1_p_voc*a1_respsensitivity)) # Update a1's probability of vocalizing
        }else{ # If a2 did not vocalize within rthresh seconds afterward
          a2toa1_r_record = c(a2toa1_r_record,0) # Record that a2 did not respond to a1
        }
      }else{ # If a1 did not vocalize at time t-rthresh
        a2toa1_r_record = c(a2toa1_r_record,NA) # Record that whether a2 responded to a1 is not applicable
      }
    }
    
    # Determine if a2 received a response from a1, and if so update a2's probability of vocalizing according to a2_respsensitivity
    if (t>rthresh){
      if (a2_voc_record[t-rthresh] == 1){ # If a2 vocalized at time t-rthresh
        if (sum(a1_voc_record[(t-rthresh+1):t])>0){ # If a1 vocalized within rthresh seconds afterward
          a1toa2_r_record = c(a1toa2_r_record,1) # Record that a1 responded to a2
          a2_p_voc = max(a2_minp,min(a2_maxp,a2_p_voc*a2_respsensitivity)) # Update a2's probability of vocalizing
        }else{ # If a1 did not vocalize within rthresh seconds afterward
          a1toa2_r_record = c(a1toa2_r_record,0) # Record that a1 did not respond to a2
        }
      }else{ # If a2 did not vocalize at time t-rthresh
        a1toa2_r_record = c(a1toa2_r_record,NA) # Record that whether a1 responded to a2 is not applicable
      }
    }
    
  }
  
  return(list(a1_voc_record,a2_voc_record,a2toa1_r_record,a1toa2_r_record))
  
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

analyze_ivis <- function(ivi_records,ivi_response_records,previvi_resids,prev3ivi_resids,simIDs){
  
  # Compare IVI with vs. without response, without controlling for previous IVI
  uncontrolled_response_model = lmer(scale(log(ivi_records))~ivi_response_records+(ivi_response_records|simIDs))
  #uncontrolled_response_model = lm(scale(log(ivi_records))~ivi_response_records)
  
  # Compare with vs. without response residuals of the correlation between current IVI and previous IVI
  residual_response_model = lmer(scale(previvi_resids)~ivi_response_records+(ivi_response_records|simIDs))
  #residual_response_model = lm(scale(previvi_resids)~ivi_response_records)
  
  # Compare with vs. without response residuals of the correlation between current IVI and previous 3 IVIs
  prev3residual_response_model = lmer(scale(prev3ivi_resids)~ivi_response_records+(ivi_response_records|simIDs))
  #prev3residual_response_model = lm(scale(prev3ivi_resids)~ivi_response_records)
  
  return(list(uncontrolled_response_model,residual_response_model,prev3residual_response_model))
  
}

analyze_timessince_logistic <- function(timessince_df){
  
  timessince_a1_model = glm(a1_voc_record ~ scale(log(timessince_last_a1)) + scale(log(timessince_2ndToLast_a1)) + scale(log(timessince_3rdToLast_a1)), family = "binomial", data = timessince_df, na.action = na.exclude)
  timessince_df$a1_resids = resid(timessince_a1_model)
  timessince_a2_model = lm(a1_resids ~ scale(log(timessince_last_a2)) + scale(log(timessince_2ndToLast_a2)) + scale(log(timessince_3rdToLast_a2)), data = timessince_df, na.action = na.exclude)
  timessince_df$a2_resids = resid(timessince_a2_model)
  if ((any(!is.na(timessince_df$timessince_last_a2toa1_r))) && (any(!is.na(timessince_df$timessince_2ndToLast_a2toa1_r))) && (any(!is.na(timessince_df$timessince_3rdToLast_a2toa1_r)))){
    timessince_r_model = lm(a2_resids ~ scale(log(timessince_last_a2toa1_r)) + scale(log(timessince_2ndToLast_a2toa1_r)) + scale(log(timessince_3rdToLast_a2toa1_r)), data = timessince_df)
    timessince_simultaneous_model = glm(a1_voc_record ~ scale(log(timessince_last_a1)) + scale(log(timessince_2ndToLast_a1)) + scale(log(timessince_3rdToLast_a1)) + scale(log(timessince_last_a2)) + scale(log(timessince_2ndToLast_a2)) + scale(log(timessince_3rdToLast_a2)) + scale(log(timessince_last_a2toa1_r)) + scale(log(timessince_2ndToLast_a2toa1_r)) + scale(log(timessince_3rdToLast_a2toa1_r)), family = "binomial", data = timessince_df, na.action = na.exclude)
    return(list(timessince_a1_model,timessince_a2_model,timessince_r_model,timessince_simultaneous_model))
  } else{
    return(list(timessince_a1_model,timessince_a2_model))
  }
  
}

