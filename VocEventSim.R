# Anne S. Warlaumont

library(lme4)
library(lmerTest)
library(ggplot2)
#library(reservoirnet)

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

analyze_ivis <- function(ivi_records,ivi_response_records,previvi_resids,simIDs){
  
  # Compare IVI with vs. without response, without controlling for previous IVI
  #uncontrolled_response_model = lmer(scale(log(ivi_records))~ivi_response_records+(ivi_response_records|simIDs))
  uncontrolled_response_model = lm(scale(log(ivi_records))~ivi_response_records)
  
  # Compare with vs. without response residuals of the correlation between current IVI and previous IVI
  #residual_response_model = lmer(scale(previvi_resids)~ivi_response_records+(ivi_response_records|simIDs))
  residual_response_model = lm(scale(previvi_resids)~ivi_response_records)
  
  return(list(uncontrolled_response_model,residual_response_model,previvi_model))
  
}

# Run simulations and analyze data:

simcounter = 0
voc_and_resp_records = list()
ivi_records = list()

sink(file="VocEventSimOutput.txt")
sim_length = 10*60*60

for (rthresh in c(1,5)){
  for (a2_othersensitivity in c(1,1.5,2,3,100)){
    for (a2_respsensitivity in c(1,1.5,2,3,100)){
      for (a1_othersensitivity in c(1,1.5,2,3,100)){
        for (a1_respsensitivity in c(1,1.5,2,3,100)){

          simcounter = simcounter+1

          print(paste("Simulation number:",simcounter))
          print("Simulation parameters:")
          print(paste("response threshold:",rthresh))
          print(paste("agent 1 (i.e. infant) other sensitivity:",
                      a1_othersensitivity))
          print(paste("agent 1 response sensitivity:",
                      a1_respsensitivity))
          print(paste("agent 2 (i.e. adult) other sensitivity:",
                      a2_othersensitivity))
          print(paste("agent 2 response sensitivity:",
                      a2_respsensitivity))
          print("* The values indicate the factor that probability of vocalizing is multiplied by after the other agent vocalizes or responds. 1 is thus completely insensitive and higher values indicate higher degrees of sensitivity.")

          simID = 0
          a1_ivi_records = c()
          a1_ivi_response_records = c()
          simIDs=c()
          previvi_resids = c()

          for (i in 1:5){ #200){

            simID = simID+1

            a1_meanlog = 0
            a2_meanlog = 0
            a1_sdlog = .2
            a2_sdlog = .2
            a1_minp = .000001
            a2_minp = .1
            a1_maxp = .4
            a2_maxp = .95
            a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
            a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)

            voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
            a1_voc_record = voc_records[[1]]
            a2_voc_record = voc_records[[2]]
            a2toa1_r_record = voc_records[[3]]
            a2toa2_r_record = voc_records[[4]]

            ivi_records = get_ivis_and_responses(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record) # This is not yet updated to handle a1toa2_r_record

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
          #a1_uncontrolled_response_model_summary = summary(a1_uncontrolled_response_model)
          print(summary(a1_uncontrolled_response_model))
          print(summary(a1_residual_response_model))
          hist(a1_ivi_record)

          voc_and_resp_records = c(voc_and_resp_records,list(data.frame(a1_voc_record,a2_voc_record,c(a2toa1_r_record,rep(NA,times=rthresh)))))
          ivi_records = c(ivi_records,list(data.frame(a1_ivi_record,a1_ivi_response_record,c(NA,previvi_resid))))
        }
      }
    }
  }
}

# # Run a single simulation and visualize the data.
# sim_length = 10*60*60
# rthresh = 1
# a2_othersensitivity = 1
# a2_respsensitivity = 1
# a1_othersensitivity = 1
# a1_respsensitivity = 1
# # a1_meanlog = 0
# # a2_meanlog = 0
# # a1_sdlog = .2
# # a2_sdlog = .2
# # a1_minp = .000001
# # a2_minp = .000001
# # a1_maxp = .4
# # a2_maxp = .4
# a1_meanlog = 0
# a2_meanlog = 0
# a1_sdlog = .2
# a2_sdlog = .2
# a1_minp = .000001
# a2_minp = .1
# a1_maxp = .4
# a2_maxp = .95
# a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
# a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)
# 
# voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
# a1_voc_record = voc_records[[1]]
# a2_voc_record = voc_records[[2]]
# a2toa1_r_record = voc_records[[3]]
# a2toa2_r_record = voc_records[[4]]
# 
# ivi_records = get_ivis_and_responses(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record) # This is not yet updated to handle a1toa2_r_record
# 
# a1_ivi_record = ivi_records[[1]]
# n_ivi = length(a1_ivi_record)
# 
# hist(a1_ivi_record)
# ggplot(data.frame(a1_ivi_record), aes(x = a1_ivi_record)) +
#   geom_histogram(bins = 30) +
#   scale_x_log10() +
#   scale_y_log10()
# 
# barcode_data = data.frame(
#   index = seq_along(a1_voc_record),
#   value = a1_voc_record
# )
# ggplot(barcode_data, aes(x = index, xend = index, y = 0, yend = 1, color = factor(value))) +
#   geom_segment(size = 5) +
#   scale_color_manual(values = c("white","black")) +
#   labs("Agent 1 vocalization times",
#        x = "Time (s)",
#        y = "") +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# #pdf(file="SingleSim_MultiscaleClusters.pdf")
# 
# par(mfrow=c(3,1),cex=1,mar=c(1,1,2,1),oma=c(0,0,4,0))
# stripchart(which(a1_voc_record==1),xaxt="n",main="Full day simulation (10 hours)",pch=19,ylim=c(.5,1.5),xlim=c(0,sim_length))
# 
# hr_offset = 2
# rect(xleft=3600+hr_offset*60*60,xright=7200+hr_offset*60*60,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
# a1_voc_record_1hr = a1_voc_record[(3601+hr_offset*60*60):(7200+hr_offset*60*60)]
# stripchart(which(a1_voc_record_1hr==1),xaxt="n",main="1 hour within the day",pch=19,ylim=c(.5,1.5),xlim=c(0,3600))
# 
# fivemin_offset = 5
# rect(xleft=0+fivemin_offset*60*5,xright=300+fivemin_offset*60*5,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
# a1_voc_record_5min = a1_voc_record_1hr[(1+fivemin_offset*60*5):(300+fivemin_offset*60*5)]
# stripchart(which(a1_voc_record_5min==1),xaxt="n",main="5 minutes within the hour",pch=19,ylim=c(.5,1.5),xlim=c(0,300))
# 
# mtext("Onsets of child vocalizations",side=3, line = 1, outer=TRUE, cex=2)
# 
# #dev.off()
# 
# sum(a1_voc_record)
# sum(a2_voc_record)


# To-do:
# * Add code to create visualizations of vocalization probability and vocalization events over time.
# * Write the simulation data and records, or a subset of it, to a data frame to help to check that this code is computing everything as expected. I haven't done any substantive checking yet for accuracy.
# * Decouple computing response for increase in a1's voc probability from computing response for analysis purposes. Then test the effects of analyzing the data with different time windows, assuming different infant response sensitivities
# * Develop an analysis method that can identify sensitivity to other vocalization regardless of contingency
# * Try using a more flexible prediction of next vocalization onset, e.g. using reservoir computing for time/event series prediction

