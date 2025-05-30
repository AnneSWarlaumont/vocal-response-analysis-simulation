setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

#recordings = c("0344_000913","0833_010606","0054_000603")
recordings = c("0196_000902")
simTypes = c("nonInteractive","a2interactive","bidirectional")

#nSims = 3 # for debugging
nSims = 10000 # for an actual batch


for (recording in recordings){
  for (simType in simTypes){
    
    ######################################
    # Read in and format the human sample
    ######################################
    
    recording_dir = paste("data/",recording,sep="")
    lena_segments = read.csv(paste(recording_dir,"/",recording,"_segments.csv",sep=""))
    chn_segments = subset(lena_segments,segtype=="CHNSP")
    adu_segments = subset(lena_segments,segtype=="FAN"|segtype=="MAN")
    rec_length = floor(lena_segments$endsec[nrow(lena_segments)])
    chn_voc_record = integer(rec_length)
    chn_voc_record[as.integer(floor(chn_segments$startsec))] = 1
    adu_voc_record = integer(rec_length)
    adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1
    
    ####################################
    # Get key measures for human sample 
    ####################################
    
    # Number of child and adult vocalizations
    chn_total = sum(chn_voc_record)
    adu_total = sum(adu_voc_record)
    
    # Number of conversational turns based on a 5 s threshold
    # Initialize leader ID and time
    # Loop through each second of the recording
    # When a turn is identified, update the turn count
    turnCount = 0
    leadID = 'none'
    leadT = 0
    for (t in 1:rec_length){
      leadT = leadT+1
      if (leadID=='chn' && adu_voc_record[t]==1) {
        if (leadT <= 5) {
          turnCount = turnCount+1
        }
        leadID = 'adu'
        leadT = 0
      } else if (leadID=='adu' && chn_voc_record[t]==1) {
        if (leadT <= 5) {
          turnCount = turnCount+1
        }
        leadID = 'chn'
        leadT = 0
      } else if (leadID=='none' && adu_voc_record[t]==1) {
        leadID = 'adu'
        leadT = 0
      } else if (leadID=='none' && chn_voc_record[t]==1) {
        leadID = 'chn'
        leadT = 0
      }
    }
    
    # 10th, 50th, and 90th percentiles of the child and adult ivi distributions
    # try the quantile function
    get_ivis <- function(voc_record){
      ivi_record = c()
      voc1_t = NA # initialize the time of the first vocalization in the current IVI pair
      voc2_t = NA # initialize the time of the second vocalization in the current IVI pair 
      for (t in 1:length(voc_record)){
        if (is.na(voc1_t) & voc_record[t]==1){
          voc1_t = t
        } else if (voc_record[t]==1){
          voc2_t = t
          ivi = voc2_t - voc1_t
          ivi_record = c(ivi_record,ivi)
          voc1_t = voc2_t
        }
      }
      return(ivi_record)
    }
    
    chn_ivi_record = get_ivis(chn_voc_record)
    adu_ivi_record = get_ivis(adu_voc_record)
    chn_ivi_10 = quantile(chn_ivi_record,probs=.1,names=FALSE)
    chn_ivi_50 = quantile(chn_ivi_record,probs=.5,names=FALSE)
    chn_ivi_90 = quantile(chn_ivi_record,probs=.9,names=FALSE)
    chn_ivi_99 = quantile(chn_ivi_record,probs=.99,names=FALSE)
    adu_ivi_10 = quantile(adu_ivi_record,probs=.1,names=FALSE)
    adu_ivi_50 = quantile(adu_ivi_record,probs=.5,names=FALSE)
    adu_ivi_90 = quantile(adu_ivi_record,probs=.9,names=FALSE)
    adu_ivi_99 = quantile(adu_ivi_record,probs=.99,names=FALSE)
    
    #########################################################################
    # Randomly choose simulation parameters and run simulation
    #########################################################################
    
    source('VocEventSim.R')
    
    sim_length = length(chn_voc_record)
    
    rthresh = 5
    
    meanlog_range = c(-.05,0)
    minp_range = c(.000000001,.00001)
    maxp_range = c(.05,.5)
    sdlog_range = c(.01,1)
    othersensitivity_range = c(1,10)
    
    sims_df = data.frame(chn_sim_minp = double(),
                         chn_sim_maxp = double(),
                         chn_sim_meanlog = double(),
                         chn_sim_sdlog = double(),
                         chn_sim_othersensitivity = double(),
                         chn_sim_respsensitivity = double(),
                         chn_sim_total = double(),
                         chn_sim_n_ivi_1 = double(),
                         chn_sim_ivi_50 = double(),
                         chn_sim_ivi_max = double(),
                         adu_sim_minp = double(),
                         adu_sim_maxp = double(),
                         adu_sim_meanlog = double(),
                         adu_sim_sdlog = double(),
                         adu_sim_othersensitivity = double(),
                         adu_sim_respsensitivity = double(),
                         adu_sim_total = double(),
                         adu_sim_n_ivi_1 = double(),
                         adu_sim_ivi_50 = double(),
                         adu_sim_ivi_max = double(),
                         sim_turncount = double(),
                         stringsAsFactors = FALSE)
    sims_chn_voc_records <- list()
    sims_adu_voc_records <- list()
    sims_chn_ivi_records <- list()
    sims_adu_ivi_records <- list()
    
    for (simNum in 1:nSims){
      
      a1_minp = runif(1,min=minp_range[1],max=minp_range[2])
      a2_minp = runif(1,min=minp_range[1],max=minp_range[2])
      
      a1_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
      a2_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
      
      a1_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
      a2_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
      
      if (simType=="nonInteractive") {
        a1_meanlog = 0
        a2_meanlog = 0
        a1_othersensitivity = 1
        a2_othersensitivity = 1
      } else if (simType=="a2interactive") {
        a1_meanlog = 0
        a2_meanlog = runif(1,min=meanlog_range[1],max=meanlog_range[2])
        a1_othersensitivity = 1
        a2_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
      } else if (simType=="bidirectional") {
        a1_meanlog = runif(1,min=meanlog_range[1],max=meanlog_range[2])
        a2_meanlog = runif(1,min=meanlog_range[1],max=meanlog_range[2])
        a1_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
        a2_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
      }
      
      a1_respsensitivity = 1
      a2_respsensitivity = 1
      
      a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
      a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)
      
      # run the sim
      voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
      a1_voc_record = voc_records[[1]]
      a2_voc_record = voc_records[[2]]
      
      # compute the measures to be used for comparison to human data
      
      a1_total = sum(a1_voc_record)
      a2_total = sum(a2_voc_record)
      
      sim_turnCount = 0
      leadID = 'none'
      leadT = 0
      for (t in 1:sim_length){
        leadT = leadT+1
        if (leadID=='a1' && a2_voc_record[t]==1) {
          if (leadT <= 5) {
            sim_turnCount = sim_turnCount+1
          }
          leadID = 'a2'
          leadT = 0
        } else if (leadID=='a2' && a1_voc_record[t]==1) {
          if (leadT <= 5) {
            sim_turnCount = sim_turnCount+1
          }
          leadID = 'a1'
          leadT = 0
        } else if (leadID=='none' && a2_voc_record[t]==1) {
          leadID = 'a2'
          leadT = 0
        } else if (leadID=='none' && a1_voc_record[t]==1) {
          leadID = 'a1'
          leadT = 0
        }
      }
      
      a1_ivi_record = get_ivis(a1_voc_record)
      a2_ivi_record = get_ivis(a2_voc_record)
      a1_n_ivi_1 = sum(a1_ivi_record==1)
      a2_n_ivi_1 = sum(a2_ivi_record==1)
      a1_ivi_50 = quantile(a1_ivi_record,probs=.5,names=FALSE)
      a2_ivi_50 = quantile(a2_ivi_record,probs=.5,names=FALSE)
      a1_ivi_max = max(a1_ivi_record)
      a2_ivi_max = max(a2_ivi_record)
      
      # match a1 and a2 with chn or adu based on which pairing provides the better fit
      # as an initial step, can try matching based on vocalization count
      # save the simulation parameters and measures to a data frame
      if ( (simType=="a2interactive") || (chn_total < adu_total && a1_total < a2_total) || (chn_total > adu_total && a2_total > a1_total) ){
        df_row = data.frame(chn_sim_minp = a1_minp,
                            chn_sim_maxp = a1_maxp,
                            chn_sim_meanlog = a1_meanlog,
                            chn_sim_sdlog = a1_sdlog,
                            chn_sim_othersensitivity = a1_othersensitivity,
                            chn_sim_respsensitivity = a1_respsensitivity,
                            chn_sim_total = a1_total,
                            chn_sim_n_ivi_1 = a1_n_ivi_1,
                            chn_sim_ivi_50 = a1_ivi_50,
                            chn_sim_ivi_max = a1_ivi_max,
                            adu_sim_minp = a2_minp,
                            adu_sim_maxp = a2_maxp,
                            adu_sim_meanlog = a2_meanlog,
                            adu_sim_sdlog = a2_sdlog,
                            adu_sim_othersensitivity = a2_othersensitivity,
                            adu_sim_respsensitivity = a2_respsensitivity,
                            adu_sim_total = a2_total,
                            adu_sim_n_ivi_1 = a2_n_ivi_1,
                            adu_sim_ivi_50 = a2_ivi_50,
                            adu_sim_ivi_max = a2_ivi_max,
                            sim_turnCount)
        sims_chn_voc_records[[simNum]] <- a1_voc_record
        sims_adu_voc_records[[simNum]] <- a2_voc_record
        sims_chn_ivi_records[[simNum]] <- a1_ivi_record
        sims_adu_ivi_records[[simNum]] <- a2_ivi_record
      } else {
        df_row = data.frame(chn_sim_minp = a2_minp,
                            chn_sim_maxp = a2_maxp,
                            chn_sim_meanlog = a2_meanlog,
                            chn_sim_sdlog = a2_sdlog,
                            chn_sim_othersensitivity = a2_othersensitivity,
                            chn_sim_respsensitivity = a2_respsensitivity,
                            chn_sim_total = a2_total,
                            chn_sim_n_ivi_1 = a2_n_ivi_1,
                            chn_sim_ivi_50 = a2_ivi_50,
                            chn_sim_ivi_max = a2_ivi_max,
                            adu_sim_minp = a1_minp,
                            adu_sim_maxp = a1_maxp,
                            adu_sim_meanlog = a1_meanlog,
                            adu_sim_sdlog = a1_sdlog,
                            adu_sim_othersensitivity = a1_othersensitivity,
                            adu_sim_respsensitivity = a1_respsensitivity,
                            adu_sim_total = a1_total,
                            adu_sim_n_ivi_1 = a1_n_ivi_1,
                            adu_sim_ivi_50 = a1_ivi_50,
                            adu_sim_ivi_max = a1_ivi_max,
                            sim_turnCount)
        sims_chn_voc_records[[simNum]] <- a2_voc_record
        sims_adu_voc_records[[simNum]] <- a1_voc_record
        sims_chn_ivi_records[[simNum]] <- a2_ivi_record
        sims_adu_ivi_records[[simNum]] <- a1_ivi_record
        
      }
      sims_df = rbind(sims_df,df_row)
    }
    
    #################################################
    # Assess the simulations' fits to the human data
    # (In the future, move this into the for loop,
    # so as to be computed after each simulation)
    #################################################
    
    chn_total_logplus1 = log(chn_total+1)
    sims_df$chn_sim_total_logplus1 = log(sims_df$chn_sim_total+1)
    sims_df$chn_total_diff = (sims_df$chn_sim_total_logplus1-chn_total_logplus1)/chn_total_logplus1
    
    chn_ivi_50_logplus1 = log(chn_ivi_50+1)
    sims_df$chn_sim_ivi_50_logplus1 = log(sims_df$chn_sim_ivi_50+1)
    sims_df$chn_ivi_50_diff = (sims_df$chn_sim_ivi_50_logplus1-chn_ivi_50_logplus1)/chn_ivi_50_logplus1
    
    chn_ivi_max_logplus1 = log(max(chn_ivi_record)+1)
    sims_df$chn_sim_ivi_max_logplus1 = log(sims_df$chn_sim_ivi_max+1)
    sims_df$chn_ivi_max_diff = (sims_df$chn_sim_ivi_max_logplus1-chn_ivi_max_logplus1)/chn_ivi_max_logplus1
    
    chn_n_ivi_1 = sum(chn_ivi_record==1)
    chn_n_ivi_1_logplus1 = log(chn_n_ivi_1+1)
    sims_df$chn_sim_n_ivi_1_logplus1 = log(sims_df$chn_sim_n_ivi_1+1)
    sims_df$chn_n_ivi_1_diff = (sims_df$chn_sim_n_ivi_1_logplus1-chn_n_ivi_1_logplus1)/chn_n_ivi_1_logplus1
    
    adu_total_logplus1 = log(adu_total+1)
    sims_df$adu_sim_total_logplus1 = log(sims_df$adu_sim_total+1)
    sims_df$adu_total_diff = (sims_df$adu_sim_total_logplus1-adu_total_logplus1)/adu_total_logplus1
    
    adu_ivi_50_logplus1 = log(adu_ivi_50+1)
    sims_df$adu_sim_ivi_50_logplus1 = log(sims_df$adu_sim_ivi_50+1)
    sims_df$adu_ivi_50_diff = (sims_df$adu_sim_ivi_50_logplus1-adu_ivi_50_logplus1)/adu_ivi_50_logplus1
    
    adu_ivi_max_logplus1 = log(max(adu_ivi_record)+1)
    sims_df$adu_sim_ivi_max_logplus1 = log(sims_df$adu_sim_ivi_max+1)
    sims_df$adu_ivi_max_diff = (sims_df$adu_sim_ivi_max_logplus1-adu_ivi_max_logplus1)/adu_ivi_max_logplus1
    
    adu_n_ivi_1 = sum(adu_ivi_record==1)
    adu_n_ivi_1_logplus1 = log(adu_n_ivi_1+1)
    sims_df$adu_sim_n_ivi_1_logplus1 = log(sims_df$adu_sim_n_ivi_1+1)
    sims_df$adu_n_ivi_1_diff = (sims_df$adu_sim_n_ivi_1_logplus1-adu_n_ivi_1_logplus1)/adu_n_ivi_1_logplus1
    
    turnCount_logplus1 = log(turnCount+1)
    sims_df$sim_turnCount_logplus1 = log(sims_df$sim_turnCount+1)
    sims_df$turnCount_diff = (sims_df$sim_turnCount_logplus1-turnCount_logplus1)/turnCount_logplus1
    
    # Get fit without considering turn count
    sims_df$simDist_noTurns = sqrt((sims_df$chn_total_diff)^2
                                   +(sims_df$chn_ivi_50_diff)^2
                                   +(sims_df$chn_ivi_max_diff)^2
                                   +(sims_df$chn_n_ivi_1_diff)^2
                                   +(sims_df$adu_total_diff)^2
                                   +(sims_df$adu_ivi_50_diff)^2
                                   +(sims_df$adu_ivi_max_diff)^2
                                   +(sims_df$adu_n_ivi_1_diff)^2)
    
    # Get fit including turn count
    sims_df$simDist_wTurns = sqrt((sims_df$chn_total_diff)^2
                                  +(sims_df$chn_ivi_50_diff)^2
                                  +(sims_df$chn_ivi_max_diff)^2
                                  +(sims_df$chn_n_ivi_1_diff)^2
                                  +(sims_df$adu_total_diff)^2
                                  +(sims_df$adu_ivi_50_diff)^2
                                  +(sims_df$adu_ivi_max_diff)^2
                                  +(sims_df$adu_n_ivi_1_diff)^2
                                  +(sims_df$turnCount_diff)^2)
    
    # Get fit based only on the turn count
    sims_df$simDist_onlyTurns = sqrt((sims_df$turnCount_diff)^2)
    
    fitOrder_noTurns = order(sims_df$simDist_noTurns)
    fitOrder_wTurns = order(sims_df$simDist_wTurns)
    fitOrder_onlyTurns = order(sims_df$simDist_onlyTurns)
    
    ####################################################
    # save the simulation data then reset the workspace
    ####################################################
    
    results_dir = paste(recording_dir,"/",simType,sep="")
    
    if (!dir.exists(results_dir)){
      dir.create(results_dir)
    }
    
    # save the full workspace
    save.image(paste(results_dir,"/VocEventSim_Optimize",".RData",sep=""))
    
    rm(list=ls()[! ls() %in% c("recording","recordings","simType","simTypes","nSims")])
  }
}

