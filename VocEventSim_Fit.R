setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

recordings = c("0054_000603","0196_000902","0274_000221","0300_000607","0344_000913","0437_010603","0833_010606") # once all the recordings I have queued up have their simulations completed.
simTypes = c("nonInteractive","a2interactive","bidirectional")
warmup <- 1*60*60 # a one-hour simulation warmup period before matching to the human recording begins

simfits_df <- data.frame(simType = character(),
                         recording = character(),
                         simNum = integer(),
                         simDist = double(),
                         fitRank = integer(),
                         stringsAsFactors = FALSE)

# Create a function that calculates IVIs based on a vocal events series
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

# Create a function that counts conversational turns based on a 5 s threshold
# Initialize leader ID and time
# Loop through each second of the recording
# When a turn is identified, update the turn count
getTurnCount <- function(a1_voc_record,a2_voc_record){
  rec_length = length(a1_voc_record)
  turnCount = 0
  leadID = 'none'
  leadT = 0
  for (t in 1:rec_length){
    leadT = leadT+1
    if (leadID=='a1' && a2_voc_record[t]==1) {
      if (leadT <= 5) {
        turnCount = turnCount+1
      }
      leadID = 'a2'
      leadT = 0
    } else if (leadID=='a2' && a1_voc_record[t]==1) {
      if (leadT <= 5) {
        turnCount = turnCount+1
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
  return(turnCount)
}

##############################################################################
# For each simType, assess the fit between each simulation and each recording
##############################################################################

for (simType in simTypes){
  
  # Load in the simulation data
  simdata_file = paste("data/",simType,"/VocEventSim_Batch.RData",sep="")
  load(simdata_file)
  
  nSims = nrow(sims_df)
  
  for (recording in recordings){
    
    simfits_df_temp <- data.frame(simType = character(),
                                  recording = character(),
                                  simNum = integer(),
                                  simDist = double(),
                                  stringsAsFactors = FALSE)
    
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
    
    # get ivi records
    chn_ivi_record = get_ivis(chn_voc_record)
    adu_ivi_record = get_ivis(adu_voc_record)

    # Number of child and adult vocalizations
    chn_total = sum(chn_voc_record)
    adu_total = sum(adu_voc_record)

    # Turn count
    turnCount = getTurnCount(chn_voc_record,adu_voc_record)

    # Three ivi measures (anchors at the distribution's left, middle, and tail)
    chn_n_ivi_1 = sum(chn_ivi_record==1)
    adu_n_ivi_1 = sum(adu_ivi_record==1)
    chn_ivi_50 = quantile(chn_ivi_record,probs=.5,names=FALSE)
    adu_ivi_50 = quantile(adu_ivi_record,probs=.5,names=FALSE)
    chn_ivi_max = max(chn_ivi_record)
    adu_ivi_max = max(adu_ivi_record)
    
    # add 1 then take the log of all the measures
    turnCount_logplus1 = log(turnCount+1)
    chn_total_logplus1 = log(chn_total+1)
    adu_total_logplus1 = log(adu_total+1)
    chn_n_ivi_1_logplus1 = log(chn_n_ivi_1+1)
    adu_n_ivi_1_logplus1 = log(adu_n_ivi_1+1)
    chn_ivi_50_logplus1 = log(chn_ivi_50+1)
    adu_ivi_50_logplus1 = log(adu_ivi_50+1)
    chn_ivi_max_logplus1 = log(chn_ivi_max+1)
    adu_ivi_max_logplus1 = log(adu_ivi_max+1)
    # (these measures will later get put on the same scale that turns the
    # simulations' measures into z-scores across simulations of a given type)
    
    # Get the same measures for each simulation,
    # scale the values,
    # compute scaled differences compared to the reference human recording,
    # then use those differences to compute simDists,
    # then get fit ranks.
    simDists <- c()
    simDists_noTurns <- c()
    simDists_wTurns <- c()
    for (sim in 1:nrow(sims_df)){
      
      sim_chn_voc_record = sims_chn_voc_records[[sim]][(warmup+1):(warmup+length(chn_voc_record))]
      sim_adu_voc_record = sims_adu_voc_records[[sim]][(warmup+1):(warmup+length(adu_voc_record))]
      
      sim_turnCount = getTurnCount(sim_chn_voc_record,sim_adu_voc_record)
      sim_chn_total = sum(sim_chn_voc_record)
      sim_adu_total = sum(sim_adu_voc_record)
      
      sim_chn_ivi_record = get_ivis(sim_chn_voc_record)
      sim_adu_ivi_record = get_ivis(sim_adu_voc_record)
      sim_chn_n_ivi_1 = sum(sim_chn_ivi_record==1)
      sim_adu_n_ivi_1 = sum(sim_adu_ivi_record==1)
      sim_chn_ivi_50 = quantile(sim_chn_ivi_record,probs=.5,names=FALSE)
      sim_adu_ivi_50 = quantile(sim_adu_ivi_record,probs=.5,names=FALSE)
      sim_chn_ivi_max = max(sim_chn_ivi_record)
      sim_adu_ivi_max = max(sim_adu_ivi_record)
      
      # add 1 then take the log of all the measures
      sim_turnCount_logplus1 = log(sim_turnCount+1)
      sim_chn_total_logplus1 = log(sim_chn_total+1)
      sim_adu_total_logplus1 = log(sim_adu_total+1)
      sim_chn_n_ivi_1_logplus1 = log(sim_chn_n_ivi_1+1)
      sim_adu_n_ivi_1_logplus1 = log(sim_adu_n_ivi_1+1)
      sim_chn_ivi_50_logplus1 = log(sim_chn_ivi_50+1)
      sim_adu_ivi_50_logplus1 = log(sim_adu_ivi_50+1)
      sim_chn_ivi_max_logplus1 = log(sim_chn_ivi_max+1)
      sim_adu_ivi_max_logplus1 = log(sim_adu_ivi_max+1)
      
      # to get a scaled difference value,
      # subtract and then divide by the human metric
      turnCount_scaleddiff <- (sim_turnCount_logplus1-turnCount_logplus1)/turnCount_logplus1
      chn_total_scaleddiff <- (sim_chn_total_logplus1-chn_total_logplus1)/chn_total_logplus1
      adu_total_scaleddiff <- (sim_adu_total_logplus1-adu_total_logplus1)/adu_total_logplus1
      chn_n_ivi_1_scaleddiff <- (sim_chn_n_ivi_1_logplus1-chn_n_ivi_1_logplus1)/chn_n_ivi_1_logplus1
      adu_n_ivi_1_scaleddiff <- (sim_adu_n_ivi_1_logplus1-adu_n_ivi_1_logplus1)/adu_n_ivi_1_logplus1
      chn_ivi_50_scaleddiff <- (sim_chn_ivi_50_logplus1-chn_ivi_50_logplus1)/chn_ivi_50_logplus1
      adu_ivi_50_scaleddiff <- (sim_adu_ivi_50_logplus1-adu_ivi_50_logplus1)/adu_ivi_50_logplus1
      chn_ivi_max_scaleddiff <- (sim_chn_ivi_max_logplus1-chn_ivi_max_logplus1)/chn_ivi_max_logplus1
      adu_ivi_max_scaleddiff <- (sim_adu_ivi_max_logplus1-adu_ivi_max_logplus1)/adu_ivi_max_logplus1
      
      # Get simDists for the current recording
      # Get fit without considering turn count
      simDist_noTurns = sqrt((chn_total_scaleddiff)^2
                             +(adu_total_scaleddiff)^2
                             +(chn_n_ivi_1_scaleddiff)^2
                             +(adu_n_ivi_1_scaleddiff)^2
                             +(chn_ivi_50_scaleddiff)^2
                             +(adu_ivi_50_scaleddiff)^2
                             +(chn_ivi_max_scaleddiff)^2
                             +(adu_ivi_max_scaleddiff)^2)
      # Get fit including turn count
      simDist_wTurns = sqrt((chn_total_scaleddiff)^2
                            +(adu_total_scaleddiff)^2
                            +(chn_n_ivi_1_scaleddiff)^2
                            +(adu_n_ivi_1_scaleddiff)^2
                            +(chn_ivi_50_scaleddiff)^2
                            +(adu_ivi_50_scaleddiff)^2
                            +(chn_ivi_max_scaleddiff)^2
                            +(adu_ivi_max_scaleddiff)^2
                            +(turnCount_scaleddiff)^2)
      if (simType=="nonInteractive"){
        simDist = simDist_noTurns
      } else if ((simType=="a2interactive")||(simType=="bidirectional")){
        simDist = simDist_wTurns 
      }
      simDists = c(simDists,simDist)
      simDists_noTurns = c(simDists_noTurns,simDist_noTurns)
      simDists_wTurns = c(simDists_wTurns,simDist_wTurns)
      simfits_df_newrow = data.frame(simType = simType,
                                     recording = recording,
                                     simNum = sim,
                                     simDist = simDist,
                                     stringsAsFactors = FALSE)
      simfits_df_temp = rbind(simfits_df_temp,simfits_df_newrow)
    }
    fitRank = rank(simDists)
    simfits_df_temp = cbind(simfits_df_temp,fitRank)
    simfits_df = rbind(simfits_df,simfits_df_temp)
    write.csv(simfits_df,file="data/simfits.csv")
  }
}
