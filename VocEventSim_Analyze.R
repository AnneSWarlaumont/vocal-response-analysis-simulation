setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

recordingsToAnalyze = c("0054_000603","0196_000902","0274_000221","0300_000607","0344_000913","0437_010603","0833_010606") # once all the recordings I have queued up have their simulations completed.
simTypesToAnalyze = c("nonInteractive","a2interactive","bidirectional")

allSimFits = read.csv("data/simfits.csv")

bestFitDists_stats = data.frame(simType = character(),
                                mean = double(),
                                sd = double())
bestFitDists_details = data.frame(simType = character(),
                                  recording = character(),
                                  bestFitDist = double())

for (simTypeToA in simTypesToAnalyze){
  simFits_subset = subset(allSimFits,simType==simTypeToA)
  bestFitDists = double()
  for (recordingToA in recordingsToAnalyze){
    simFits_subsubset = subset(simFits_subset,recording==recordingToA)
    bestFitRow = subset(simFits_subsubset,fitRank==1)
    bestFitDist = bestFitRow$simDist
    details_row = data.frame(simType = simTypeToA,
                             recording = recordingToA,
                             bestFitDist = bestFitDist)
    bestFitDists_details = rbind(bestFitDists_details,details_row)
    bestFitDists = c(bestFitDists,bestFitDist)
  }
  bestFitDists_mean = mean(bestFitDists)
  bestFitDists_sd = sd(bestFitDists)
  stats_row = data.frame(simType = simTypeToA,
                         mean = bestFitDists_mean,
                         sd = bestFitDists_sd)
  bestFitDists_stats = rbind(bestFitDists_stats,stats_row)
}
write.csv(bestFitDists_details, file = "data/bestFitDists_details.csv")
write.csv(bestFitDists_stats, file = "data/bestFitDists_stats.csv")

################################################################################
# For the best-matched simulations,
# analyze the chn and adu ivis to test for response effects on ivis using:
# 1. no control for previous ivi,
# 2. control for previous ivi, or
# 3. control for 3 previous ivis.
################################################################################

chi_response_results = data.frame(simType = character(),
                                  rBeta0 = double(),
                                  rBeta0Lower = double(),
                                  rBeta0Upper = double(),
                                  rP0 = double(),
                                  rBeta1 = double(),
                                  rBeta1Lower = double(),
                                  rBeta1Upper = double(),
                                  rP1 = double(),
                                  rBeta3 = double(),
                                  rBeta3Lower = double(),
                                  rBeta3Upper = double(),
                                  rP3 = double())

adu_response_results = data.frame(simType = character(),
                                  rBeta0 = double(),
                                  rBeta0Lower = double(),
                                  rBeta0Upper = double(),
                                  rP0 = double(),
                                  rBeta1 = double(),
                                  rBeta1Lower = double(),
                                  rBeta1Upper = double(),
                                  rP1 = double(),
                                  rBeta3 = double(),
                                  rBeta3Lower = double(),
                                  rBeta3Upper = double(),
                                  rP3 = double())

##############################################################################
# The code below is old; I am adapting it to the all-in-one batch simulations
# fitOrder has been renamed fitRank
# and there are no longer multiple different fitTypes
# and I'm in the process of removing the assumption that chi is default
# (adding the "chi_" prefix to various variable names)
##############################################################################

for (simTypeToA in simTypesToAnalyze){
  
  ivi_records = c()
  ivi_r_records = c()
  simIDs = c()
  previvi_resids = c()
  prev3ivi_resids = c()
  adu_ivi_records = c()
  adu_ivi_r_records = c()
  adu_simIDs = c()
  adu_previvi_resids = c()
  adu_prev3ivi_resids = c()
  
  simID = 0
  
  for (recordingToA in recordingsToAnalyze){
    
    # Load in the simulation data
    recording_file = paste("data/",recordingToA,"/",simTypeToA,"/VocEventSim_Optimize.RData",sep="")
    load(recording_file)
    
    if (simTypeToA == "nonInteractive"){
      fitOrder = fitOrder_noTurns
    } else{
      fitOrder = fitOrder_wTurns
    }
    
    for (fitRank in 1:20){

      simID = simID + 1
      
      ivi_series = sims_chn_ivi_records[[fitOrder[fitRank]]]
      adu_ivi_series = sims_adu_ivi_records[[fitOrder[fitRank]]]
      chn_voc_record = sims_chn_voc_records[[fitOrder[fitRank]]]
      adu_voc_record = sims_adu_voc_records[[fitOrder[fitRank]]]
      
      n_ivi = length(ivi_series)
      adu_n_ivi = length(adu_ivi_series)
      
      ivi_r_record = c() # initialize the record of when adu responded to chn
      t = 1
      for (i in 1:n_ivi){
        if (ivi_series[i]<=rthresh){
          ivi_r_record[i] = NA
        } else if (sum(adu_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
          ivi_r_record[i] = 1
        } else{
          ivi_r_record[i] = 0
        }
        t = t+ivi_series[i]
      }
      
      adu_ivi_r_record = c() # initialize the record of when chn responded to adu
      t = 1
      for (i in 1:adu_n_ivi){
        if (adu_ivi_series[i]<=rthresh){
          adu_ivi_r_record[i] = NA
        } else if (sum(chn_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
          adu_ivi_r_record[i] = 1
        } else{
          adu_ivi_r_record[i] = 0
        }
        t = t+adu_ivi_series[i]
      }
      
      ivi_records = c(ivi_records,ivi_series[4:n_ivi])
      ivi_r_records = c(ivi_r_records,ivi_r_record[4:n_ivi])
      simIDs = c(simIDs,rep(simID,(n_ivi-3)))
      
      adu_ivi_records = c(adu_ivi_records,adu_ivi_series[4:adu_n_ivi])
      adu_ivi_r_records = c(adu_ivi_r_records,adu_ivi_r_record[4:adu_n_ivi])
      adu_simIDs = c(adu_simIDs,rep(simID,(adu_n_ivi-3)))
      
      # Correlate current IVI with previous IVI and get the residuals of that correlation
      previvi_model = lm(scale(log(ivi_series[2:n_ivi]))~scale(log(ivi_series[1:(n_ivi-1)])))
      previvi_resid = resid(previvi_model)
      previvi_resids = c(previvi_resids,previvi_resid[3:(n_ivi-1)])
      adu_previvi_model = lm(scale(log(adu_ivi_series[2:adu_n_ivi]))~scale(log(adu_ivi_series[1:(adu_n_ivi-1)])))
      adu_previvi_resid = resid(adu_previvi_model)
      adu_previvi_resids = c(adu_previvi_resids,adu_previvi_resid[3:(adu_n_ivi-1)])
      
      # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
      prev3ivi_model = lm(scale(log(ivi_series[4:n_ivi]))~scale(log(ivi_series[3:(n_ivi-1)]))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])+log(ivi_series[1:(n_ivi-3)]))))
      prev3ivi_resid = resid(prev3ivi_model)
      prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid)
      adu_prev3ivi_model = lm(scale(log(adu_ivi_series[4:adu_n_ivi]))~scale(log(adu_ivi_series[3:(adu_n_ivi-1)]))+(scale(log(adu_ivi_series[3:(adu_n_ivi-1)])+log(adu_ivi_series[2:(adu_n_ivi-2)])))+(scale(log(adu_ivi_series[3:(adu_n_ivi-1)])+log(adu_ivi_series[2:(adu_n_ivi-2)])+log(adu_ivi_series[1:(adu_n_ivi-3)]))))
      adu_prev3ivi_resid = resid(adu_prev3ivi_model)
      adu_prev3ivi_resids = c(adu_prev3ivi_resids,adu_prev3ivi_resid)
      
    }
  }
  
  # Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
  source('VocEventSim.R') # Need to reload this here, as a different version of the function is contained in the workspace that were loaded for each simulation batch
  ivi_models = analyze_ivis(ivi_records,ivi_r_records,previvi_resids,prev3ivi_resids,simIDs)
  uncontrolled_response_model = ivi_models[[1]]
  residual_response_model = ivi_models[[2]]
  prev3residual_response_model = ivi_models[[3]]
  
  adu_ivi_models = analyze_ivis(adu_ivi_records,adu_ivi_r_records,adu_previvi_resids,adu_prev3ivi_resids,adu_simIDs)
  adu_uncontrolled_response_model = adu_ivi_models[[1]]
  adu_residual_response_model = adu_ivi_models[[2]]
  adu_prev3residual_response_model = adu_ivi_models[[3]]
  
  rSummary0 = summary(uncontrolled_response_model)
  rBeta0 = rSummary0$coefficients["ivi_response_records","Estimate"]
  rP0 = rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(uncontrolled_response_model)=="lm"){
    rCI0 = confint(uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    rCI0 = confint.merMod(uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_rSummary0 = summary(adu_uncontrolled_response_model)
  adu_rBeta0 = adu_rSummary0$coefficients["ivi_response_records","Estimate"]
  adu_rP0 = adu_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_uncontrolled_response_model)=="lm"){
    adu_rCI0 = confint(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_rCI0 = confint.merMod(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  
  rSummary1 = summary(residual_response_model)
  rBeta1 = rSummary1$coefficients["ivi_response_records","Estimate"]
  rP1 = rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(residual_response_model)=="lm"){
    rCI1 = confint(residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    rCI1 = confint.merMod(residual_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_rSummary1 = summary(adu_residual_response_model)
  adu_rBeta1 = adu_rSummary1$coefficients["ivi_response_records","Estimate"]
  adu_rP1 = adu_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_residual_response_model)=="lm"){
    adu_rCI1 = confint(adu_residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_rCI1 = confint.merMod(adu_residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  rSummary3 = summary(prev3residual_response_model)
  rBeta3 = rSummary3$coefficients["ivi_response_records","Estimate"]
  rP3 = rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(prev3residual_response_model)=="lm"){
    rCI3 = confint(prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    rCI3 = confint.merMod(prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_rSummary3 = summary(adu_prev3residual_response_model)
  adu_rBeta3 = adu_rSummary3$coefficients["ivi_response_records","Estimate"]
  adu_rP3 = adu_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_prev3residual_response_model)=="lm"){
    adu_rCI3 = confint(adu_prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_rCI3 = confint.merMod(adu_prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  newrow = data.frame(simType = simTypeToA,
                      rBeta0 = rBeta0,
                      rBeta0Lower = rCI0[1,1],
                      rBeta0Upper = rCI0[1,2],
                      rP0 = rP0,
                      rBeta1 = rBeta1,
                      rBeta1Lower = rCI1[1,1],
                      rBeta1Upper = rCI1[1,2],
                      rP1 = rP1,
                      rBeta3 = rBeta3,
                      rBeta3Lower = rCI3[1,1],
                      rBeta3Upper = rCI3[1,2],
                      rP3 = rP3)
  response_results = rbind(response_results,newrow)
  
  adu_newrow = data.frame(simType = simTypeToA,
                      adu_rBeta0 = adu_rBeta0,
                      rBeta0Lower = adu_rCI0[1,1],
                      rBeta0Upper = adu_rCI0[1,2],
                      rP0 = adu_rP0,
                      rBeta1 = adu_rBeta1,
                      rBeta1Lower = adu_rCI1[1,1],
                      rBeta1Upper = adu_rCI1[1,2],
                      rP1 = adu_rP1,
                      rBeta3 = adu_rBeta3,
                      rBeta3Lower = adu_rCI3[1,1],
                      rBeta3Upper = adu_rCI3[1,2],
                      rP3 = adu_rP3)
  adu_response_results = rbind(adu_response_results,adu_newrow)
}

write.csv(response_results, file = "data/chn_response_results.csv")
write.csv(adu_response_results, file = "data/adu_response_results.csv")


################################################################################
# For the human recordings,
# analyze the chn and adu ivis to test for response effects on ivis using:
# 1. control for 3 previous ivis
# 2. control for previous ivi, and
# 3. no control for previous ivi.
#
# We may eventually want to add controlled response beta as a criterion for
# simulation optimization.
################################################################################

human_response_results = data.frame(recording = character(),
                                    rBeta0 = double(),
                                    rBeta0Lower = double(),
                                    rBeta0Upper = double(),
                                    rP0 = double(),
                                    rBeta1 = double(),
                                    rBeta1Lower = double(),
                                    rBeta1Upper = double(),
                                    rP1 = double(),
                                    rBeta3 = double(),
                                    rBeta3Lower = double(),
                                    rBeta3Upper = double(),
                                    rP3 = double())

adu_human_response_results = data.frame(recording = character(),
                                    rBeta0 = double(),
                                    rBeta0Lower = double(),
                                    rBeta0Upper = double(),
                                    rP0 = double(),
                                    rBeta1 = double(),
                                    rBeta1Lower = double(),
                                    rBeta1Upper = double(),
                                    rP1 = double(),
                                    rBeta3 = double(),
                                    rBeta3Lower = double(),
                                    rBeta3Upper = double(),
                                    rP3 = double())

chn_cont_human_response_results = data.frame(recording = character(),
                                        rBeta0 = double(),
                                        rBeta0Lower = double(),
                                        rBeta0Upper = double(),
                                        rP0 = double(),
                                        rBeta1 = double(),
                                        rBeta1Lower = double(),
                                        rBeta1Upper = double(),
                                        rP1 = double(),
                                        rBeta3 = double(),
                                        rBeta3Lower = double(),
                                        rBeta3Upper = double(),
                                        rP3 = double())


ivi_records = c()
ivi_r_records = c()
recordingIDs = c()
previvi_resids = c()
prev3ivi_resids = c()

adu_ivi_records = c()
adu_ivi_r_records = c()
adu_recordingIDs = c()
adu_previvi_resids = c()
adu_prev3ivi_resids = c()

chn_cont_ivi_records = c()
chn_cont_ivi_r_records = c()
chn_cont_recordingIDs = c()
chn_cont_previvi_resids = c()
chn_cont_prev3ivi_resids = c()

rthresh = 5

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

get_ivis_continuous <- function(a1_segments,a2_segments,rthresh){
  ivi_record = c()
  ivi_r_record = c()
  voc1_endt = NA
  voc2_startt = NA
  for (n in 1:(nrow(a1_segments)-1)){
    voc1_endt = a1_segments[n,]$endsec
    voc2_startt = a1_segments[n+1,]$startsec
    ivi = voc2_startt - voc1_endt
    if (ivi!=0){
      if (ivi<=rthresh){
        ivi_r = NA 
      } else if (any(between(a2_segments$startsec,voc1_endt,voc1_endt+rthresh))){
        ivi_r = 1
      } else{
        ivi_r = 0
      }
      ivi_record = c(ivi_record,ivi)
      ivi_r_record = c(ivi_r_record,ivi_r)
    }
  }
  return(list(ivi_record,ivi_r_record))
}

for (recordingToA in recordingsToAnalyze){
  
  # Read in and format the human sample
  recording_dir = paste("data/",recordingToA,sep="")
  lena_segments = read.csv(paste(recording_dir,"/",recordingToA,"_segments.csv",sep=""))
  chn_segments = subset(lena_segments,segtype=="CHNSP")
  adu_segments = subset(lena_segments,segtype=="FAN"|segtype=="MAN")
  rec_length = floor(lena_segments$endsec[nrow(lena_segments)])
  chn_voc_record = integer(rec_length)
  chn_voc_record[as.integer(floor(chn_segments$startsec))] = 1
  adu_voc_record = integer(rec_length)
  adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1
  
  chn_ivi_record = get_ivis(chn_voc_record)
  n_ivi = length(chn_ivi_record)
  adu_ivi_record = get_ivis(adu_voc_record)
  adu_n_ivi = length(adu_ivi_record)
  
  ivi_r_record = c() # initialize the record of when adu responded to chn
  t = 1
  for (i in 1:n_ivi){
    if (chn_ivi_record[i]<=rthresh){
      ivi_r_record[i] = NA
    } else if (sum(adu_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
      ivi_r_record[i] = 1
    } else{
      ivi_r_record[i] = 0
    }
    t = t+chn_ivi_record[i]
  }
  adu_ivi_r_record = c() # initialize the record of when adu responded to chn
  t = 1
  for (i in 1:adu_n_ivi){
    if (adu_ivi_record[i]<=rthresh){
      adu_ivi_r_record[i] = NA
    } else if (sum(chn_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
      adu_ivi_r_record[i] = 1
    } else{
      adu_ivi_r_record[i] = 0
    }
    t = t+adu_ivi_record[i]
  }
  
  ivi_records = c(ivi_records,chn_ivi_record[4:n_ivi])
  ivi_r_records = c(ivi_r_records,ivi_r_record[4:n_ivi])
  recordingIDs = c(recordingIDs,rep(recordingToA,(n_ivi-3)))
  
  adu_ivi_records = c(adu_ivi_records,adu_ivi_record[4:adu_n_ivi])
  adu_ivi_r_records = c(adu_ivi_r_records,adu_ivi_r_record[4:adu_n_ivi])
  adu_recordingIDs = c(adu_recordingIDs,rep(recordingToA,(adu_n_ivi-3)))
  
  chn_continuous = get_ivis_continuous(chn_segments,adu_segments,rthresh)
  chn_cont_ivi_record = chn_continuous[[1]]
  chn_cont_ivi_r_record = chn_continuous[[2]]
  chn_n_ivi_cont = length(chn_cont_ivi_record)
  chn_cont_ivi_records = c(chn_cont_ivi_records,chn_cont_ivi_record[4:chn_n_ivi_cont])
  chn_cont_ivi_r_records = c(chn_cont_ivi_r_records,chn_cont_ivi_r_record[4:chn_n_ivi_cont])
  chn_cont_recordingIDs = c(chn_cont_recordingIDs,rep(recordingToA,(chn_n_ivi_cont-3)))
  
  adu_continuous = get_ivis_continuous(adu_segments,chn_segments,rthresh)
  adu_ivi_record_continuous = adu_continuous[[1]]
  adu_ivi_r_record_continuous = adu_continuous[[2]]
  
  # Correlate current IVI with previous IVI and get the residuals of that correlation
  previvi_model = lm(scale(log(chn_ivi_record[2:n_ivi]))~scale(log(chn_ivi_record[1:(n_ivi-1)])))
  previvi_resid = resid(previvi_model)
  previvi_resids = c(previvi_resids,previvi_resid[3:(n_ivi-1)])
  adu_previvi_model = lm(scale(log(adu_ivi_record[2:adu_n_ivi]))~scale(log(adu_ivi_record[1:(adu_n_ivi-1)])))
  adu_previvi_resid = resid(adu_previvi_model)
  adu_previvi_resids = c(adu_previvi_resids,adu_previvi_resid[3:(adu_n_ivi-1)])
  chn_cont_previvi_model = lm(scale(log(chn_cont_ivi_record[2:chn_n_ivi_cont]+1))~scale(log(chn_cont_ivi_record[1:(chn_n_ivi_cont-1)]+1)))
  chn_cont_previvi_resid = resid(chn_cont_previvi_model)
  chn_cont_previvi_resids = c(chn_cont_previvi_resids,chn_cont_previvi_resid[3:(chn_n_ivi_cont-1)])
  
  # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
  prev3ivi_model = lm(scale(log(chn_ivi_record[4:n_ivi]))~scale(log(chn_ivi_record[3:(n_ivi-1)]))+(scale(log(chn_ivi_record[3:(n_ivi-1)])+log(chn_ivi_record[2:(n_ivi-2)])))+(scale(log(chn_ivi_record[3:(n_ivi-1)])+log(chn_ivi_record[2:(n_ivi-2)])+log(chn_ivi_record[1:(n_ivi-3)]))))
  prev3ivi_resid = resid(prev3ivi_model)
  prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid)
  adu_prev3ivi_model = lm(scale(log(adu_ivi_record[4:adu_n_ivi]))~scale(log(adu_ivi_record[3:(adu_n_ivi-1)]))+(scale(log(adu_ivi_record[3:(adu_n_ivi-1)])+log(adu_ivi_record[2:(adu_n_ivi-2)])))+(scale(log(adu_ivi_record[3:(adu_n_ivi-1)])+log(adu_ivi_record[2:(adu_n_ivi-2)])+log(adu_ivi_record[1:(adu_n_ivi-3)]))))
  adu_prev3ivi_resid = resid(adu_prev3ivi_model)
  adu_prev3ivi_resids = c(adu_prev3ivi_resids,adu_prev3ivi_resid)
  chn_cont_prev3ivi_model = lm(scale(log(chn_cont_ivi_record[4:chn_n_ivi_cont]+1))~scale(log(chn_cont_ivi_record[3:(chn_n_ivi_cont-1)]+1))+(scale(log(chn_cont_ivi_record[3:(chn_n_ivi_cont-1)]+1)+log(chn_cont_ivi_record[2:(chn_n_ivi_cont-2)]+1)))+(scale(log(chn_cont_ivi_record[3:(chn_n_ivi_cont-1)]+1)+log(chn_cont_ivi_record[2:(chn_n_ivi_cont-2)]+1)+log(chn_cont_ivi_record[1:(chn_n_ivi_cont-3)]+1))))
  chn_cont_prev3ivi_resid = resid(chn_cont_prev3ivi_model)
  chn_cont_prev3ivi_resids = c(chn_cont_prev3ivi_resids,chn_cont_prev3ivi_resid)
  
  # Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
  uncontrolled_response_model = lm(scale(log(chn_ivi_record[4:n_ivi]))~ivi_r_record[4:n_ivi])
  residual_response_model = lm(scale(previvi_resid[3:(n_ivi-1)])~ivi_r_record[4:n_ivi])
  prev3residual_response_model = lm(scale(prev3ivi_resid)~ivi_r_record[4:n_ivi])
  adu_uncontrolled_response_model = lm(scale(log(adu_ivi_record[4:adu_n_ivi]))~adu_ivi_r_record[4:adu_n_ivi])
  adu_residual_response_model = lm(scale(adu_previvi_resid[3:(adu_n_ivi-1)])~adu_ivi_r_record[4:adu_n_ivi])
  adu_prev3residual_response_model = lm(scale(adu_prev3ivi_resid)~adu_ivi_r_record[4:adu_n_ivi])
  chn_cont_uncontrolled_response_model = lm(scale(log(chn_cont_ivi_record[4:chn_n_ivi_cont]+1))~chn_cont_ivi_r_record[4:chn_n_ivi_cont])
  chn_cont_residual_response_model = lm(scale(chn_cont_previvi_resid[3:(chn_n_ivi_cont-1)])~chn_cont_ivi_r_record[4:chn_n_ivi_cont])
  chn_cont_prev3residual_response_model = lm(scale(chn_cont_prev3ivi_resid)~chn_cont_ivi_r_record[4:chn_n_ivi_cont])
  
  rSummary0 = summary(uncontrolled_response_model)
  rBeta0 = rSummary0$coefficients["ivi_r_record[4:n_ivi]","Estimate"]
  rP0 = rSummary0$coefficients["ivi_r_record[4:n_ivi]","Pr(>|t|)"]
  if (class(uncontrolled_response_model)=="lm"){
    rCI0 = confint(uncontrolled_response_model, "ivi_r_record[4:n_ivi]", level = 0.99)
  } else{
    rCI0 = confint.merMod(uncontrolled_response_model, "ivi_r_record[4:n_ivi]", level = 0.99) 
  }
  adu_rSummary0 = summary(adu_uncontrolled_response_model)
  adu_rBeta0 = adu_rSummary0$coefficients["adu_ivi_r_record[4:adu_n_ivi]","Estimate"]
  adu_rP0 = adu_rSummary0$coefficients["adu_ivi_r_record[4:adu_n_ivi]","Pr(>|t|)"]
  if (class(adu_uncontrolled_response_model)=="lm"){
    adu_rCI0 = confint(adu_uncontrolled_response_model, "adu_ivi_r_record[4:adu_n_ivi]", level = 0.99)
  } else{
    adu_rCI0 = confint.merMod(adu_uncontrolled_response_model, "adu_ivi_r_record[4:adu_n_ivi]", level = 0.99) 
  }
  chn_cont_rSummary0 = summary(chn_cont_uncontrolled_response_model)
  chn_cont_rBeta0 = chn_cont_rSummary0$coefficients["chn_cont_ivi_r_record[4:chn_n_ivi_cont]","Estimate"]
  chn_cont_rP0 = chn_cont_rSummary0$coefficients["chn_cont_ivi_r_record[4:chn_n_ivi_cont]","Pr(>|t|)"]
  if (class(chn_cont_uncontrolled_response_model)=="lm"){
    chn_cont_rCI0 = confint(chn_cont_uncontrolled_response_model, "chn_cont_ivi_r_record[4:chn_n_ivi_cont]", level = 0.99)
  } else{
    chn_cont_rCI0 = confint.merMod(chn_cont_uncontrolled_response_model, "chn_cont_ivi_r_record[4:chn_n_ivi_cont]", level = 0.99) 
  }
  
  rSummary1 = summary(residual_response_model)
  rBeta1 = rSummary1$coefficients["ivi_r_record[4:n_ivi]","Estimate"]
  rP1 = rSummary1$coefficients["ivi_r_record[4:n_ivi]","Pr(>|t|)"]
  if (class(residual_response_model)=="lm"){
    rCI1 = confint(residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99)
  } else{
    rCI1 = confint.merMod(residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99) 
  }
  adu_rSummary1 = summary(adu_residual_response_model)
  adu_rBeta1 = adu_rSummary1$coefficients["adu_ivi_r_record[4:adu_n_ivi]","Estimate"]
  adu_rP1 = adu_rSummary1$coefficients["adu_ivi_r_record[4:adu_n_ivi]","Pr(>|t|)"]
  if (class(adu_residual_response_model)=="lm"){
    adu_rCI1 = confint(adu_residual_response_model, "adu_ivi_r_record[4:adu_n_ivi]", level = 0.99)
  } else{
    adu_rCI1 = confint.merMod(adu_residual_response_model, "adu_ivi_r_record[4:adu_n_ivi]", level = 0.99) 
  }
  chn_cont_rSummary1 = summary(chn_cont_residual_response_model)
  chn_cont_rBeta1 = chn_cont_rSummary1$coefficients["chn_cont_ivi_r_record[4:chn_n_ivi_cont]","Estimate"]
  chn_cont_rP1 = chn_cont_rSummary1$coefficients["chn_cont_ivi_r_record[4:chn_n_ivi_cont]","Pr(>|t|)"]
  if (class(chn_cont_residual_response_model)=="lm"){
    chn_cont_rCI1 = confint(chn_cont_residual_response_model, "chn_cont_ivi_r_record[4:chn_n_ivi_cont]", level = 0.99)
  } else{
    chn_cont_rCI1 = confint.merMod(chn_cont_residual_response_model, "chn_cont_ivi_r_record[4:chn_n_ivi_cont]", level = 0.99) 
  }
  
  rSummary3 = summary(prev3residual_response_model)
  rBeta3 = rSummary3$coefficients["ivi_r_record[4:n_ivi]","Estimate"]
  rP3 = rSummary3$coefficients["ivi_r_record[4:n_ivi]","Pr(>|t|)"]
  if (class(prev3residual_response_model)=="lm"){
    rCI3 = confint(prev3residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99)
  } else{
    rCI3 = confint.merMod(prev3residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99) 
  }
  adu_rSummary3 = summary(adu_prev3residual_response_model)
  adu_rBeta3 = adu_rSummary3$coefficients["adu_ivi_r_record[4:adu_n_ivi]","Estimate"]
  adu_rP3 = adu_rSummary3$coefficients["adu_ivi_r_record[4:adu_n_ivi]","Pr(>|t|)"]
  if (class(adu_prev3residual_response_model)=="lm"){
    adu_rCI3 = confint(adu_prev3residual_response_model, "adu_ivi_r_record[4:adu_n_ivi]", level = 0.99)
  } else{
    adu_rCI3 = confint.merMod(adu_prev3residual_response_model, "adu_ivi_r_record[4:adu_n_ivi]", level = 0.99) 
  }
  chn_cont_rSummary3 = summary(chn_cont_prev3residual_response_model)
  chn_cont_rBeta3 = chn_cont_rSummary3$coefficients["chn_cont_ivi_r_record[4:chn_n_ivi_cont]","Estimate"]
  chn_cont_rP3 = chn_cont_rSummary3$coefficients["chn_cont_ivi_r_record[4:chn_n_ivi_cont]","Pr(>|t|)"]
  if (class(chn_cont_prev3residual_response_model)=="lm"){
    chn_cont_rCI3 = confint(chn_cont_prev3residual_response_model, "chn_cont_ivi_r_record[4:chn_n_ivi_cont]", level = 0.99)
  } else{
    chn_cont_rCI3 = confint.merMod(chn_cont_prev3residual_response_model, "chn_cont_ivi_r_record[4:chn_n_ivi_cont]", level = 0.99) 
  }
  
  newrow = data.frame(recording = recordingToA,
                      rBeta0 = rBeta0,
                      rBeta0Lower = rCI0[1,1],
                      rBeta0Upper = rCI0[1,2],
                      rP0 = rP0,
                      rBeta1 = rBeta1,
                      rBeta1Lower = rCI1[1,1],
                      rBeta1Upper = rCI1[1,2],
                      rP1 = rP1,
                      rBeta3 = rBeta3,
                      rBeta3Lower = rCI3[1,1],
                      rBeta3Upper = rCI3[1,2],
                      rP3 = rP3)
  human_response_results = rbind(human_response_results,newrow)
  
  adu_newrow = data.frame(recording = recordingToA,
                      rBeta0 = adu_rBeta0,
                      rBeta0Lower = adu_rCI0[1,1],
                      rBeta0Upper = adu_rCI0[1,2],
                      rP0 = adu_rP0,
                      rBeta1 = adu_rBeta1,
                      rBeta1Lower = adu_rCI1[1,1],
                      rBeta1Upper = adu_rCI1[1,2],
                      rP1 = adu_rP1,
                      rBeta3 = adu_rBeta3,
                      rBeta3Lower = adu_rCI3[1,1],
                      rBeta3Upper = adu_rCI3[1,2],
                      rP3 = adu_rP3)
  adu_human_response_results = rbind(adu_human_response_results,adu_newrow)
  
  chn_cont_newrow = data.frame(recording = recordingToA,
                          rBeta0 = chn_cont_rBeta0,
                          rBeta0Lower = chn_cont_rCI0[1,1],
                          rBeta0Upper = chn_cont_rCI0[1,2],
                          rP0 = chn_cont_rP0,
                          rBeta1 = chn_cont_rBeta1,
                          rBeta1Lower = chn_cont_rCI1[1,1],
                          rBeta1Upper = chn_cont_rCI1[1,2],
                          rP1 = chn_cont_rP1,
                          rBeta3 = chn_cont_rBeta3,
                          rBeta3Lower = chn_cont_rCI3[1,1],
                          rBeta3Upper = chn_cont_rCI3[1,2],
                          rP3 = chn_cont_rP3)
  chn_cont_human_response_results = rbind(chn_cont_human_response_results,chn_cont_newrow)
  
}

# Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
source('VocEventSim.R') # Need to reload this here, as a different version of the function is contained in the workspace that were loaded for each simulation batch
ivi_models = analyze_ivis(ivi_records,ivi_r_records,previvi_resids,prev3ivi_resids,as.factor(recordingIDs))
uncontrolled_response_model = ivi_models[[1]]
residual_response_model = ivi_models[[2]]
prev3residual_response_model = ivi_models[[3]]
adu_ivi_models = analyze_ivis(adu_ivi_records,adu_ivi_r_records,adu_previvi_resids,adu_prev3ivi_resids,as.factor(adu_recordingIDs))
adu_uncontrolled_response_model = adu_ivi_models[[1]]
adu_residual_response_model = adu_ivi_models[[2]]
adu_prev3residual_response_model = adu_ivi_models[[3]]

rSummary0 = summary(uncontrolled_response_model)
rBeta0 = rSummary0$coefficients["ivi_response_records","Estimate"]
rP0 = rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(uncontrolled_response_model)=="lm"){
  rCI0 = confint(uncontrolled_response_model, "ivi_response_records", level = 0.99)
} else{
  rCI0 = confint.merMod(uncontrolled_response_model, "ivi_response_records", level = 0.99) 
}
adu_rSummary0 = summary(adu_uncontrolled_response_model)
adu_rBeta0 = adu_rSummary0$coefficients["ivi_response_records","Estimate"]
adu_rP0 = adu_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(adu_uncontrolled_response_model)=="lm"){
  adu_rCI0 = confint(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99)
} else{
  adu_rCI0 = confint.merMod(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
}

rSummary1 = summary(residual_response_model)
rBeta1 = rSummary1$coefficients["ivi_response_records","Estimate"]
rP1 = rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(residual_response_model)=="lm"){
  rCI0 = confint(residual_response_model, "ivi_response_records", level = 0.99)
} else{
  rCI0 = confint.merMod(residual_response_model, "ivi_response_records", level = 0.99) 
}
adu_rSummary1 = summary(adu_residual_response_model)
adu_rBeta1 = adu_rSummary1$coefficients["ivi_response_records","Estimate"]
adu_rP1 = adu_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(adu_residual_response_model)=="lm"){
  adu_rCI0 = confint(adu_residual_response_model, "ivi_response_records", level = 0.99)
} else{
  adu_rCI0 = confint.merMod(adu_residual_response_model, "ivi_response_records", level = 0.99) 
}

rSummary3 = summary(prev3residual_response_model)
rBeta3 = rSummary3$coefficients["ivi_response_records","Estimate"]
rP3 = rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(prev3residual_response_model)=="lm"){
  rCI0 = confint(prev3residual_response_model, "ivi_response_records", level = 0.99)
} else{
  rCI0 = confint.merMod(prev3residual_response_model, "ivi_response_records", level = 0.99) 
}
adu_rSummary3 = summary(adu_prev3residual_response_model)
adu_rBeta3 = adu_rSummary3$coefficients["ivi_response_records","Estimate"]
adu_rP3 = adu_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(adu_prev3residual_response_model)=="lm"){
  adu_rCI0 = confint(adu_prev3residual_response_model, "ivi_response_records", level = 0.99)
} else{
  adu_rCI0 = confint.merMod(adu_prev3residual_response_model, "ivi_response_records", level = 0.99) 
}

newrow = data.frame(recording = "all",
                    rBeta0 = rBeta0,
                    rBeta0Lower = rCI0[1,1],
                    rBeta0Upper = rCI0[1,2],
                    rP0 = rP0,
                    rBeta1 = rBeta1,
                    rBeta1Lower = rCI1[1,1],
                    rBeta1Upper = rCI1[1,2],
                    rP1 = rP1,
                    rBeta3 = rBeta3,
                    rBeta3Lower = rCI3[1,1],
                    rBeta3Upper = rCI3[1,2],
                    rP3 = rP3)
human_response_results = rbind(human_response_results,newrow)

adu_newrow = data.frame(recording = "all",
                    rBeta0 = adu_rBeta0,
                    rBeta0Lower = adu_rCI0[1,1],
                    rBeta0Upper = adu_rCI0[1,2],
                    rP0 = adu_rP0,
                    rBeta1 = adu_rBeta1,
                    rBeta1Lower = adu_rCI1[1,1],
                    rBeta1Upper = adu_rCI1[1,2],
                    rP1 = adu_rP1,
                    rBeta3 = adu_rBeta3,
                    rBeta3Lower = adu_rCI3[1,1],
                    rBeta3Upper = adu_rCI3[1,2],
                    rP3 = adu_rP3)
adu_human_response_results = rbind(adu_human_response_results,adu_newrow)

write.csv(human_response_results, file = "data/chn_human_response_results.csv")
write.csv(adu_human_response_results, file = "data/adu_human_response_results.csv")

################################################################################
# For the best-matched simulations,
# run other tests for adult effects on child from the extant literature
# 1. Su et al.'s lag sequential method (RD index)
# 2. Pretzer et al.'s sequential prediction approach
# 3. Harbison et al.'s child vocal reciprocity test
# 4. de Barbaro et al.'s Granger causality approach
################################################################################

################################################################################
# TODO: Add code to get response beta for the three analysis methods for the top
#       fit for each simType for each recording. Then find which simType has
#       the response beta set that best matches the response beta set for the
#       human recording. That may be a nice indicator of whether the real human
#       data are best explained as non-interactive, one agent interactive, or
#       both agents interactive.
# TODO: Update all the above to also analyze adult ivis!
################################################################################
