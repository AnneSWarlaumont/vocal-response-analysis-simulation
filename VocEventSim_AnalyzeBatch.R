setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

recordingsToAnalyze = c("0054_000603","0344_000913","0833_010606")
# recordingsToAnalyze = c("0054_000603","0196_000902","0344_000913","0833_010606")
# recordingsToAnalyze = c("0054_000603","0196_000902","0274_000221","0300_000607","0344_000913","0437_010603","0833_010606") # once all the recordings I have queued up have their simulations completed.
simTypesToAnalyze = c("nonInteractive","a2interactive","bidirectional")

allSimFits <- data.frame(recording = character(),
                      simType = character(),
                      simNum = integer(),
                      fitType = character(),
                      simDist = double(),
                      fitOrder = integer(),
                      stringsAsFactors = FALSE)

#############################################################
# Assess the fit between each simulation and each recording,
# with and without turn count included.
#
# It would be good to add tests for statistical significance
# between simulation types.
#
# It may also be good to add an evaluation of fit based on
# match of the 3 prev ivi controlled ivi response beta.
#############################################################

for (recordingToA in recordingsToAnalyze){
  for (simTypeToA in simTypesToAnalyze){
    
    # Load in the simulation data
    recording_file = paste("data/",recordingToA,"/",simTypeToA,"/VocEventSim_Optimize.RData",sep="")
    load(recording_file)
    
    simFits = data.frame(recording = rep(recordingToA,nSims),
                         simType = rep(simTypeToA,nSims),
                         simNum = seq(1,nSims),
                         fitType = rep("noTurns",nSims),
                         simDist = sims_df$simDist_noTurns,
                         fitOrder = fitOrder_noTurns)
    allSimFits = rbind(allSimFits,simFits)
    
    simFits = data.frame(recording = rep(recordingToA,nSims),
                         simType = rep(simTypeToA,nSims),
                         simNum = seq(1,nSims),
                         fitType = rep("wTurns",nSims),
                         simDist = sims_df$simDist_wTurns,
                         fitOrder = fitOrder_wTurns)
    allSimFits = rbind(allSimFits,simFits)
    
    simFits = data.frame(recording = rep(recordingToA,nSims),
                         simType = rep(simTypeToA,nSims),
                         simNum = seq(1,nSims),
                         fitType = rep("onlyTurns",nSims),
                         simDist = sims_df$simDist_onlyTurns,
                         fitOrder = fitOrder_onlyTurns)
    allSimFits = rbind(allSimFits,simFits)
    
    
  }
}

# export allSimFits to csv
write.csv(allSimFits, file = "data/allSimFits.csv")

bestFitDists_stats = data.frame(simType = character(),
                                fitType = character(),
                                mean = double(),
                                sd = double())
bestFitDists_details = data.frame(simType = character(),
                                  fitType = character(),
                                  recording = character(),
                                  bestFitDist = double())
for (simTypeToA in simTypesToAnalyze){
  for (fitTypeToA in c("noTurns","wTurns","onlyTurns")){
    simFits_subset = subset(allSimFits,((simType==simTypeToA)&(fitType==fitTypeToA)))
    bestFitDists = double()
    for (recordingToA in recordingsToAnalyze){
      simFits_subsubset = subset(simFits_subset,recording==recordingToA)
      bestFitRow = subset(simFits_subsubset,fitOrder==1)
      bestFitDist = bestFitRow$simDist
      details_row = data.frame(simType = simTypeToA,
                                        fitType = fitTypeToA,
                                        recording = recordingToA,
                                        bestFitDist = bestFitDist)
      bestFitDists_details = rbind(bestFitDists_details,details_row)
      bestFitDists = c(bestFitDists,bestFitDist)
    }
    bestFitDists_mean = mean(bestFitDists)
    bestFitDists_sd = sd(bestFitDists)
    stats_row = data.frame(simType = simTypeToA,
                                    fitType = fitTypeToA,
                                    mean = bestFitDists_mean,
                                    sd = bestFitDists_sd)
    bestFitDists_stats = rbind(bestFitDists_stats,stats_row)
  }
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

response_results = data.frame(simType = character(),
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

for (simTypeToA in simTypesToAnalyze){
  
  ivi_records = c()
  ivi_r_records = c()
  simIDs = c()
  previvi_resids = c()
  prev3ivi_resids = c()
  
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
      chn_voc_record = sims_chn_voc_records[[fitOrder[fitRank]]]
      adu_voc_record = sims_adu_voc_records[[fitOrder[fitRank]]]
      
      n_ivi = length(ivi_series)
      
      ivi_r_record = c() # initialize the record of when adu responded to chn
      t = 1
      for (i in 1:n_ivi){
        if (ivi_series[i]==1){
          ivi_r_record[i] = NA
        } else if (sum(adu_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
          ivi_r_record[i] = 1
        } else{
          ivi_r_record[i] = 0
        }
        t = t+ivi_series[i]
      }
      
      ivi_records = c(ivi_records,ivi_series[4:n_ivi])
      ivi_r_records = c(ivi_r_records,ivi_r_record[4:n_ivi])
      simIDs = c(simIDs,rep(simID,(n_ivi-3)))
      
      # Correlate current IVI with previous IVI and get the residuals of that correlation
      previvi_model = lm(scale(log(ivi_series[2:n_ivi]))~scale(log(ivi_series[1:(n_ivi-1)])))
      previvi_resid = resid(previvi_model)
      previvi_resids = c(previvi_resids,previvi_resid[3:(n_ivi-1)])
      
      # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
      prev3ivi_model = lm(scale(log(ivi_series[4:n_ivi]))~scale(log(ivi_series[3:(n_ivi-1)]))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])+log(ivi_series[1:(n_ivi-3)]))))
      prev3ivi_resid = resid(prev3ivi_model)
      prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid)
      
    }
  }
  
  # Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
  source('VocEventSim.R') # Need to reload this here, as a different version of the function is contained in the workspace that were loaded for each simulation batch
  ivi_models = analyze_ivis(ivi_records,ivi_r_records,previvi_resids,prev3ivi_resids,simIDs)
  uncontrolled_response_model = ivi_models[[1]]
  residual_response_model = ivi_models[[2]]
  prev3residual_response_model = ivi_models[[3]]
  
  rSummary0 = summary(uncontrolled_response_model)
  rBeta0 = rSummary0$coefficients["ivi_response_records","Estimate"]
  rP0 = rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(uncontrolled_response_model)=="lm"){
    rCI0 = confint(uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    rCI0 = confint.merMod(uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  
  rSummary1 = summary(residual_response_model)
  rBeta1 = rSummary1$coefficients["ivi_response_records","Estimate"]
  rP1 = rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(residual_response_model)=="lm"){
    rCI1 = confint(residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    rCI1 = confint.merMod(residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  rSummary3 = summary(prev3residual_response_model)
  rBeta3 = rSummary3$coefficients["ivi_response_records","Estimate"]
  rP3 = rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(prev3residual_response_model)=="lm"){
    rCI3 = confint(prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    rCI3 = confint.merMod(prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  # print(simTypeToA)
  # print(rSummary0)
  # print(rSummary1)
  # print(rSummary3)
  
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
  
  
  
}

write.csv(response_results, file = "data/response_results.csv")


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

ivi_records = c()
ivi_r_records = c()
recordingIDs = c()
previvi_resids = c()
prev3ivi_resids = c()

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
  n_ivi = length(ivi_series)
  
  print(n_ivi)
  
  ivi_r_record = c() # initialize the record of when adu responded to chn
  t = 1
  for (i in 1:n_ivi){
    if (ivi_series[i]==1){
      ivi_r_record[i] = NA
    } else if (sum(adu_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
      ivi_r_record[i] = 1
    } else{
      ivi_r_record[i] = 0
    }
    t = t+ivi_series[i]
  }
  
  ivi_records = c(ivi_records,ivi_series[4:n_ivi])
  ivi_r_records = c(ivi_r_records,ivi_r_record[4:n_ivi])
  recordingIDs = c(recordingIDs,rep(recordingToA,(n_ivi-3)))
  
  # Correlate current IVI with previous IVI and get the residuals of that correlation
  previvi_model = lm(scale(log(ivi_series[2:n_ivi]))~scale(log(ivi_series[1:(n_ivi-1)])))
  previvi_resid = resid(previvi_model)
  previvi_resids = c(previvi_resids,previvi_resid[3:(n_ivi-1)])
  
  # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
  prev3ivi_model = lm(scale(log(ivi_series[4:n_ivi]))~scale(log(ivi_series[3:(n_ivi-1)]))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])+log(ivi_series[1:(n_ivi-3)]))))
  prev3ivi_resid = resid(prev3ivi_model)
  prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid)
  
  # Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
  uncontrolled_response_model = lm(scale(log(chn_ivi_record[4:n_ivi]))~ivi_r_record[4:n_ivi])
  residual_response_model = lm(scale(previvi_resid[3:(n_ivi-1)])~ivi_r_record[4:n_ivi])
  prev3residual_response_model = lm(scale(prev3ivi_resid)~ivi_r_record[4:n_ivi])
  
  rSummary0 = summary(uncontrolled_response_model)
  rBeta0 = rSummary0$coefficients["ivi_r_record[4:n_ivi]","Estimate"]
  rP0 = rSummary0$coefficients["ivi_r_record[4:n_ivi]","Pr(>|t|)"]
  if (class(uncontrolled_response_model)=="lm"){
    rCI0 = confint(uncontrolled_response_model, "ivi_r_record[4:n_ivi]", level = 0.99)
  } else{
    rCI0 = confint.merMod(uncontrolled_response_model, "ivi_r_record[4:n_ivi]", level = 0.99) 
  }
  
  rSummary1 = summary(residual_response_model)
  rBeta1 = rSummary1$coefficients["ivi_r_record[4:n_ivi]","Estimate"]
  rP1 = rSummary1$coefficients["ivi_r_record[4:n_ivi]","Pr(>|t|)"]
  if (class(residual_response_model)=="lm"){
    rCI1 = confint(residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99)
  } else{
    rCI1 = confint.merMod(residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99) 
  }
  
  rSummary3 = summary(prev3residual_response_model)
  rBeta3 = rSummary3$coefficients["ivi_r_record[4:n_ivi]","Estimate"]
  rP3 = rSummary3$coefficients["ivi_r_record[4:n_ivi]","Pr(>|t|)"]
  if (class(prev3residual_response_model)=="lm"){
    rCI3 = confint(prev3residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99)
  } else{
    rCI3 = confint.merMod(prev3residual_response_model, "ivi_r_record[4:n_ivi]", level = 0.99) 
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
  
}

# Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
source('VocEventSim.R') # Need to reload this here, as a different version of the function is contained in the workspace that were loaded for each simulation batch
ivi_models = analyze_ivis(ivi_records,ivi_r_records,previvi_resids,prev3ivi_resids,recordingIDs)
uncontrolled_response_model = ivi_models[[1]]
residual_response_model = ivi_models[[2]]
prev3residual_response_model = ivi_models[[3]]

rSummary0 = summary(uncontrolled_response_model)
rBeta0 = rSummary0$coefficients["ivi_response_records","Estimate"]
rP0 = rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(uncontrolled_response_model)=="lm"){
  rCI0 = confint(uncontrolled_response_model, "ivi_response_records", level = 0.99)
} else{
  rCI0 = confint.merMod(uncontrolled_response_model, "ivi_response_records", level = 0.99) 
}

rSummary1 = summary(residual_response_model)
rBeta1 = rSummary1$coefficients["ivi_response_records","Estimate"]
rP1 = rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(residual_response_model)=="lm"){
  rCI0 = confint(residual_response_model, "ivi_response_records", level = 0.99)
} else{
  rCI0 = confint.merMod(residual_response_model, "ivi_response_records", level = 0.99) 
}

rSummary3 = summary(prev3residual_response_model)
rBeta3 = rSummary3$coefficients["ivi_response_records","Estimate"]
rP3 = rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(prev3residual_response_model)=="lm"){
  rCI0 = confint(prev3residual_response_model, "ivi_response_records", level = 0.99)
} else{
  rCI0 = confint.merMod(prev3residual_response_model, "ivi_response_records", level = 0.99) 
}

# print(simTypeToA)
# print(rSummary0)
# print(rSummary1)
# print(rSummary3)

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

write.csv(human_response_results, file = "data/human_response_results.csv")

################################################################################
# TODO: Add code to get response beta for the three analysis methods for the top
#       fit for each simType for each recording. Then find which simType has
#       the response beta set that best matches the response beta set for the
#       human recording. That may be a nice indicator of whether the real human
#       data are best explained as non-interactive, one agent interactive, or
#       both agents interactive.
# TODO: Update all the above to also analyze adult ivis!
################################################################################
