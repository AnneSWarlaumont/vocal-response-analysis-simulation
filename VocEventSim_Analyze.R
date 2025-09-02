setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

recordingsToAnalyze = c("0054_000603","0196_000902","0274_000221","0300_000607","0344_000913","0437_010603","0833_010606")
simTypesToAnalyze = c("nonInteractive","a2interactive","bidirectional")
simsPerRecording = 20

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
# Identify the best-matched simulations whose data should be analyzed.
# We should take the top 20 best-fitting simulations for each human infant
# daylong recording. However, we don't want repeats of any simulation within
# the dataset. So if there are cases of the same simulation being among the 20
# best-fitting for multiple human recordings, we need a protocol for selecting
# which human recording gets to keep its match with that simulation, and for
# finding a replacement simulation for the other human recording(s).
# One way to do this would be to start with the first best simulation.
# Then randomly order the human recordings and then in that random order, get
# the best fit simulation. If the best fit simulation in any case has already
# been used, replace it with the next best fit simulation. Do this recursively
# as needed. And then repeat 20 more times.
################################################################################

matches = data.frame(simType = character(),
                    recording = character(),
                    simNum = integer(),
                    fitRank = integer())

for (simTypeToMatch in simTypesToAnalyze){
  
  simFits_subset = subset(allSimFits,simType==simTypeToMatch)
  used_sims = c()
  nextFitRanks = setNames(rep(1,length(recordingsToAnalyze)),recordingsToAnalyze)
  
  for (i in 1:simsPerRecording){
    
    rand_sort_recordings = sample(recordingsToAnalyze)
    
    for (current_recording in rand_sort_recordings){
      
      simFits_subsubset = subset(simFits_subset,recording==current_recording)
      
      proposed_sim=NA
      while(!any(proposed_sim %in% used_sims)){
        nextFitRank = nextFitRanks[current_recording]
        simFits_subsubsubset = subset(simFits_subsubset,fitRank==nextFitRank)
        proposed_sim = simFits_subsubsubset$simNum
        used_sims = c(used_sims,proposed_sim)
        nextFitRanks[current_recording]<-nextFitRanks[current_recording]+1
      }
      matches_row = data.frame(simType = simTypeToMatch,
                               recording = current_recording,
                               simNum = proposed_sim,
                               fitRank = nextFitRank)
      matches = rbind(matches,matches_row)
      
    }

  }
  
}

################################################################################
# Now for those best-matched simulations, save out their vocalization records,
# and IVI records, and IVI response records to for upcoming statistical analyses.
# The matching human recording can be found in the matches data frame and
# determines what portion of the simulation to include.
################################################################################

source("get_ivis.R")
source("get_ivi_responses.R")
rthresh = 5
warmup <- 1*60*60 # from VocEventSim_RunBatch.R

rec_lengths <- setNames(vector("list",length(recordingsToAnalyze)),recordingsToAnalyze)
for (recording in recordingsToAnalyze){
  recording_dir = paste("data/",recording,sep="")
  lena_segments = read.csv(paste(recording_dir,"/",recording,"_segments.csv",sep=""))
  rec_lengths[recording] = floor(lena_segments$endsec[nrow(lena_segments)])
}

chi_sim_voc_records <- vector("list",nrow(matches))
adu_sim_voc_records <- vector("list",nrow(matches))
chi_sim_ivi_records <- vector("list",nrow(matches))
adu_sim_ivi_records <- vector("list",nrow(matches))
chi_sim_ivi_r_records <- vector("list",nrow(matches))
adu_sim_ivi_r_records <- vector("list",nrow(matches))
currentSimType = "not set"

for (i in 1:nrow(matches)){
  
  recording <- matches$recording[i]
  rec_length <- as.integer(rec_lengths[recording]) # Application of this (plus warmup time still needs to be added in / implemented)
  if (matches$simType[i]!=currentSimType){
    currentSimType <- matches$simType[i]
    simdata_file <- paste("data/",currentSimType,"/VocEventSim_Batch.RData",sep="")
    load(simdata_file)
  }
  
  chi_sim_voc_records[[i]] <- sims_chn_voc_records[[matches$simNum[i]]][(warmup+1):(warmup+rec_length)]
  adu_sim_voc_records[[i]] <- sims_adu_voc_records[[matches$simNum[i]]][(warmup+1):(warmup+rec_length)]
  chi_sim_ivi_records[[i]] <- get_ivis(chi_sim_voc_records[[i]])
  adu_sim_ivi_records[[i]] <- get_ivis(adu_sim_voc_records[[i]])
  
  chi_sim_ivi_r_records[[i]] <- get_ivi_responses(chi_sim_ivi_records[[i]],
                                                  adu_sim_voc_records[[i]],
                                                  rthresh)
  adu_sim_ivi_r_records[[i]] <- get_ivi_responses(adu_sim_ivi_records[[i]],
                                                  chi_sim_voc_records[[i]],
                                                  rthresh)
}

################################################################################
# Now for those best-matched simulations' IVI records, predict the current IVI
# based on the previous 1 or previous 3 IVIs, and save out both sets of
# predictions and their residuals.
################################################################################

chi_sim_previvi_resids <- vector("list",nrow(matches))
adu_sim_previvi_resids <- vector("list",nrow(matches))
chi_sim_prev3ivi_resids <- vector("list",nrow(matches))
adu_sim_prev3ivi_resids <- vector("list",nrow(matches))

for (i in 1:nrow(matches)){
  
  chi_n_ivi = length(chi_sim_ivi_records[[i]])
  adu_n_ivi = length(adu_sim_ivi_records[[i]])
  
  # Correlate current IVI with previous IVI and get the residuals of that
  # correlation
  chi_sim_previvi_model <- lm(scale(log(chi_sim_ivi_records[[i]][2:chi_n_ivi]))
                              ~scale(log(chi_sim_ivi_records[[i]][1:(chi_n_ivi-1)])))
  chi_sim_previvi_resids[[i]] <- resid(chi_sim_previvi_model)[3:(chi_n_ivi-1)]
  adu_sim_previvi_model <- lm(scale(log(adu_sim_ivi_records[[i]][2:adu_n_ivi]))
                              ~scale(log(adu_sim_ivi_records[[i]][1:(adu_n_ivi-1)])))
  adu_sim_previvi_resids[[i]] <- resid(adu_sim_previvi_model)[3:(adu_n_ivi-1)]
  
  # Correlate current IVI with time since the past 3 IVIs and get the residuals
  # of that correlation 
  chi_sim_prev3ivi_model<-lm(scale(log(chi_sim_ivi_records[[i]][4:chi_n_ivi]))
                             ~scale(log(chi_sim_ivi_records[[i]][3:(chi_n_ivi-1)]))
                             +(scale(log(chi_sim_ivi_records[[i]][3:(chi_n_ivi-1)])
                                     +log(chi_sim_ivi_records[[i]][2:(chi_n_ivi-2)])))
                             +(scale(log(chi_sim_ivi_records[[i]][3:(chi_n_ivi-1)])
                                     +log(chi_sim_ivi_records[[i]][2:(chi_n_ivi-2)])
                                     +log(chi_sim_ivi_records[[i]][1:(chi_n_ivi-3)]))))
  chi_sim_prev3ivi_resids[[i]] <- resid(chi_sim_prev3ivi_model)
  adu_sim_prev3ivi_model<-lm(scale(log(adu_sim_ivi_records[[i]][4:adu_n_ivi]))
                             ~scale(log(adu_sim_ivi_records[[i]][3:(adu_n_ivi-1)]))
                             +(scale(log(adu_sim_ivi_records[[i]][3:(adu_n_ivi-1)])
                                     +log(adu_sim_ivi_records[[i]][2:(adu_n_ivi-2)])))
                             +(scale(log(adu_sim_ivi_records[[i]][3:(adu_n_ivi-1)])
                                     +log(adu_sim_ivi_records[[i]][2:(adu_n_ivi-2)])
                                     +log(adu_sim_ivi_records[[i]][1:(adu_n_ivi-3)]))))
  adu_sim_prev3ivi_resids[[i]] <- resid(adu_sim_prev3ivi_model)
  
}

################################################################################
# For the best-matched simulations,
# analyze the chn and adu ivis to test for response effects on ivis using:
# 1. no control for previous ivi,
# 2. control for previous ivi, or
# 3. control for 3 previous ivis.
################################################################################

source('VocEventSim.R') # Need to reload this here, as a different version of
# the function is contained in the workspace that were loaded for each
# simulation batch
source('analyze_ivis.R')

chi_sim_response_results = data.frame(simType = character(),
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

adu_sim_response_results = data.frame(simType = character(),
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
  
  # set up the child simulation data
  ivis = integer()
  ivi_rs = integer()
  previvi_resids = double()
  prev3ivi_resids = double()
  simIDs = character()
  for (i in which(matches$simType==simTypeToA)){
    n <- length(chi_sim_ivi_records[[i]])
    ivis <- c(ivis,chi_sim_ivi_records[[i]][4:n])
    ivi_rs = c(ivi_rs,chi_sim_ivi_r_records[[i]][4:n])
    previvi_resids = c(previvi_resids,chi_sim_previvi_resids[[i]])
    prev3ivi_resids = c(prev3ivi_resids,chi_sim_prev3ivi_resids[[i]])
    simIDs = c(simIDs,rep(as.character(matches$simNum[[i]]),n-3))
  }
  # Analyze child sim data using IVI-based approaches
  # (our original and its variants controlling for previous IVIs)
  chi_sim_ivi_models = analyze_ivis(ivis,ivi_rs,previvi_resids,prev3ivi_resids,simIDs)
  chi_sim_uncontrolled_response_model = chi_sim_ivi_models[[1]]
  chi_sim_residual_response_model = chi_sim_ivi_models[[2]]
  chi_sim_prev3residual_response_model = chi_sim_ivi_models[[3]]
  
  # set up the adult simulation data
  ivis = integer()
  ivi_rs = integer()
  previvi_resids = double()
  prev3ivi_resids = double()
  simIDs = character()
  for (i in which(matches$simType==simTypeToA)){
    n <- length(adu_sim_ivi_records[[i]])
    ivis <- c(ivis,adu_sim_ivi_records[[i]][4:n])
    ivi_rs = c(ivi_rs,adu_sim_ivi_r_records[[i]][4:n])
    previvi_resids = c(previvi_resids,adu_sim_previvi_resids[[i]])
    prev3ivi_resids = c(prev3ivi_resids,adu_sim_prev3ivi_resids[[i]])
    simIDs = c(simIDs,rep(as.character(matches$simNum[[i]]),n-3))
  }
  # Analyze adult sim data using IVI-based approaches
  # (our original and its variants controlling for previous IVIs)
  adu_sim_ivi_models = analyze_ivis(ivis,ivi_rs,previvi_resids,prev3ivi_resids,simIDs)
  adu_sim_uncontrolled_response_model = adu_sim_ivi_models[[1]]
  adu_sim_residual_response_model = adu_sim_ivi_models[[2]]
  adu_sim_prev3residual_response_model = adu_sim_ivi_models[[3]]
  
  chi_sim_rSummary0 = summary(chi_sim_uncontrolled_response_model)
  chi_sim_rBeta0 = chi_sim_rSummary0$coefficients["ivi_response_records","Estimate"]
  chi_sim_rP0 = chi_sim_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_sim_uncontrolled_response_model)=="lm"){
    chi_sim_rCI0 = confint(chi_sim_uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_sim_rCI0 = confint.merMod(chi_sim_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_sim_rSummary0 = summary(adu_sim_uncontrolled_response_model)
  adu_sim_rBeta0 = adu_sim_rSummary0$coefficients["ivi_response_records","Estimate"]
  adu_sim_rP0 = adu_sim_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_sim_uncontrolled_response_model)=="lm"){
    adu_sim_rCI0 = confint(adu_sim_uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_sim_rCI0 = confint.merMod(adu_sim_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  
  # old code to be adapted is below
  
  chi_sim_rSummary1 = summary(chi_sim_residual_response_model)
  chi_sim_rBeta1 = chi_sim_rSummary1$coefficients["ivi_response_records","Estimate"]
  chi_sim_rP1 = chi_sim_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_sim_residual_response_model)=="lm"){
    chi_sim_rCI1 = confint(chi_sim_residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_sim_rCI1 = confint.merMod(chi_sim_residual_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_sim_rSummary1 = summary(adu_sim_residual_response_model)
  adu_sim_rBeta1 = adu_sim_rSummary1$coefficients["ivi_response_records","Estimate"]
  adu_sim_rP1 = adu_sim_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_sim_residual_response_model)=="lm"){
    adu_sim_rCI1 = confint(adu_sim_residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_sim_rCI1 = confint.merMod(adu_sim_residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  chi_sim_rSummary3 = summary(chi_sim_prev3residual_response_model)
  chi_sim_rBeta3 = chi_sim_rSummary3$coefficients["ivi_response_records","Estimate"]
  chi_sim_rP3 = chi_sim_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_sim_prev3residual_response_model)=="lm"){
    chi_sim_rCI3 = confint(chi_sim_prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_sim_rCI3 = confint.merMod(chi_sim_prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_sim_rSummary3 = summary(adu_sim_prev3residual_response_model)
  adu_sim_rBeta3 = adu_sim_rSummary3$coefficients["ivi_response_records","Estimate"]
  adu_sim_rP3 = adu_sim_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_sim_prev3residual_response_model)=="lm"){
    adu_sim_rCI3 = confint(adu_sim_prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_sim_rCI3 = confint.merMod(adu_sim_prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  chi_sim_newrow = data.frame(simType = simTypeToA,
                              chi_sim_rBeta0 = chi_sim_rBeta0,
                              rBeta0Lower = chi_sim_rCI0[1,1],
                              rBeta0Upper = chi_sim_rCI0[1,2],
                              rP0 = chi_sim_rP0,
                              rBeta1 = chi_sim_rBeta1,
                              rBeta1Lower = chi_sim_rCI1[1,1],
                              rBeta1Upper = chi_sim_rCI1[1,2],
                              rP1 = chi_sim_rP1,
                              rBeta3 = chi_sim_rBeta3,
                              rBeta3Lower = chi_sim_rCI3[1,1],
                              rBeta3Upper = chi_sim_rCI3[1,2],
                              rP3 = chi_sim_rP3)
  chi_sim_response_results = rbind(chi_sim_response_results,chi_sim_newrow)
  
  adu_sim_newrow = data.frame(simType = simTypeToA,
                          adu_sim_rBeta0 = adu_sim_rBeta0,
                          rBeta0Lower = adu_sim_rCI0[1,1],
                          rBeta0Upper = adu_sim_rCI0[1,2],
                          rP0 = adu_sim_rP0,
                          rBeta1 = adu_sim_rBeta1,
                          rBeta1Lower = adu_sim_rCI1[1,1],
                          rBeta1Upper = adu_sim_rCI1[1,2],
                          rP1 = adu_sim_rP1,
                          rBeta3 = adu_sim_rBeta3,
                          rBeta3Lower = adu_sim_rCI3[1,1],
                          rBeta3Upper = adu_sim_rCI3[1,2],
                          rP3 = adu_sim_rP3)
  adu_sim_response_results = rbind(adu_sim_response_results,adu_sim_newrow)
  
}

write.csv(chi_sim_response_results, file = "data/chi_sim_response_results.csv")
write.csv(adu_sim_response_results, file = "data/adu_sim_response_results.csv")

################################################################################
# For the best-matched simulations,
# run other tests for adult effects on child from the extant literature
# 1. Su et al.'s lag sequential method (RD index)
# 2. Pretzer et al.'s sequential prediction approach
# 3. Harbison et al.'s child vocal reciprocity test
# 4. de Barbaro et al.'s Granger causality approach
################################################################################
