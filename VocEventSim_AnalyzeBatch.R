setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

recordingsToAnalyze = c("0344_000913","0833_010606","0054_000603")
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
# with and without turn count included
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
# To-do: for the best-matched simulations,
# analyze the chn and adu ivis to test for response effects on ivis using:
# 1. no control for previous ivi,
# 2. control for previous ivi, or
# 3. control for 3 previous ivis
################################################################################

for (simTypeToA in simTypesToAnalyze){
  for (recordingToA in recordingsToAnalyze){
    
    # Load in the simulation data
    recording_file = paste("data/",recordingToA,"/",simTypeToA,"/VocEventSim_Optimize.RData",sep="")
    load(recording_file)
    
    simID = 0
    ivi_records = c()
    ivi_r_records = c()
    simIDs=c()
    previvi_resids = c()
    prev3ivi_resids = c()
    
    for (fitRank in 1:20){

      simID = simID + 1
      #ivi_series = sims_chn_ivi_records[[fitOrder[fitRank]]] # HERE in adaptation of old code below.
    }
    
  }
}

#### Code below is not yet adapted. Will eventually be adapted or deleted.
#
# 
# for (fitRank in 1:20){
#   
#   ivi_series = sims_chn_ivi_records[[fitOrder[fitRank]]]
#   chn_voc_record = sims_chn_voc_records[[fitOrder[fitRank]]]
#   adu_voc_record = sims_adu_voc_records[[fitOrder[fitRank]]]
#   
#   n_ivi = length(ivi_series)
#   
#   
#   ivi_r_record = c() # initialize the record of when adu responded to chn
#   t = 1
#   for (i in 1:n_ivi){
#     if (ivi_series[i]==1){
#       ivi_r_record[i] = NA
#     } else if (sum(adu_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
#       ivi_r_record[i] = 1
#     } else{
#       ivi_r_record[i] = 0
#     }
#     t = t+ivi_series[i]
#   }
#   
#   ivi_records = c(ivi_records,ivi_series[4:n_ivi])
#   ivi_r_records = c(ivi_r_records,ivi_r_record[4:n_ivi])
#   simIDs = c(simIDs,rep(simID,(n_ivi-3)))
#   
#   # Correlate current IVI with previous IVI and get the residuals of that correlation
#   previvi_model = lm(scale(log(ivi_series[2:n_ivi]))~scale(log(ivi_series[1:(n_ivi-1)])))
#   previvi_resid = resid(previvi_model)
#   previvi_resids = c(previvi_resids,previvi_resid[3:(n_ivi-1)])
#   
#   # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
#   prev3ivi_model = lm(scale(log(ivi_series[4:n_ivi]))~scale(log(ivi_series[3:(n_ivi-1)]))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])+log(ivi_series[1:(n_ivi-3)]))))
#   prev3ivi_resid = resid(prev3ivi_model)
#   prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid)
#   
# }
# 
# # Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
# ivi_models = analyze_ivis(ivi_records,ivi_r_records,previvi_resids,prev3ivi_resids,simIDs)
# uncontrolled_response_model = ivi_models[[1]]
# residual_response_model = ivi_models[[2]]
# prev3residual_response_model = ivi_models[[3]]
# 
# summary(uncontrolled_response_model)
# summary(residual_response_model)
# summary(prev3residual_response_model)
# 
# #################################################################################
# # Plot fit to human data as a function of parameter value
# #################################################################################
# 
# library(lattice)
# library(latticeExtra) 
# 
# # showing data points on the same color scale 
# nona_sims_df = subset(sims_df, !is.na(simDist) & !is.infinite(simDist))
# levelplot(simDist ~ chn_sim_minp * chn_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
# levelplot(simDist ~ adu_sim_minp * adu_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
# levelplot(simDist ~ chn_sim_sdlog * chn_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
# levelplot(simDist ~ adu_sim_sdlog * adu_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))