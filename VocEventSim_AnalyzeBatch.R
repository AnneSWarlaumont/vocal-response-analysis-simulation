setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

# Will need to revise the data loading since we're no longer going to be merging across simulations
# Load in data:
sims_df = read.csv("~/Documents/GitHub/vocal-response-analysis-simulation/sims_df_merged_20250521.csv")

recordings = c("0344_000913","0833_010606","0054_000603")
simTypes = c("nonInteractive","a2interactive","bidirectional")

#############################################################
# Assess the fit between each simulation and each recording,
# with and without turn count included
#############################################################

# initialize two matrices, which will store the fits of each simulation (in rows) to each recording (in columns)
simFits_noTurns = matrix(data = NA, nrow = nrow(sims_df), ncol = length(recordings))
simFits_wTurns = matrix(data = NA, nrow = nrow(sims_df), ncol = length(recordings))
simFits_onlyTurns = matrix(data = NA, nrow = nrow(sims_df), ncol = length(recordings))
fitOrder_noTurns = matrix(data = NA, nrow = nrow(sims_df), ncol = length(recordings))
fitOrder_wTurns = matrix(data = NA, nrow = nrow(sims_df), ncol = length(recordings))
fitOrder_onlyTurns = matrix(data = NA, nrow = nrow(sims_df), ncol = length(recordings))
colnames(simFits_noTurns) = recordings
colnames(simFits_wTurns) = recordings
colnames(simFits_onlyTurns) = recordings
colnames(fitOrder_noTurns) = recordings
colnames(fitOrder_wTurns) = recordings
colnames(fitOrder_onlyTurns) = recordings

for (recording in recordings){
  recording_dir = paste("data/",recording,sep="")
  lena_segments = read.csv(paste(recording_dir,"/",recording,"_segments.csv",sep=""))
  chn_segments = subset(lena_segments,segtype=="CHNSP")
  adu_segments = subset(lena_segments,segtype=="FAN"|segtype=="MAN")
  rec_length = floor(lena_segments$endsec[nrow(lena_segments)])
  chn_voc_record = integer(rec_length)
  chn_voc_record[as.integer(floor(chn_segments$startsec))] = 1
  adu_voc_record = integer(rec_length)
  adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1
  
  ###################################################
  # Scale the simulation and human recording measures
  ###################################################
  
  # To-do: Load in the simulation data
  
  fitOrder_noTurns[,recording] = order(simFits_noTurns[,recording])
  fitOrder_wTurns[,recording] = order(simFits_wTurns[,recording])
  fitOrder_onlyTurns[,recording] = order(simFits_onlyTurns[,recording])
  
}

### subset sims_df by type of simulation
sims_df_nonInteractive = subset(sims_df,(chn_sim_othersensitivity==1&adu_sim_othersensitivity==1))
sims_df_a2Interactive = subset(sims_df,(chn_sim_othersensitivity==1&adu_sim_othersensitivity!=1))
sims_df_bidirectional = subset(sims_df,(chn_sim_othersensitivity!=1&adu_sim_othersensitivity!=1))

### subset fitOrder_noTurns by type of simulation
fitOrder_noTurns_nonInteractive = subset(fitOrder_noTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity==1))
fitOrder_noTurns_a2Interactive = subset(fitOrder_noTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity!=1))
fitOrder_noTurns_bidirectional = subset(fitOrder_noTurns,(sims_df$chn_sim_othersensitivity!=1&sims_df$adu_sim_othersensitivity!=1))

### Find out what the fits are for the three simulation types for each recording (lower is better)
colMeans(cbind(simFits_noTurns[fitOrder_noTurns_nonInteractive[1:5,1],1],simFits_noTurns[fitOrder_noTurns_nonInteractive[1:5,2],2],simFits_noTurns[fitOrder_noTurns_nonInteractive[1:5,3],3]))
colMeans(cbind(simFits_noTurns[fitOrder_noTurns_a2Interactive[1:5,1],1],simFits_noTurns[fitOrder_noTurns_a2Interactive[1:5,2],2],simFits_noTurns[fitOrder_noTurns_a2Interactive[1:5,3],3]))
colMeans(cbind(simFits_noTurns[fitOrder_noTurns_bidirectional[1:5,1],1],simFits_noTurns[fitOrder_noTurns_bidirectional[1:5,2],2],simFits_noTurns[fitOrder_noTurns_bidirectional[1:5,3],3]))

### subset fitOrder_wTurns by type of simulation
fitOrder_wTurns_nonInteractive = subset(fitOrder_wTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity==1))
fitOrder_wTurns_a2Interactive = subset(fitOrder_wTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity!=1))
fitOrder_wTurns_bidirectional = subset(fitOrder_wTurns,(sims_df$chn_sim_othersensitivity!=1&sims_df$adu_sim_othersensitivity!=1))

### Find out what the fits are for the three simulation types for each recording (lower is better)
colMeans(cbind(simFits_wTurns[fitOrder_wTurns_nonInteractive[1:5,1],1],simFits_wTurns[fitOrder_wTurns_nonInteractive[1:5,2],2],simFits_wTurns[fitOrder_wTurns_nonInteractive[1:5,3],3]))
colMeans(cbind(simFits_wTurns[fitOrder_wTurns_a2Interactive[1:5,1],1],simFits_wTurns[fitOrder_wTurns_a2Interactive[1:5,2],2],simFits_wTurns[fitOrder_wTurns_a2Interactive[1:5,3],3]))
colMeans(cbind(simFits_wTurns[fitOrder_wTurns_bidirectional[1:5,1],1],simFits_wTurns[fitOrder_wTurns_bidirectional[1:5,2],2],simFits_wTurns[fitOrder_wTurns_bidirectional[1:5,3],3]))

### subset fitOrder_onlyTurns by type of simulation
fitOrder_onlyTurns_nonInteractive = subset(fitOrder_onlyTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity==1))
fitOrder_onlyTurns_a2Interactive = subset(fitOrder_onlyTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity!=1))
fitOrder_onlyTurns_bidirectional = subset(fitOrder_onlyTurns,(sims_df$chn_sim_othersensitivity!=1&sims_df$adu_sim_othersensitivity!=1))

### Find out what the fits are for the three simulation types for each recording (lower is better)
colMeans(cbind(simFits_onlyTurns[fitOrder_onlyTurns_nonInteractive[1:5,1],1],simFits_onlyTurns[fitOrder_onlyTurns_nonInteractive[1:5,2],2],simFits_onlyTurns[fitOrder_onlyTurns_nonInteractive[1:5,3],3]))
colMeans(cbind(simFits_onlyTurns[fitOrder_onlyTurns_a2Interactive[1:5,1],1],simFits_onlyTurns[fitOrder_onlyTurns_a2Interactive[1:5,2],2],simFits_onlyTurns[fitOrder_onlyTurns_a2Interactive[1:5,3],3]))
colMeans(cbind(simFits_onlyTurns[fitOrder_onlyTurns_bidirectional[1:5,1],1],simFits_onlyTurns[fitOrder_onlyTurns_bidirectional[1:5,2],2],simFits_onlyTurns[fitOrder_onlyTurns_bidirectional[1:5,3],3]))

#################################################################################
# for the best-matched simulations,
# analyze the chn and adu ivis to test for response effects on ivis
# using 1. no control for previous ivi and 2. control for previous ivi or 2 ivis
#################################################################################

load("~/Documents/GitHub/vocal-response-analysis-simulation/mergedSimData_20250521.Rdat")

# simID = 0
# ivi_records = c()
# ivi_r_records = c()
# simIDs=c()
# previvi_resids = c()
# prev3ivi_resids = c()
# 
# for (fitRank in 1:20){
#   
#   simID = simID + 1
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