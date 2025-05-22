setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

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
  
  #######################################
  # Get key measures for the human sample 
  #######################################
  
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
  
  # 50th percentiles of the child and adult ivi distributions using the quantile function
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
  chn_ivi_50 = quantile(chn_ivi_record,probs=.5,names=FALSE)
  adu_ivi_50 = quantile(adu_ivi_record,probs=.5,names=FALSE)
  
  ###################################################
  # Scale the simulation and human recording measures
  ###################################################
  
  chn_sim_total_scaled = scale(log(sims_df$chn_sim_total+1))
  chn_total_scaled = (log(chn_total+1)-attr(chn_sim_total_scaled,"scaled:center"))/attr(chn_sim_total_scaled,"scaled:scale")
  
  chn_sim_ivi_50_scaled = scale(log(sims_df$chn_sim_ivi_50))
  chn_ivi_50_scaled = (log(chn_ivi_50)-attr(chn_sim_ivi_50_scaled,"scaled:center"))/attr(chn_sim_ivi_50_scaled,"scaled:scale")
  
  chn_ivi_max = max(chn_ivi_record)
  sims_df$chn_sim_ivi_max[which(sims_df$chn_sim_ivi_max==-Inf)] = 0
  chn_sim_ivi_max_scaled = scale(log(sims_df$chn_sim_ivi_max+1))
  chn_ivi_max_scaled = (log(chn_ivi_max+1)-attr(chn_sim_ivi_max_scaled,"scaled:center"))/attr(chn_sim_ivi_max_scaled,"scaled:scale")
  
  chn_n_ivi_1 = sum(chn_ivi_record==1)
  chn_sim_n_ivi_1_scaled = scale(log(sims_df$chn_sim_n_ivi_1+1))
  chn_n_ivi_1_scaled = (log(chn_n_ivi_1+1)-attr(chn_sim_n_ivi_1_scaled,"scaled:center"))/attr(chn_sim_n_ivi_1_scaled,"scaled:scale")
  
  adu_sim_total_scaled = scale(log(sims_df$adu_sim_total+1))
  adu_total_scaled = (log(adu_total+1)-attr(adu_sim_total_scaled,"scaled:center"))/attr(adu_sim_total_scaled,"scaled:scale")
  
  adu_sim_ivi_50_scaled = scale(log(sims_df$adu_sim_ivi_50))
  adu_ivi_50_scaled = (log(adu_ivi_50)-attr(adu_sim_ivi_50_scaled,"scaled:center"))/attr(adu_sim_ivi_50_scaled,"scaled:scale")
  
  adu_ivi_max = max(adu_ivi_record)
  sims_df$adu_sim_ivi_max[which(sims_df$adu_sim_ivi_max==-Inf)] = 0
  adu_sim_ivi_max_scaled = scale(log(sims_df$adu_sim_ivi_max+1))
  adu_ivi_max_scaled = (log(adu_ivi_max+1)-attr(adu_sim_ivi_max_scaled,"scaled:center"))/attr(adu_sim_ivi_max_scaled,"scaled:scale")
  
  adu_n_ivi_1 = sum(adu_ivi_record==1)
  adu_sim_n_ivi_1_scaled = scale(log(sims_df$adu_sim_n_ivi_1+1))
  adu_n_ivi_1_scaled = (log(adu_n_ivi_1+1)-attr(adu_sim_n_ivi_1_scaled,"scaled:center"))/attr(adu_sim_n_ivi_1_scaled,"scaled:scale")
  
  sim_turnCount_scaled = scale(log(sims_df$sim_turnCount+1))
  turnCount_scaled = (log(turnCount+1)-attr(sim_turnCount_scaled,"scaled:center"))/attr(sim_turnCount_scaled,"scaled:scale")
  
  # Get fit without considering turn count
  simFits_noTurns[,recording] = sqrt((chn_sim_total_scaled-chn_total_scaled)^2+(chn_sim_n_ivi_1_scaled-chn_n_ivi_1_scaled)^2+(chn_sim_ivi_50_scaled-chn_ivi_50_scaled)^2+(chn_sim_ivi_max_scaled-chn_ivi_max_scaled)^2+(adu_sim_total_scaled-adu_total_scaled)^2+(adu_sim_n_ivi_1_scaled-adu_n_ivi_1_scaled)^2+(adu_sim_ivi_50_scaled-adu_ivi_50_scaled)^2+(adu_sim_ivi_max_scaled-adu_ivi_max_scaled)^2)
  
  # Get fit including turn count
  simFits_wTurns[,recording] = sqrt((chn_sim_total_scaled-chn_total_scaled)^2+(chn_sim_n_ivi_1_scaled-chn_n_ivi_1_scaled)^2+(chn_sim_ivi_50_scaled-chn_ivi_50_scaled)^2+(chn_sim_ivi_max_scaled-chn_ivi_max_scaled)^2+(adu_sim_total_scaled-adu_total_scaled)^2+(adu_sim_n_ivi_1_scaled-adu_n_ivi_1_scaled)^2+(adu_sim_ivi_50_scaled-adu_ivi_50_scaled)^2+(adu_sim_ivi_max_scaled-adu_ivi_max_scaled)^2+(sim_turnCount_scaled-turnCount_scaled)^2)
  
  # Get fit focusing only on the turn count
  simFits_onlyTurns[,recording] = sqrt((sim_turnCount_scaled-turnCount_scaled)^2)
  
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

### Find out what the fitRanks are for the three simulation types for each recording
colMeans(cbind(simFits_noTurns[fitOrder_noTurns_nonInteractive[1:5,1],1],simFits_noTurns[fitOrder_noTurns_nonInteractive[1:5,2],2],simFits_noTurns[fitOrder_noTurns_nonInteractive[1:5,3],3]))
colMeans(cbind(simFits_noTurns[fitOrder_noTurns_a2Interactive[1:5,1],1],simFits_noTurns[fitOrder_noTurns_a2Interactive[1:5,2],2],simFits_noTurns[fitOrder_noTurns_a2Interactive[1:5,3],3]))
colMeans(cbind(simFits_noTurns[fitOrder_noTurns_bidirectional[1:5,1],1],simFits_noTurns[fitOrder_noTurns_bidirectional[1:5,2],2],simFits_noTurns[fitOrder_noTurns_bidirectional[1:5,3],3]))

### subset fitOrder_wTurns by type of simulation
fitOrder_wTurns_nonInteractive = subset(fitOrder_wTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity==1))
fitOrder_wTurns_a2Interactive = subset(fitOrder_wTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity!=1))
fitOrder_wTurns_bidirectional = subset(fitOrder_wTurns,(sims_df$chn_sim_othersensitivity!=1&sims_df$adu_sim_othersensitivity!=1))

### Find out what the fitRanks are for the three simulation types for each recording
colMeans(cbind(simFits_wTurns[fitOrder_wTurns_nonInteractive[1:5,1],1],simFits_wTurns[fitOrder_wTurns_nonInteractive[1:5,2],2],simFits_wTurns[fitOrder_wTurns_nonInteractive[1:5,3],3]))
colMeans(cbind(simFits_wTurns[fitOrder_wTurns_a2Interactive[1:5,1],1],simFits_wTurns[fitOrder_wTurns_a2Interactive[1:5,2],2],simFits_wTurns[fitOrder_wTurns_a2Interactive[1:5,3],3]))
colMeans(cbind(simFits_wTurns[fitOrder_wTurns_bidirectional[1:5,1],1],simFits_wTurns[fitOrder_wTurns_bidirectional[1:5,2],2],simFits_wTurns[fitOrder_wTurns_bidirectional[1:5,3],3]))

### subset fitOrder_onlyTurns by type of simulation
fitOrder_onlyTurns_nonInteractive = subset(fitOrder_onlyTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity==1))
fitOrder_onlyTurns_a2Interactive = subset(fitOrder_onlyTurns,(sims_df$chn_sim_othersensitivity==1&sims_df$adu_sim_othersensitivity!=1))
fitOrder_onlyTurns_bidirectional = subset(fitOrder_onlyTurns,(sims_df$chn_sim_othersensitivity!=1&sims_df$adu_sim_othersensitivity!=1))

### Find out what the fitRanks are for the three simulation types for each recording
colMeans(cbind(simFits_onlyTurns[fitOrder_onlyTurns_nonInteractive[1:5,1],1],simFits_onlyTurns[fitOrder_onlyTurns_nonInteractive[1:5,2],2],simFits_onlyTurns[fitOrder_onlyTurns_nonInteractive[1:5,3],3]))
colMeans(cbind(simFits_onlyTurns[fitOrder_onlyTurns_a2Interactive[1:5,1],1],simFits_onlyTurns[fitOrder_onlyTurns_a2Interactive[1:5,2],2],simFits_onlyTurns[fitOrder_onlyTurns_a2Interactive[1:5,3],3]))
colMeans(cbind(simFits_onlyTurns[fitOrder_onlyTurns_bidirectional[1:5,1],1],simFits_onlyTurns[fitOrder_onlyTurns_bidirectional[1:5,2],2],simFits_onlyTurns[fitOrder_onlyTurns_bidirectional[1:5,3],3]))


#### Left off here

# #################################################################################
# # for the best-matched simulations,
# # analyze the chn and adu ivis to test for response effects on ivis
# # using 1. no control for previous ivi and 2. control for previous ivi or 2 ivis
# #################################################################################
# 
# rthresh = 5
# 
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