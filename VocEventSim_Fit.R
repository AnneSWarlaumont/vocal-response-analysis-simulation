# get_ivis <- function(voc_record){
#   ivi_record = c()
#   voc1_t = NA # initialize the time of the first vocalization in the current IVI pair
#   voc2_t = NA # initialize the time of the second vocalization in the current IVI pair 
#   for (t in 1:length(voc_record)){
#     if (is.na(voc1_t) & voc_record[t]==1){
#       voc1_t = t
#     } else if (voc_record[t]==1){
#       voc2_t = t
#       ivi = voc2_t - voc1_t
#       ivi_record = c(ivi_record,ivi)
#       voc1_t = voc2_t
#     }
#   }
#   return(ivi_record)
# }

########################################################################
# compute the simulations' measures to be used for fitting to human data
########################################################################
# 
# a1_total = sum(a1_voc_record)
# a2_total = sum(a2_voc_record)
# 
# sim_turnCount = 0
# leadID = 'none'
# leadT = 0
# for (t in 1:sim_length){
#   leadT = leadT+1
#   if (leadID=='a1' && a2_voc_record[t]==1) {
#     if (leadT <= 5) {
#       sim_turnCount = sim_turnCount+1
#     }
#     leadID = 'a2'
#     leadT = 0
#   } else if (leadID=='a2' && a1_voc_record[t]==1) {
#     if (leadT <= 5) {
#       sim_turnCount = sim_turnCount+1
#     }
#     leadID = 'a1'
#     leadT = 0
#   } else if (leadID=='none' && a2_voc_record[t]==1) {
#     leadID = 'a2'
#     leadT = 0
#   } else if (leadID=='none' && a1_voc_record[t]==1) {
#     leadID = 'a1'
#     leadT = 0
#   }
#   
# }
# 
# a1_ivi_record = get_ivis(a1_voc_record)
# a2_ivi_record = get_ivis(a2_voc_record)
# a1_n_ivi_1 = sum(a1_ivi_record==1)
# a2_n_ivi_1 = sum(a2_ivi_record==1)
# a1_ivi_50 = quantile(a1_ivi_record,probs=.5,names=FALSE)
# a2_ivi_50 = quantile(a2_ivi_record,probs=.5,names=FALSE)
# a1_ivi_max = max(a1_ivi_record)
# a2_ivi_max = max(a2_ivi_record)

recordings = c("0344_000913","0833_010606","0054_000603","0196_000902","0274_000221","0300_000607","0437_010603")

# ######################################
# # Read in and format the human sample
# ######################################
# 
# recording_dir = paste("data/",recording,sep="")
# lena_segments = read.csv(paste(recording_dir,"/",recording,"_segments.csv",sep=""))
# chn_segments = subset(lena_segments,segtype=="CHNSP")
# adu_segments = subset(lena_segments,segtype=="FAN"|segtype=="MAN")
# rec_length = floor(lena_segments$endsec[nrow(lena_segments)])
# chn_voc_record = integer(rec_length)
# chn_voc_record[as.integer(floor(chn_segments$startsec))] = 1
# adu_voc_record = integer(rec_length)
# adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1
# 
# ####################################
# # Get key measures for human sample 
# ####################################
# 
# # Number of child and adult vocalizations
# chn_total = sum(chn_voc_record)
# adu_total = sum(adu_voc_record)
# 
# # Number of conversational turns based on a 5 s threshold
# # Initialize leader ID and time
# # Loop through each second of the recording
# # When a turn is identified, update the turn count
# turnCount = 0
# leadID = 'none'
# leadT = 0
# for (t in 1:rec_length){
#   leadT = leadT+1
#   if (leadID=='chn' && adu_voc_record[t]==1) {
#     if (leadT <= 5) {
#       turnCount = turnCount+1
#     }
#     leadID = 'adu'
#     leadT = 0
#   } else if (leadID=='adu' && chn_voc_record[t]==1) {
#     if (leadT <= 5) {
#       turnCount = turnCount+1
#     }
#     leadID = 'chn'
#     leadT = 0
#   } else if (leadID=='none' && adu_voc_record[t]==1) {
#     leadID = 'adu'
#     leadT = 0
#   } else if (leadID=='none' && chn_voc_record[t]==1) {
#     leadID = 'chn'
#     leadT = 0
#   }
# }
# 
# # 10th, 50th, and 90th percentiles of the child and adult ivi distributions
# # try the quantile function
# 
# chn_ivi_record = get_ivis(chn_voc_record)
# adu_ivi_record = get_ivis(adu_voc_record)
# chn_ivi_10 = quantile(chn_ivi_record,probs=.1,names=FALSE)
# chn_ivi_50 = quantile(chn_ivi_record,probs=.5,names=FALSE)
# chn_ivi_90 = quantile(chn_ivi_record,probs=.9,names=FALSE)
# chn_ivi_99 = quantile(chn_ivi_record,probs=.99,names=FALSE)
# adu_ivi_10 = quantile(adu_ivi_record,probs=.1,names=FALSE)
# adu_ivi_50 = quantile(adu_ivi_record,probs=.5,names=FALSE)
# adu_ivi_90 = quantile(adu_ivi_record,probs=.9,names=FALSE)
# adu_ivi_99 = quantile(adu_ivi_record,probs=.99,names=FALSE)

# #################################################
# # Assess the simulations' fits to the human data
# # (In the future, move this into the for loop,
# # so as to be computed after each simulation)
# #################################################
# 
# chn_total_logplus1 = log(chn_total+1)
# sims_df$chn_sim_total_logplus1 = log(sims_df$chn_sim_total+1)
# sims_df$chn_total_diff = (sims_df$chn_sim_total_logplus1-chn_total_logplus1)/chn_total_logplus1
# 
# chn_ivi_50_logplus1 = log(chn_ivi_50+1)
# sims_df$chn_sim_ivi_50_logplus1 = log(sims_df$chn_sim_ivi_50+1)
# sims_df$chn_ivi_50_diff = (sims_df$chn_sim_ivi_50_logplus1-chn_ivi_50_logplus1)/chn_ivi_50_logplus1
# 
# chn_ivi_max_logplus1 = log(max(chn_ivi_record)+1)
# sims_df$chn_sim_ivi_max_logplus1 = log(sims_df$chn_sim_ivi_max+1)
# sims_df$chn_ivi_max_diff = (sims_df$chn_sim_ivi_max_logplus1-chn_ivi_max_logplus1)/chn_ivi_max_logplus1
# 
# chn_n_ivi_1 = sum(chn_ivi_record==1)
# chn_n_ivi_1_logplus1 = log(chn_n_ivi_1+1)
# sims_df$chn_sim_n_ivi_1_logplus1 = log(sims_df$chn_sim_n_ivi_1+1)
# sims_df$chn_n_ivi_1_diff = (sims_df$chn_sim_n_ivi_1_logplus1-chn_n_ivi_1_logplus1)/chn_n_ivi_1_logplus1
# 
# adu_total_logplus1 = log(adu_total+1)
# sims_df$adu_sim_total_logplus1 = log(sims_df$adu_sim_total+1)
# sims_df$adu_total_diff = (sims_df$adu_sim_total_logplus1-adu_total_logplus1)/adu_total_logplus1
# 
# adu_ivi_50_logplus1 = log(adu_ivi_50+1)
# sims_df$adu_sim_ivi_50_logplus1 = log(sims_df$adu_sim_ivi_50+1)
# sims_df$adu_ivi_50_diff = (sims_df$adu_sim_ivi_50_logplus1-adu_ivi_50_logplus1)/adu_ivi_50_logplus1
# 
# adu_ivi_max_logplus1 = log(max(adu_ivi_record)+1)
# sims_df$adu_sim_ivi_max_logplus1 = log(sims_df$adu_sim_ivi_max+1)
# sims_df$adu_ivi_max_diff = (sims_df$adu_sim_ivi_max_logplus1-adu_ivi_max_logplus1)/adu_ivi_max_logplus1
# 
# adu_n_ivi_1 = sum(adu_ivi_record==1)
# adu_n_ivi_1_logplus1 = log(adu_n_ivi_1+1)
# sims_df$adu_sim_n_ivi_1_logplus1 = log(sims_df$adu_sim_n_ivi_1+1)
# sims_df$adu_n_ivi_1_diff = (sims_df$adu_sim_n_ivi_1_logplus1-adu_n_ivi_1_logplus1)/adu_n_ivi_1_logplus1
# 
# turnCount_logplus1 = log(turnCount+1)
# sims_df$sim_turnCount_logplus1 = log(sims_df$sim_turnCount+1)
# sims_df$turnCount_diff = (sims_df$sim_turnCount_logplus1-turnCount_logplus1)/turnCount_logplus1
# 
# # Get fit without considering turn count
# sims_df$simDist_noTurns = sqrt((sims_df$chn_total_diff)^2
#                                +(sims_df$chn_ivi_50_diff)^2
#                                +(sims_df$chn_ivi_max_diff)^2
#                                +(sims_df$chn_n_ivi_1_diff)^2
#                                +(sims_df$adu_total_diff)^2
#                                +(sims_df$adu_ivi_50_diff)^2
#                                +(sims_df$adu_ivi_max_diff)^2
#                                +(sims_df$adu_n_ivi_1_diff)^2)
# 
# # Get fit including turn count
# sims_df$simDist_wTurns = sqrt((sims_df$chn_total_diff)^2
#                               +(sims_df$chn_ivi_50_diff)^2
#                               +(sims_df$chn_ivi_max_diff)^2
#                               +(sims_df$chn_n_ivi_1_diff)^2
#                               +(sims_df$adu_total_diff)^2
#                               +(sims_df$adu_ivi_50_diff)^2
#                               +(sims_df$adu_ivi_max_diff)^2
#                               +(sims_df$adu_n_ivi_1_diff)^2
#                               +(sims_df$turnCount_diff)^2)
# 
# # Get fit based only on the turn count
# sims_df$simDist_onlyTurns = sqrt((sims_df$turnCount_diff)^2)
# 
# fitOrder_noTurns = order(sims_df$simDist_noTurns)
# fitOrder_wTurns = order(sims_df$simDist_wTurns)
# fitOrder_onlyTurns = order(sims_df$simDist_onlyTurns)
