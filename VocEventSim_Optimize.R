setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

######################################
# Read in and format the human sample
######################################

lena_segments = read.csv("0344_000913_segments.csv")
chn_segments = subset(lena_segments,segtype=="CHNSP"|segtype=="CHNNSP")
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
adu_ivi_10 = quantile(adu_ivi_record,probs=.1,names=FALSE)
adu_ivi_50 = quantile(adu_ivi_record,probs=.5,names=FALSE)
adu_ivi_90 = quantile(adu_ivi_record,probs=.9,names=FALSE)

#########################################################################
# Randomly search simulation parameters to find best match to human data
#########################################################################

source('VocEventSim.R')

sim_length = length(chn_voc_record)

rthresh = 5
a1_meanlog = 0
a2_meanlog = 0

minp_range = c(.000001,.01)
maxp_range = c(.01,.5)
sdlog_range = c(.01,.5)
othersensitivity_range = c(.5,3)
respsensitivity_range = c(.5,3)

sims_df = data.frame(chn_sim_minp = double(),
                     chn_sim_maxp = double(),
                     chn_sim_sdlog = double(),
                     chn_sim_othersensitivity = double(),
                     chn_sim_respsensitivity = double(),
                     chn_sim_total = double(),
                     chn_sim_ivi_10 = double(),
                     chn_sim_ivi_50 = double(),
                     chn_sim_ivi_90 = double(),
                     adu_sim_minp = double(),
                     adu_sim_maxp = double(),
                     adu_sim_sdlog = double(),
                     adu_sim_othersensitivity = double(),
                     adu_sim_respsensitivity = double(),
                     adu_sim_total = double(),
                     adu_sim_ivi_10 = double(),
                     adu_sim_ivi_50 = double(),
                     adu_sim_ivi_90 = double(),
                     sim_turncount = double(),
                     stringsAsFactors = FALSE)

for (simNum in 1:10000){
  
  a1_minp = runif(1,min=minp_range[1],max=minp_range[2])
  a1_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
  a1_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
  a1_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
  a1_respsensitivity = runif(1,min=respsensitivity_range[1],max=respsensitivity_range[2])
  
  a2_minp = runif(1,min=minp_range[1],max=minp_range[2])
  a2_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
  a2_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
  a2_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
  a2_respsensitivity = runif(1,min=respsensitivity_range[1],max=respsensitivity_range[2])
  
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
  a1_ivi_10 = quantile(a1_ivi_record,probs=.1,names=FALSE)
  a1_ivi_50 = quantile(a1_ivi_record,probs=.5,names=FALSE)
  a1_ivi_90 = quantile(a1_ivi_record,probs=.9,names=FALSE)
  a2_ivi_10 = quantile(a2_ivi_record,probs=.1,names=FALSE)
  a2_ivi_50 = quantile(a2_ivi_record,probs=.5,names=FALSE)
  a2_ivi_90 = quantile(a2_ivi_record,probs=.9,names=FALSE)
  
  # match a1 and a2 with chn or adu based on which pairing provides the better fit
  # as an initial step, can try matching based on vocalization count
  # save the simulation parameters and measures to a data frame
  if ((chn_total < adu_total && a1_total < a2_total) || (chn_total > adu_total && a2_total > a1_total)){
    df_row = data.frame(chn_sim_minp = a1_minp,
                        chn_sim_maxp = a1_maxp,
                        chn_sim_sdlog = a1_sdlog,
                        chn_sim_othersensitivity = a1_othersensitivity,
                        chn_sim_respsensitivity = a1_respsensitivity,
                        chn_sim_total = a1_total,
                        chn_sim_ivi_10 = a1_ivi_10,
                        chn_sim_ivi_50 = a1_ivi_50,
                        chn_sim_ivi_90 = a1_ivi_90,
                        adu_sim_minp = a2_minp,
                        adu_sim_maxp = a2_maxp,
                        adu_sim_sdlog = a2_sdlog,
                        adu_sim_othersensitivity = a2_othersensitivity,
                        adu_sim_respsensitivity = a2_respsensitivity,
                        adu_sim_total = a2_total,
                        adu_sim_ivi_10 = a2_ivi_10,
                        adu_sim_ivi_50 = a2_ivi_50,
                        adu_sim_ivi_90 = a2_ivi_90,
                        sim_turnCount)
  } else {
    df_row = data.frame(chn_sim_minp = a2_minp,
                        chn_sim_maxp = a2_maxp,
                        chn_sim_sdlog = a2_sdlog,
                        chn_sim_othersensitivity = a2_othersensitivity,
                        chn_sim_respsensitivity = a2_respsensitivity,
                        chn_sim_total = a2_total,
                        chn_sim_ivi_10 = a2_ivi_10,
                        chn_sim_ivi_50 = a2_ivi_50,
                        chn_sim_ivi_90 = a2_ivi_90,
                        adu_sim_minp = a1_minp,
                        adu_sim_maxp = a1_maxp,
                        adu_sim_sdlog = a1_sdlog,
                        adu_sim_othersensitivity = a1_othersensitivity,
                        adu_sim_respsensitivity = a1_respsensitivity,
                        adu_sim_total = a1_total,
                        adu_sim_ivi_10 = a1_ivi_10,
                        adu_sim_ivi_50 = a1_ivi_50,
                        adu_sim_ivi_90 = a1_ivi_90,
                        sim_turnCount)
  }
  sims_df = rbind(sims_df,df_row)
}

# find out what parameter combinations are the best (and maybe 2nd through 5th best match).
# Consider different weightings for different measures?
# Plot the best match(es)'s ivi distributions (on log-log), compared to human data
# And plot raw data points at three timescales, compared to human data
