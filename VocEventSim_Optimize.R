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

sim_length = length(chn_voc_record)

# Start with a child-only simulation
minp_range = c(.000001,.01)
maxp_range = c(.01,1)
sdlog_range = c(.01,1)

for (sim in 1:50){
  a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
  a1_minp = runif(1,min=minp_range[1],max=minp_range[2])
  a1_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
  a1_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
  
  # run the sim
  # save the parameters in a data frame
}

# find out what parameter combinations are the best (and maybe 2nd through 5th best match).
# Consider different weightings for different measures?
# Plot the best match(es)'s ivi distributions (on log-log), compared to human data
# And plot raw data points at three timescales, compared to human data
