#########################################################################
# Randomly search simulation parameters to find best match to human data
#########################################################################

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

# 10th, 50th, and 90th percentiles of the child and adult IEI distributions