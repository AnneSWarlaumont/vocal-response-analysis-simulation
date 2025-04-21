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
chn_voc_record = numeric(rec_length)
chn_voc_record[as.integer(floor(chn_segments$startsec))] = 1
adu_voc_record = numeric(rec_length)
adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1

####################################
# Get key measures for human sample 
####################################

# Number of child and adult vocalizations
chn_total = sum(chn_voc_record)
adu_total = sum(adu_voc_record)

# Number of conversational turns based on a 5 s threshold
# Initialize leader and follower ID and time
# Loop through each second of simulation
# When a turn is identified, update the turn count

# 10th, 50th, and 90th percentiles of the child and adult IEI distributions