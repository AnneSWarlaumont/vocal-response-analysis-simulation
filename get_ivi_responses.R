# A function that, given a series of IVIs by one agent, a record of
# vocalizations by a second agent, and a response time threshold, returns a
# record of which agent 1 IVIs corresponded to the first vocalization in the IVI
# were followed by an agent 2 vocalization with onset within the response time
# threshold
get_ivi_responses <- function(a1_ivis,a2_voc_record,rthresh){
  ivi_responses <- c()
  t <- 1
  for (i in 1:length(a1_ivis)){
    if (a1_ivis[i]<=rthresh){
      ivi_responses[i] <- NA
    } else if (sum(a2_voc_record[(t+1):(t+rthresh)],na.rm = TRUE)>0){
      ivi_responses[i] <- 1
    } else{
      ivi_responses[i] <- 0
    }
    t <- t+a1_ivis[i]
  }
  return(ivi_responses)
}