# A function that calculates IVIs based on a vocal events series
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