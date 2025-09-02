library(dplyr)
get_ivis_continuous <- function(a1_segments,a2_segments,rthresh){
  ivi_record = c()
  ivi_r_record = c()
  voc1_endt = NA
  voc2_startt = NA
  for (n in 1:(nrow(a1_segments)-1)){
    voc1_endt = a1_segments[n,]$endsec
    voc2_startt = a1_segments[n+1,]$startsec
    ivi = voc2_startt - voc1_endt
    if (ivi!=0){
      if (ivi<=rthresh){
        ivi_r = NA 
      } else if (any(between(a2_segments$startsec,voc1_endt,voc1_endt+rthresh))){
        ivi_r = 1
      } else{
        ivi_r = 0
      }
      ivi_record = c(ivi_record,ivi)
      ivi_r_record = c(ivi_r_record,ivi_r)
    }
  }
  return(list(ivi_record,ivi_r_record))
}