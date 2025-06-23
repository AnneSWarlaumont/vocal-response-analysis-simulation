library(dplyr)

recordings = c("0054_000603","0196_000902","0274_000221","0300_000607","0344_000913","0437_010603","0833_010606")
recordings = c("0054_000603") # debugging

for (recording in recordings){
  
  ######################################
  # Read in and format the human sample
  ######################################
  
  recording_dir = paste("data/",recording,sep="")
  lena_segments = read.csv(paste(recording_dir,"/",recording,"_segments.csv",sep=""))
  chn_segments = subset(lena_segments,segtype=="CHNSP")
  adu_segments = subset(lena_segments,segtype=="FAN"|segtype=="MAN")
  rec_length = floor(lena_segments$endsec[nrow(lena_segments)])
  chn_voc_record = integer(rec_length)
  chn_voc_record[as.integer(floor(chn_segments$startsec))] = 1
  adu_voc_record = integer(rec_length)
  adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1 
  
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
  
  get_ivis_continuous <- function(segments,rthresh){
    ivi_record = c()
    ivi_r_record = c()
    voc1_endt = NA
    voc2_startt = NA
    for (n in 1:(nrow(segments)-1)){
      voc1_endt = segments[n,]$endsec
      voc2_startt = segments[n+1,]$startsec
      ivi = voc2_startt - voc1_endt
      ivi_record = c(ivi_record,ivi)
      if (ivi<=rthresh){
        ivi_r_record[n] = NA 
      } else if (any(between(segments$startsec,voc1_endt,voc1_endt+rthresh))){
        ivi_r_record[n] = 1
      } else{
        ivi_r_record[n] = 0
      }
    }
    return(ivi_record,ivi_r_record)
  }

  chn_ivi_record = get_ivis(chn_voc_record)
  adu_ivi_record = get_ivis(adu_voc_record)
  
  chn_continuous = get_ivis_continuous(chn_segments)
  chn_ivi_record_continuous = chn_continuous$ivi_record
  chn_ivi_r_record_continuous = chn_continuous$ivi_r_record
  
  adu_continuous = get_ivis_continuous(adu_segments)
  adu_ivi_record_continuous = adu_continuous$ivi_record
  adu_ivi_r_record_continuous = adu_continuous$ivi_r_record
  
  hist(log(chn_ivi_record+1))
  hist(log(chn_ivi_record_continuous+1))
  
}
