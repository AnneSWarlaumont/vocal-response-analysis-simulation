two_agent_vocal_sim <- function(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity) {

  a1_voc_record = c()
  a2_voc_record = c()

  for (t in 1:sim_length){

    a1_voc_event = rbinom(1,1,a1_p_voc)
    a1_voc_record = c(a1_voc_record,a1_voc_event)
    a1_p_voc_multiplier = rlnorm(1,meanlog=a1_meanlog,sdlog=a1_sdlog)
    a1_p_voc = max(a1_minp,min(a1_maxp,a1_p_voc_multiplier*a1_p_voc))

    if(a1_voc_event==1){ # For now, it is assumed that if a1 vocalizes, a2 does not.
      a2_voc_event = 0
      a1_p_voc = a1_p_voc*a1_p_voc_multiplier
    } else{
      a2_voc_event = rbinom(1,1,a2_p_voc)
      a2_p_voc_multiplier = rlnorm(1,meanlog=a2_meanlog,sdlog=a2_sdlog)
    }
    a2_voc_record = c(a2_voc_record,a2_voc_event)
    a2_p_voc = max(a2_minp,min(a2_maxp,a2_p_voc_multiplier*a2_p_voc))

  }

  return(list(a1_voc_record,a2_voc_record))

}

get_responses <- function(a1_voc_record,a2_voc_record,rthresh){
  
  sim_length = length(a1_voc_record)
  a2toa1_r_record = c()
  ## I left off here.
  
  return(a2toa1_r_record)
  
}

voc_records = two_agent_vocal_sim(10*60*60,.02,.02,.02,.02,.2,.2,.002,.002,.2,.2,1,100,1,1)
a1_voc_record = voc_records[[1]]
a2_voc_record = voc_records[[2]]
rthresh = 1
a2toa1_r_record = get_responses(a1_voc_record,a2_voc_record,rthresh)
