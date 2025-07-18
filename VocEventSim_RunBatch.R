setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

simTypes = c("nonInteractive","a2interactive","bidirectional")

nSims = 10000


for (simType in simTypes){
    
    #########################################################################
    # Randomly choose simulation parameters and run simulation
    #########################################################################
    
    source('VocEventSim.R')
    
    sim_length = (1+16)*60*60 # 16 hours + an hour of warm-up
    
    meanlog_range = c(-.05,0)
    minp_range = c(.000000001,.00001)
    maxp_range = c(.05,.5)
    sdlog_range = c(.01,1)
    othersensitivity_range = c(1,10)
    
    sims_df = data.frame(chn_sim_minp = double(),
                         chn_sim_maxp = double(),
                         chn_sim_meanlog = double(),
                         chn_sim_sdlog = double(),
                         chn_sim_othersensitivity = double(),
                         adu_sim_minp = double(),
                         adu_sim_maxp = double(),
                         adu_sim_meanlog = double(),
                         adu_sim_sdlog = double(),
                         adu_sim_othersensitivity = double(),
                         stringsAsFactors = FALSE)
    sims_chn_voc_records <- list()
    sims_adu_voc_records <- list()
    
    for (simNum in 1:nSims){
      
      a1_minp = runif(1,min=minp_range[1],max=minp_range[2])
      a2_minp = runif(1,min=minp_range[1],max=minp_range[2])
      
      a1_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
      a2_maxp = runif(1,min=maxp_range[1],max=maxp_range[2])
      
      a1_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
      a2_sdlog = runif(1,min=sdlog_range[1],max=sdlog_range[2])
      
      if (simType=="nonInteractive") {
        a1_meanlog = 0
        a2_meanlog = 0
        a1_othersensitivity = 1
        a2_othersensitivity = 1
      } else if (simType=="a2interactive") {
        a1_meanlog = 0
        a2_meanlog = runif(1,min=meanlog_range[1],max=meanlog_range[2])
        a1_othersensitivity = 1
        a2_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
      } else if (simType=="bidirectional") {
        a1_meanlog = runif(1,min=meanlog_range[1],max=meanlog_range[2])
        a2_meanlog = runif(1,min=meanlog_range[1],max=meanlog_range[2])
        a1_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
        a2_othersensitivity = runif(1,min=othersensitivity_range[1],max=othersensitivity_range[2])
      }
      
      a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
      a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)
      
      # run the sim
      voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,1,1,5)
      a1_voc_record = voc_records[[1]]
      a2_voc_record = voc_records[[2]]
      
      # save the simulation parameters and measures to a data frame
      df_row = data.frame(chn_sim_minp = a1_minp,
                          chn_sim_maxp = a1_maxp,
                          chn_sim_meanlog = a1_meanlog,
                          chn_sim_sdlog = a1_sdlog,
                          chn_sim_othersensitivity = a1_othersensitivity,
                          adu_sim_minp = a2_minp,
                          adu_sim_maxp = a2_maxp,
                          adu_sim_meanlog = a2_meanlog,
                          adu_sim_sdlog = a2_sdlog,
                          adu_sim_othersensitivity = a2_othersensitivity)
      
      # save the vocalization and ivi records to lists
      sims_chn_voc_records[[simNum]] <- a1_voc_record
      sims_adu_voc_records[[simNum]] <- a2_voc_record
      sims_df = rbind(sims_df,df_row)
    
    }
        
    ###########################
    # save the simulation data
    ###########################
    
    results_dir = paste("data/",simType,sep="")
    
    if (!dir.exists(results_dir)){
      dir.create(results_dir)
    }
    
    # save the full workspace
    save(list = c("sims_df","sims_chn_voc_records","sims_adu_voc_records"),
               file = paste(results_dir,"/VocEventSim_Batch",".RData",sep=""))
}

