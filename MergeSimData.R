# After running the above simulations, I realized we don't need separate
# optimization runs for each daylong recording, since at this point the
# optimization is pure random search. So we can merge the three sets of
# simulation runs and then re-do the search for best fits.

simTypes2merge = c("nonInteractive","a2interactive","bidirectional")
recordings2merge = c("0344_000913","0833_010606","0054_000603")

for (simType2merge in simTypes2merge){
  for (recording2merge in recordings2merge){
    workspace_fileAndPath = Sys.glob(paste("~/Documents/GitHub/vocal-response-analysis-simulation/data/",recording2merge,"/",simType2merge,"/*",sep=""))
    load(workspace_fileAndPath)
    if (exists("sims_df_merged")){
      sims_df_merged = rbind(sims_df_merged,sims_df)
    } else{
      sims_df_merged = sims_df
    }
  }
}
sims_df_merged <- sims_df_merged[, -which(names(sims_df_merged) == "simDist")]
write.csv(sims_df_merged, file = "~/Documents/GitHub/vocal-response-analysis-simulation/sims_df_merged_20250521.csv")
