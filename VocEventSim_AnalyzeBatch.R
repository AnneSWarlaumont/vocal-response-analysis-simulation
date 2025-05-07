# This script assumes that key data frames and variables from a run of VocEventSim_Optimize are loaded.


    fitRank = 1
    
    library(ggplot2)
    hist(log(chn_ivi_record))
    hist(log(sims_chn_ivi_records[[fitOrder[fitRank]]]))
    hist(log(adu_ivi_record))
    hist(log(sims_adu_ivi_records[[fitOrder[fitRank]]]))
    
    library(poweRlaw)
    m = displ$new(chn_ivi_record)
    est = estimate_xmin(m)
    m$setXmin(est[[2]])
    m$setPars(est[[3]])
    plot(m)
    lines(m, col=2)
    
    m = displ$new(sims_chn_ivi_records[[fitOrder[fitRank]]])
    est = estimate_xmin(m)
    m$setXmin(est[[2]])
    m$setPars(est[[3]])
    plot(m)
    lines(m, col=2)
    
    m = displ$new(adu_ivi_record)
    est = estimate_xmin(m)
    m$setXmin(est[[2]])
    m$setPars(est[[3]])
    plot(m)
    lines(m, col=2)
    
    m = displ$new(sims_adu_ivi_records[[fitOrder[fitRank]]])
    est = estimate_xmin(m)
    m$setXmin(est[[2]])
    m$setPars(est[[3]])
    plot(m)
    lines(m, col=2)
    
    par(mfrow=c(3,1),cex=1,mar=c(1,1,2,1),oma=c(0,0,4,0))
    stripchart(which(sims_chn_voc_records[[fitOrder[fitRank]]]==1),xaxt="n",main="Full day simulation (10 hours)",pch=19,ylim=c(.5,1.5),xlim=c(0,sim_length))
    
    hr_offset = 0.2
    rect(xleft=hr_offset*60*60,xright=3600+hr_offset*60*60,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
    chn_sim_voc_record_1hr = sims_chn_voc_records[[fitOrder[fitRank]]][(hr_offset*60*60+1):(hr_offset*60*60+3600)]
    stripchart(which(chn_sim_voc_record_1hr==1),xaxt="n",main="1 hour within the day",pch=19,ylim=c(.5,1.5),xlim=c(0,3600))
    
    fivemin_offset = 4.5
    rect(xleft=fivemin_offset*60*5,xright=300+fivemin_offset*60*5,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
    chn_sim_voc_record_5min = chn_sim_voc_record_1hr[(fivemin_offset*60*5+1):(fivemin_offset*60*5+300)]
    stripchart(which(chn_sim_voc_record_5min==1),xaxt="n",main="5 minutes within the hour",pch=19,ylim=c(.5,1.5),xlim=c(0,300))
    mtext("Onsets of simulated child vocalizations",side=3, line = 1, outer=TRUE, cex=2)
    
    par(mfrow=c(3,1),cex=1,mar=c(1,1,2,1),oma=c(0,0,4,0))
    stripchart(which(chn_voc_record==1),xaxt="n",main="Full day recording (10 hours)",pch=19,ylim=c(.5,1.5),xlim=c(0,sim_length))
    
    rect(xleft=10000,xright=13600,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
    chn_voc_record_1hr = chn_voc_record[10000:13600]
    stripchart(which(chn_voc_record_1hr==1),xaxt="n",main="1 hour within the day",pch=19,ylim=c(.5,1.5),xlim=c(0,3600))
    
    fivemin_offset = 0
    rect(xleft=700,xright=1000,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
    chn_voc_record_5min = chn_voc_record_1hr[700:1000]
    stripchart(which(chn_voc_record_5min==1),xaxt="n",main="5 minutes within the hour",pch=19,ylim=c(.5,1.5),xlim=c(0,300))
    mtext("Onsets of human child vocalizations",side=3, line = 1, outer=TRUE, cex=2)

#################################################################################
# for the best-matched simulation,
# analyze the chn and adu ivis to test for response effects on ivis
# using 1. no control for previous ivi and 2. control for previous ivi or 2 ivis
#################################################################################

rthresh = 5

simID = 0
ivi_records = c()
ivi_r_records = c()
simIDs=c()
previvi_resids = c()
prev3ivi_resids = c()

for (fitRank in 1:20){
  
  simID = simID + 1
  
  ivi_series = sims_chn_ivi_records[[fitOrder[fitRank]]]
  chn_voc_record = sims_chn_voc_records[[fitOrder[fitRank]]]
  adu_voc_record = sims_adu_voc_records[[fitOrder[fitRank]]]
  
  n_ivi = length(ivi_series)
  
  
  ivi_r_record = c() # initialize the record of when adu responded to chn
  t = 1
  for (i in 1:n_ivi){
    if (ivi_series[i]==1){
      ivi_r_record[i] = NA
    } else if (sum(adu_voc_record[(t+1):(t+rthresh)])>0){
      ivi_r_record[i] = 1
    } else{
      ivi_r_record[i] = 0
    }
    t = t+ivi_series[i]
  }
  
  ivi_records = c(ivi_records,ivi_series[4:n_ivi])
  ivi_r_records = c(ivi_r_records,ivi_r_record[4:n_ivi])
  simIDs = c(simIDs,rep(simID,(n_ivi-3)))
  
  # Correlate current IVI with previous IVI and get the residuals of that correlation
  previvi_model = lm(scale(log(ivi_series[2:n_ivi]))~scale(log(ivi_series[1:(n_ivi-1)])))
  previvi_resid = resid(previvi_model)
  previvi_resids = c(previvi_resids,previvi_resid[3:(n_ivi-1)])
  
  # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
  prev3ivi_model = lm(scale(log(ivi_series[4:n_ivi]))~scale(log(ivi_series[3:(n_ivi-1)]))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])))+(scale(log(ivi_series[3:(n_ivi-1)])+log(ivi_series[2:(n_ivi-2)])+log(ivi_series[1:(n_ivi-3)]))))
  prev3ivi_resid = resid(prev3ivi_model)
  prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid)
  
}

# Analyze using IVI-based approaches (our original and its variants controlling for previous IVIs)
ivi_models = analyze_ivis(ivi_records,ivi_r_records,previvi_resids,prev3ivi_resids,simIDs)
uncontrolled_response_model = ivi_models[[1]]
residual_response_model = ivi_models[[2]]
prev3residual_response_model = ivi_models[[3]]

summary(uncontrolled_response_model)
summary(residual_response_model)
summary(prev3residual_response_model)

#################################################################################
# Plot fit to human data as a function of parameter value
#################################################################################

library(lattice)
library(latticeExtra) 

# showing data points on the same color scale 
nona_sims_df = subset(sims_df, !is.na(simDist) & !is.infinite(simDist))
levelplot(simDist ~ chn_sim_minp * chn_sim_maxp, nona_sims_df, 
          panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
levelplot(simDist ~ adu_sim_minp * adu_sim_maxp, nona_sims_df, 
          panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
levelplot(simDist ~ chn_sim_sdlog * chn_sim_maxp, nona_sims_df, 
          panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
levelplot(simDist ~ adu_sim_sdlog * adu_sim_maxp, nona_sims_df, 
          panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))