setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')
load('data/0344_000913/nonInteractive/VocEventSim_Optimize.RData')

fitRank = 1

library(ggplot2)
hist(log(chn_ivi_record))
hist(log(sims_chn_ivi_records[[fitOrder[fitRank]]]))
hist(log(adu_ivi_record))
hist(log(sims_adu_ivi_records[[fitOrder[fitRank]]]))

library(poweRlaw)

# power law model of human infant ivis
pl_m = displ$new(chn_ivi_record)
est_pl = estimate_xmin(pl_m)
pl_m$setXmin(est_pl[[2]])
pl_m$setPars(est_pl[[3]])
plot(pl_m)
lines(pl_m, col=2)

# lognormal model of human infant ivis
ln_m = dislnorm$new(chn_ivi_record)
ln_m$setXmin(est_pl[[2]])
est_ln = estimate_xmin(ln_m)
ln_m$setPars(est_ln[[3]])
lines(ln_m, col=3)

comp = compare_distributions(pl_m, ln_m)
comp$test_statistic

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

pdf(file="Fig_MultiscaleClusters_Simulation.pdf")

par(mfrow=c(3,1),cex=1,mar=c(1,1,2,1),oma=c(0,0,4,0))
stripchart(which(sims_chn_voc_records[[fitOrder[fitRank]]]==1),xaxt="n",main="Full day simulation (10 hours)",pch=19,ylim=c(.5,1.5),xlim=c(0,sim_length))

hr_offset = 0.2
rect(xleft=hr_offset*60*60,xright=3600+hr_offset*60*60,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
chn_sim_voc_record_1hr = sims_chn_voc_records[[fitOrder[fitRank]]][(hr_offset*60*60+1):(hr_offset*60*60+3600)]
stripchart(which(chn_sim_voc_record_1hr==1),xaxt="n",main="1 hour within the day",pch=19,ylim=c(.5,1.5),xlim=c(0,3600))

fivemin_offset = 8.5
rect(xleft=fivemin_offset*60*5,xright=300+fivemin_offset*60*5,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
chn_sim_voc_record_5min = chn_sim_voc_record_1hr[(fivemin_offset*60*5+1):(fivemin_offset*60*5+300)]
stripchart(which(chn_sim_voc_record_5min==1),xaxt="n",main="5 minutes within the hour",pch=19,ylim=c(.5,1.5),xlim=c(0,300))
mtext("Onsets of simulated child vocalizations",side=3, line = 1, outer=TRUE, cex=2)

dev.off()

pdf(file="Fig_MultiscaleClusters_Human.pdf")

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

dev.off()

print(sims_df[fitOrder[1],]$chn_sim_minp)
print(sims_df[fitOrder[1],]$chn_sim_maxp)
print(sims_df[fitOrder[1],]$chn_sim_meanlog)
print(sims_df[fitOrder[1],]$chn_sim_sdlog)

load('data/0344_000913/a2interactive/VocEventSim_Optimize.RData')

pdf(file="Fig_MultiscaleClusters_Simulation_a2interactive.pdf")

par(mfrow=c(3,1),cex=1,mar=c(1,1,2,1),oma=c(0,0,4,0))
stripchart(which(sims_chn_voc_records[[fitOrder[fitRank]]]==1),xaxt="n",main="Full day simulation (10 hours)",pch=19,ylim=c(.5,1.5),xlim=c(0,sim_length))

hr_offset = 2.8
rect(xleft=hr_offset*60*60,xright=3600+hr_offset*60*60,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
chn_sim_voc_record_1hr = sims_chn_voc_records[[fitOrder[fitRank]]][(hr_offset*60*60+1):(hr_offset*60*60+3600)]
stripchart(which(chn_sim_voc_record_1hr==1),xaxt="n",main="1 hour within the day",pch=19,ylim=c(.5,1.5),xlim=c(0,3600))

fivemin_offset = 8.0
rect(xleft=fivemin_offset*60*5,xright=300+fivemin_offset*60*5,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
chn_sim_voc_record_5min = chn_sim_voc_record_1hr[(fivemin_offset*60*5+1):(fivemin_offset*60*5+300)]
stripchart(which(chn_sim_voc_record_5min==1),xaxt="n",main="5 minutes within the hour",pch=19,ylim=c(.5,1.5),xlim=c(0,300))
mtext("Onsets of simulated child vocalizations",side=3, line = 1, outer=TRUE, cex=2)

dev.off()

# 
# #################################################################################
# # Plot fit to human data as a function of parameter value
# #################################################################################
# 
# library(lattice)
# library(latticeExtra) 
# 
# # showing data points on the same color scale 
# nona_sims_df = subset(sims_df, !is.na(simDist) & !is.infinite(simDist))
# levelplot(simDist ~ chn_sim_minp * chn_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
# levelplot(simDist ~ adu_sim_minp * adu_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
# levelplot(simDist ~ chn_sim_sdlog * chn_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))
# levelplot(simDist ~ adu_sim_sdlog * adu_sim_maxp, nona_sims_df, 
#           panel = panel.levelplot.points, cex = 1.2) + layer_(panel.2dsmoother(..., n = 200))