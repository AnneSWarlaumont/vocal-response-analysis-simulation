
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