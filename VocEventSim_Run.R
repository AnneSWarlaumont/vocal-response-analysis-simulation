source('VocEventSim.R')

# Run a batch of simulations and analyze data:

simcounter = 0
voc_and_resp_records = list()
ivi_records = list()

sink(file="VocEventSimOutput.txt")
sim_length = 10*60*60

for (rthresh in c(1)){#,5)){
  for (a2_othersensitivity in c(1,1.5)){#,2,3,100)){
    for (a2_respsensitivity in c(1,1.5)){#,2,3,100)){
      for (a1_othersensitivity in c(1,1.5)){#,2,3,100)){
        for (a1_respsensitivity in c(1,1.5)){#,2,3,100)){
  
          print("Simulation parameters:")
          print(paste("response threshold:",rthresh))
          print(paste("agent 1 (i.e. infant) other sensitivity:",
                      a1_othersensitivity))
          print(paste("agent 1 response sensitivity:",
                      a1_respsensitivity))
          print(paste("agent 2 (i.e. adult) other sensitivity:",
                      a2_othersensitivity))
          print(paste("agent 2 response sensitivity:",
                      a2_respsensitivity))
          print("* The values indicate the factor that probability of vocalizing is multiplied by after the other agent vocalizes or responds. 1 is thus completely insensitive and higher values indicate higher degrees of sensitivity.")
          
          simID = 0
          a1_ivi_records = c()
          a1_ivi_response_records = c()
          simIDs=c()
          previvi_resids = c()
          prev3ivi_resids = c()
          
          for (i in 1:200){
            
            #simcounter = simcounter+1
            #print(paste("Simulation number:",simcounter))
            
            simID = simID+1
            
            a1_meanlog = 0
            a2_meanlog = 0
            a1_sdlog = .2
            a2_sdlog = .2
            a1_minp = .000001
            a2_minp = .00001
            a1_maxp = .4
            a2_maxp = .4
            a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
            a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)
            
            voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
            a1_voc_record = voc_records[[1]]
            a2_voc_record = voc_records[[2]]
            a2toa1_r_record = voc_records[[3]]
            a2toa2_r_record = voc_records[[4]]
            
            ivi_records = get_ivis_and_responses(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record) # This is not yet updated to handle a1toa2_r_record
            
            a1_ivi_record = ivi_records[[1]]
            n_ivi = length(a1_ivi_record)
            a1_ivi_records = c(a1_ivi_records,a1_ivi_record[1:(n_ivi-1)])
            
            # Correlate current IVI with previous IVI and get the residuals of that correlation
            previvi_model = lm(scale(log(a1_ivi_record[2:n_ivi]))~scale(log(a1_ivi_record[1:(n_ivi-1)])))
            previvi_resid = resid(previvi_model)
            previvi_resids = c(previvi_resids,previvi_resid[1:(n_ivi-1)])
            
            # Correlate current IVI with time since the past 3 a1 IVIs and get the residuals of that correlation
            prev3ivi_model = lm(scale(log(a1_ivi_record[4:n_ivi]))~scale(log(a1_ivi_record[3:(n_ivi-1)]))+(scale(log(a1_ivi_record[3:(n_ivi-1)])+log(a1_ivi_record[2:(n_ivi-2)])))+(scale(log(a1_ivi_record[3:(n_ivi-1)])+log(a1_ivi_record[2:(n_ivi-2)])+log(a1_ivi_record[1:(n_ivi-3)]))))
            prev3ivi_resid = resid(prev3ivi_model)
            prev3ivi_resids = c(prev3ivi_resids,prev3ivi_resid[1:(n_ivi-1)])
            
            a1_ivi_response_record = ivi_records[[2]]
            a1_ivi_response_records = c(a1_ivi_response_records,a1_ivi_response_record[2:length(a1_ivi_response_record)])
            
            simIDs = c(simIDs,rep(simID,n_ivi-1))
            
            voc_and_resp_records = c(voc_and_resp_records,list(data.frame(a1_voc_record,a2_voc_record,c(a2toa1_r_record,rep(NA,times=rthresh)))))
            ivi_records = c(ivi_records,list(data.frame(a1_ivi_record,a1_ivi_response_record,c(NA,previvi_resid))))
            
          }
          
          ivi_models = analyze_ivis(a1_ivi_records,a1_ivi_response_records,previvi_resids,prev3ivi_resids,simIDs)
          a1_uncontrolled_response_model = ivi_models[[1]]
          a1_residual_response_model = ivi_models[[2]]
          a1_prev3residual_response_model = ivi_models[[3]]
    
          print(summary(a1_uncontrolled_response_model))
          print(summary(a1_residual_response_model))
          print(summary(a1_prev3residual_response_model))
          hist(a1_ivi_record)
          
        }
      }
    }
  }
}

# # Run a single simulation and visualize the data.
# sim_length = 10*60*60
# rthresh = 1
# a2_othersensitivity = 1
# a2_respsensitivity = 1
# a1_othersensitivity = 1
# a1_respsensitivity = 1
# # a1_meanlog = 0
# # a2_meanlog = 0
# # a1_sdlog = .2
# # a2_sdlog = .2
# # a1_minp = .000001
# # a2_minp = .000001
# # a1_maxp = .4
# # a2_maxp = .4
# a1_meanlog = 0
# a2_meanlog = 0
# a1_sdlog = .2
# a2_sdlog = .2
# a1_minp = .000001
# a2_minp = .1
# a1_maxp = .4
# a2_maxp = .95
# a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
# a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)
# 
# voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
# a1_voc_record = voc_records[[1]]
# a2_voc_record = voc_records[[2]]
# a2toa1_r_record = voc_records[[3]]
# a2toa2_r_record = voc_records[[4]]
# 
# ivi_records = get_ivis_and_responses(a1_voc_record,a2_voc_record,rthresh,a2toa1_r_record) # This is not yet updated to handle a1toa2_r_record
# 
# a1_ivi_record = ivi_records[[1]]
# n_ivi = length(a1_ivi_record)
# 
# hist(a1_ivi_record)
# ggplot(data.frame(a1_ivi_record), aes(x = a1_ivi_record)) +
#   geom_histogram(bins = 30) +
#   scale_x_log10() +
#   scale_y_log10()
# 
# barcode_data = data.frame(
#   index = seq_along(a1_voc_record),
#   value = a1_voc_record
# )
# ggplot(barcode_data, aes(x = index, xend = index, y = 0, yend = 1, color = factor(value))) +
#   geom_segment(size = 5) +
#   scale_color_manual(values = c("white","black")) +
#   labs("Agent 1 vocalization times",
#        x = "Time (s)",
#        y = "") +
#   theme_minimal() +
#   theme(legend.position = "none")
# 
# #pdf(file="SingleSim_MultiscaleClusters.pdf")
# 
# par(mfrow=c(3,1),cex=1,mar=c(1,1,2,1),oma=c(0,0,4,0))
# stripchart(which(a1_voc_record==1),xaxt="n",main="Full day simulation (10 hours)",pch=19,ylim=c(.5,1.5),xlim=c(0,sim_length))
# 
# hr_offset = 2
# rect(xleft=3600+hr_offset*60*60,xright=7200+hr_offset*60*60,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
# a1_voc_record_1hr = a1_voc_record[(3601+hr_offset*60*60):(7200+hr_offset*60*60)]
# stripchart(which(a1_voc_record_1hr==1),xaxt="n",main="1 hour within the day",pch=19,ylim=c(.5,1.5),xlim=c(0,3600))
# 
# fivemin_offset = 5
# rect(xleft=0+fivemin_offset*60*5,xright=300+fivemin_offset*60*5,ybottom=0,ytop=2,col=rgb(0.5,0.5,0.5,.3),border=NA)
# a1_voc_record_5min = a1_voc_record_1hr[(1+fivemin_offset*60*5):(300+fivemin_offset*60*5)]
# stripchart(which(a1_voc_record_5min==1),xaxt="n",main="5 minutes within the hour",pch=19,ylim=c(.5,1.5),xlim=c(0,300))
# 
# mtext("Onsets of child vocalizations",side=3, line = 1, outer=TRUE, cex=2)
# 
# #dev.off()
# 
# sum(a1_voc_record)
# sum(a2_voc_record)