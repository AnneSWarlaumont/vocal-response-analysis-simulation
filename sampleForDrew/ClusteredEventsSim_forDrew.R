# By Anne S. Warlaumont For Drew H. Abney, Dec. 16, 2024
# Instructions: "Generate 100 x 2000 s spike trains from your cool algorithm with 400 events each (so 400 1's and 1600 0's)"

# I fiddled around to set sdlog, minp, and maxp so that we're getting about 400 events per series.
# We then take the runs that have at least 400 events and randomly remove events from those series to get to exactly 400.

single_sim <- function(sim_length,p_event,meanlog,sdlog,minp,maxp) {
  
  event_record = c() # initialize a variable that will record when agent 1 (the infant) vocalizes
  
  for (t in 1:sim_length){
    p_event_multiplier = rlnorm(1,meanlog=meanlog,sdlog=sdlog)
    p_event = max(minp,min(maxp,p_event_multiplier*p_event))
    
    # Determine whether an event occurs at time t
    event = rbinom(1,1,p_event)
    event_record = c(event_record,event)
  }
  
  return(event_record)
  
}

n_trains = 100
sim_length = 2000
target_n_events = 400
meanlog = 0
sdlog = .2
minp = .05
maxp = .5
p_event = runif(1,min=minp,max=maxp)

train_counter = 0
while(train_counter < n_trains){
  event_record = single_sim(sim_length,p_event,meanlog,sdlog,minp,maxp)
  n_events = sum(event_record)
  if (n_events >= target_n_events){
    n_over = n_events - target_n_events
    event_indices = which(event_record==1)
    indices_to_remove = sample(event_indices,n_over)
    event_record[indices_to_remove] = 0
    train_counter = train_counter + 1
    n_events = sum(event_record)
    if (train_counter == 1){
      event_trains = data.frame(event_record)
    }else{
      event_trains = cbind(event_trains,event_record) 
    }
  }
}

write.table(event_trains,'event_trains_for_Drew.csv',row.names=FALSE,col.names=FALSE,sep=',')
