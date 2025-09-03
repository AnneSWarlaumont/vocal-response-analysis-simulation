setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

recordingsToAnalyze = c("0054_000603","0196_000902","0274_000221","0300_000607","0344_000913","0437_010603","0833_010606")
rthresh = 5

################################################################################
# For the human recordings,
# analyze the chn and adu ivis to test for response effects on ivis using:
# 1. control for 3 previous ivis
# 2. control for previous ivi, and
# 3. no control for previous ivi.
#
# (Note to self: We may eventually want to add controlled response beta as a
# criterion for simulation optimization.)
################################################################################

source("get_ivis.R")
source("get_ivi_responses.R")
source("get_ivis_continuous.R")
source("analyze_ivis.R")

chi_human_response_results = data.frame(recording = character(),
                                        rBeta0 = double(),
                                        rBeta0Lower = double(),
                                        rBeta0Upper = double(),
                                        rP0 = double(),
                                        rBeta1 = double(),
                                        rBeta1Lower = double(),
                                        rBeta1Upper = double(),
                                        rP1 = double(),
                                        rBeta3 = double(),
                                        rBeta3Lower = double(),
                                        rBeta3Upper = double(),
                                        rP3 = double())

adu_human_response_results = data.frame(recording = character(),
                                        rBeta0 = double(),
                                        rBeta0Lower = double(),
                                        rBeta0Upper = double(),
                                        rP0 = double(),
                                        rBeta1 = double(),
                                        rBeta1Lower = double(),
                                        rBeta1Upper = double(),
                                        rP1 = double(),
                                        rBeta3 = double(),
                                        rBeta3Lower = double(),
                                        rBeta3Upper = double(),
                                        rP3 = double())

chi_cont_human_response_results = data.frame(recording = character(),
                                             rBeta0 = double(),
                                             rBeta0Lower = double(),
                                             rBeta0Upper = double(),
                                             rP0 = double(),
                                             rBeta1 = double(),
                                             rBeta1Lower = double(),
                                             rBeta1Upper = double(),
                                             rP1 = double(),
                                             rBeta3 = double(),
                                             rBeta3Lower = double(),
                                             rBeta3Upper = double(),
                                             rP3 = double())


chi_ivi_records = c()
chi_ivi_r_records = c()
chi_recordingIDs = c()
chi_previvi_resids = c()
chi_prev3ivi_resids = c()

adu_ivi_records = c()
adu_ivi_r_records = c()
adu_recordingIDs = c()
adu_previvi_resids = c()
adu_prev3ivi_resids = c()

chi_cont_ivi_records = c()
chi_cont_ivi_r_records = c()
chi_cont_recordingIDs = c()
chi_cont_previvi_resids = c()
chi_cont_prev3ivi_resids = c()

# Currently working on getting for the human data:
# chi_previvi_resids, chi_prev3ivi_resids, chi_recordingIDs
# and the same for adu and chi_cont
for (recording in recordingsToAnalyze){
  
  recording_dir = paste("data/",recording,sep="")
  lena_segments = read.csv(paste(recording_dir,"/",recording,"_segments.csv",sep=""))
  chn_segments = subset(lena_segments,segtype=="CHNSP")
  adu_segments = subset(lena_segments,segtype=="FAN"|segtype=="MAN")
  rec_length = floor(lena_segments$endsec[nrow(lena_segments)])
  chi_voc_record = integer(rec_length)
  chi_voc_record[as.integer(floor(chn_segments$startsec))] = 1
  adu_voc_record = integer(rec_length)
  adu_voc_record[as.integer(floor(adu_segments$startsec))] = 1
  
  # get ivi and response records
  chi_ivi_record = get_ivis(chi_voc_record)
  adu_ivi_record = get_ivis(adu_voc_record)
  chi_ivi_r_record <- get_ivi_responses(chi_ivi_record,
                                        adu_voc_record,
                                        rthresh)
  adu_ivi_r_record <- get_ivi_responses(adu_ivi_record,
                                        chi_voc_record,
                                        rthresh)
  temp <- get_ivis_continuous(chn_segments,adu_segments,rthresh)
  chi_cont_ivi_record <- temp[[1]]
  chi_cont_ivi_r_record <- temp[[2]]
  
  chi_n_ivi <- length(chi_ivi_record)
  adu_n_ivi <- length(adu_ivi_record)
  chi_cont_n_ivi <- length(chi_cont_ivi_record)
  
  # Correlate current IVI with previous IVI and get the residuals of that
  # correlation
  
  chi_previvi_model <- lm(scale(log(chi_ivi_record[2:chi_n_ivi]))
                         ~scale(log(chi_ivi_record[1:(chi_n_ivi-1)])))
  chi_previvi_resid <- resid(chi_previvi_model)[3:(chi_n_ivi-1)]
  adu_previvi_model <- lm(scale(log(adu_ivi_record[2:adu_n_ivi]))
                          ~scale(log(adu_ivi_record[1:(adu_n_ivi-1)])))
  adu_previvi_resid <- resid(adu_previvi_model)[3:(adu_n_ivi-1)]
  chi_cont_previvi_model <- lm(scale(log(chi_cont_ivi_record[2:chi_cont_n_ivi]))
                          ~scale(log(chi_cont_ivi_record[1:(chi_cont_n_ivi-1)])))
  chi_cont_previvi_resid <- resid(chi_cont_previvi_model)[3:(chi_cont_n_ivi-1)]
  
  # Correlate current IVI with time since the past 3 IVIs and get the residuals
  # of that correlation 
  chi_prev3ivi_model<-lm(scale(log(chi_ivi_record[4:chi_n_ivi]))
                             ~scale(log(chi_ivi_record[3:(chi_n_ivi-1)]))
                             +(scale(log(chi_ivi_record[3:(chi_n_ivi-1)])
                                     +log(chi_ivi_record[2:(chi_n_ivi-2)])))
                             +(scale(log(chi_ivi_record[3:(chi_n_ivi-1)])
                                     +log(chi_ivi_record[2:(chi_n_ivi-2)])
                                     +log(chi_ivi_record[1:(chi_n_ivi-3)]))))
  chi_prev3ivi_resid <- resid(chi_prev3ivi_model)
  adu_prev3ivi_model<-lm(scale(log(adu_ivi_record[4:adu_n_ivi]))
                             ~scale(log(adu_ivi_record[3:(adu_n_ivi-1)]))
                             +(scale(log(adu_ivi_record[3:(adu_n_ivi-1)])
                                     +log(adu_ivi_record[2:(adu_n_ivi-2)])))
                             +(scale(log(adu_ivi_record[3:(adu_n_ivi-1)])
                                     +log(adu_ivi_record[2:(adu_n_ivi-2)])
                                     +log(adu_ivi_record[1:(adu_n_ivi-3)]))))
  adu_prev3ivi_resid <- resid(adu_prev3ivi_model)
  chi_cont_prev3ivi_model<-lm(scale(log(chi_cont_ivi_record[4:chi_cont_n_ivi]))
                         ~scale(log(chi_cont_ivi_record[3:(chi_cont_n_ivi-1)]))
                         +(scale(log(chi_cont_ivi_record[3:(chi_cont_n_ivi-1)])
                                 +log(chi_cont_ivi_record[2:(chi_cont_n_ivi-2)])))
                         +(scale(log(chi_cont_ivi_record[3:(chi_cont_n_ivi-1)])
                                 +log(chi_cont_ivi_record[2:(chi_cont_n_ivi-2)])
                                 +log(chi_cont_ivi_record[1:(chi_cont_n_ivi-3)]))))
  chi_cont_prev3ivi_resid <- resid(chi_cont_prev3ivi_model)
  
  chi_ivi_record <- chi_ivi_record[4:chi_n_ivi]
  chi_ivi_r_record <- chi_ivi_r_record[4:chi_n_ivi]
  chi_recordingID <- rep(recording,chi_n_ivi-3)
  chi_ivi_records = c(chi_ivi_records,chi_ivi_record)
  chi_ivi_r_records = c(chi_ivi_r_records,chi_ivi_r_record)
  chi_recordingIDs = c(chi_recordingIDs,chi_recordingID)
  chi_previvi_resids = c(chi_previvi_resids,chi_previvi_resid)
  chi_prev3ivi_resids = c(chi_prev3ivi_resids,chi_prev3ivi_resid)
  
  adu_ivi_record <- adu_ivi_record[4:adu_n_ivi]
  adu_ivi_r_record <- adu_ivi_r_record[4:adu_n_ivi]
  adu_recordingID <- rep(recording,adu_n_ivi-3)
  adu_ivi_records = c(adu_ivi_records,adu_ivi_record)
  adu_ivi_r_records = c(adu_ivi_r_records,adu_ivi_r_record)
  adu_recordingIDs = c(adu_recordingIDs,adu_recordingID)
  adu_previvi_resids = c(adu_previvi_resids,adu_previvi_resid)
  adu_prev3ivi_resids = c(adu_prev3ivi_resids,adu_prev3ivi_resid)
  
  chi_cont_ivi_record <- chi_cont_ivi_record[4:chi_cont_n_ivi]
  chi_cont_ivi_r_record <- chi_cont_ivi_r_record[4:chi_cont_n_ivi]
  chi_cont_recordingID <- rep(recording,chi_cont_n_ivi-3)
  chi_cont_ivi_records = c(chi_cont_ivi_records,chi_cont_ivi_record)
  chi_cont_ivi_r_records = c(chi_cont_ivi_r_records,chi_cont_ivi_r_record)
  chi_cont_recordingIDs = c(chi_cont_recordingIDs,chi_cont_recordingID)
  chi_cont_previvi_resids = c(chi_cont_previvi_resids,chi_cont_previvi_resid)
  chi_cont_prev3ivi_resids = c(chi_cont_prev3ivi_resids,chi_cont_prev3ivi_resid)
  
  # Analyze this recording using IVI-based approaches
  
  chi_ivi_models = analyze_ivis(chi_ivi_record,chi_ivi_r_record,chi_previvi_resid,chi_prev3ivi_resid)
  chi_uncontrolled_response_model = chi_ivi_models[[1]]
  chi_residual_response_model = chi_ivi_models[[2]]
  chi_prev3residual_response_model = chi_ivi_models[[3]]
  adu_ivi_models = analyze_ivis(adu_ivi_record,adu_ivi_r_record,adu_previvi_resid,adu_prev3ivi_resid)
  adu_uncontrolled_response_model = adu_ivi_models[[1]]
  adu_residual_response_model = adu_ivi_models[[2]]
  adu_prev3residual_response_model = adu_ivi_models[[3]]
  chi_cont_ivi_models = analyze_ivis(chi_cont_ivi_record,chi_cont_ivi_r_record,chi_cont_previvi_resid,chi_cont_prev3ivi_resid)
  chi_cont_uncontrolled_response_model = chi_cont_ivi_models[[1]]
  chi_cont_residual_response_model = chi_cont_ivi_models[[2]]
  chi_cont_prev3residual_response_model = chi_cont_ivi_models[[3]]
  
  chi_rSummary0 = summary(chi_uncontrolled_response_model)
  chi_rBeta0 = chi_rSummary0$coefficients["ivi_response_records","Estimate"]
  chi_rP0 = chi_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_uncontrolled_response_model)=="lm"){
    chi_rCI0 = confint(chi_uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_rCI0 = confint.merMod(chi_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_rSummary0 = summary(adu_uncontrolled_response_model)
  adu_rBeta0 = adu_rSummary0$coefficients["ivi_response_records","Estimate"]
  adu_rP0 = adu_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_uncontrolled_response_model)=="lm"){
    adu_rCI0 = confint(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_rCI0 = confint.merMod(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  chi_cont_rSummary0 = summary(chi_cont_uncontrolled_response_model)
  chi_cont_rBeta0 = chi_cont_rSummary0$coefficients["ivi_response_records","Estimate"]
  chi_cont_rP0 = chi_cont_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_cont_uncontrolled_response_model)=="lm"){
    chi_cont_rCI0 = confint(chi_cont_uncontrolled_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_cont_rCI0 = confint.merMod(chi_cont_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
  }
  
  chi_rSummary1 = summary(chi_residual_response_model)
  chi_rBeta1 = chi_rSummary1$coefficients["ivi_response_records","Estimate"]
  chi_rP1 = chi_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_residual_response_model)=="lm"){
    chi_rCI1 = confint(chi_residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_rCI1 = confint.merMod(chi_residual_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_rSummary1 = summary(adu_residual_response_model)
  adu_rBeta1 = adu_rSummary1$coefficients["ivi_response_records","Estimate"]
  adu_rP1 = adu_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_residual_response_model)=="lm"){
    adu_rCI1 = confint(adu_residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_rCI1 = confint.merMod(adu_residual_response_model, "ivi_response_records", level = 0.99) 
  }
  chi_cont_rSummary1 = summary(chi_cont_residual_response_model)
  chi_cont_rBeta1 = chi_cont_rSummary1$coefficients["ivi_response_records","Estimate"]
  chi_cont_rP1 = chi_cont_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_cont_residual_response_model)=="lm"){
    chi_cont_rCI1 = confint(chi_cont_residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_cont_rCI1 = confint.merMod(chi_cont_residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  chi_rSummary3 = summary(chi_prev3residual_response_model)
  chi_rBeta3 = chi_rSummary3$coefficients["ivi_response_records","Estimate"]
  chi_rP3 = chi_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_prev3residual_response_model)=="lm"){
    chi_rCI3 = confint(chi_prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_rCI3 = confint.merMod(chi_prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  adu_rSummary3 = summary(adu_prev3residual_response_model)
  adu_rBeta3 = adu_rSummary3$coefficients["ivi_response_records","Estimate"]
  adu_rP3 = adu_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(adu_prev3residual_response_model)=="lm"){
    adu_rCI3 = confint(adu_prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    adu_rCI3 = confint.merMod(adu_prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  chi_cont_rSummary3 = summary(chi_cont_prev3residual_response_model)
  chi_cont_rBeta3 = chi_cont_rSummary3$coefficients["ivi_response_records","Estimate"]
  chi_cont_rP3 = chi_cont_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
  if (class(chi_cont_prev3residual_response_model)=="lm"){
    chi_cont_rCI3 = confint(chi_cont_prev3residual_response_model, "ivi_response_records", level = 0.99)
  } else{
    chi_cont_rCI3 = confint.merMod(chi_cont_prev3residual_response_model, "ivi_response_records", level = 0.99) 
  }
  
  chi_newrow = data.frame(recording = recording,
                          rBeta0 = chi_rBeta0,
                          rBeta0Lower = chi_rCI0[1,1],
                          rBeta0Upper = chi_rCI0[1,2],
                          rP0 = chi_rP0,
                          rBeta1 = chi_rBeta1,
                          rBeta1Lower = chi_rCI1[1,1],
                          rBeta1Upper = chi_rCI1[1,2],
                          rP1 = chi_rP1,
                          rBeta3 = chi_rBeta3,
                          rBeta3Lower = chi_rCI3[1,1],
                          rBeta3Upper = chi_rCI3[1,2],
                          rP3 = chi_rP3)
  chi_human_response_results = rbind(chi_human_response_results,chi_newrow)
  
  adu_newrow = data.frame(recording = recording,
                          rBeta0 = adu_rBeta0,
                          rBeta0Lower = adu_rCI0[1,1],
                          rBeta0Upper = adu_rCI0[1,2],
                          rP0 = adu_rP0,
                          rBeta1 = adu_rBeta1,
                          rBeta1Lower = adu_rCI1[1,1],
                          rBeta1Upper = adu_rCI1[1,2],
                          rP1 = adu_rP1,
                          rBeta3 = adu_rBeta3,
                          rBeta3Lower = adu_rCI3[1,1],
                          rBeta3Upper = adu_rCI3[1,2],
                          rP3 = adu_rP3)
  adu_human_response_results = rbind(adu_human_response_results,adu_newrow)
  
  chi_cont_newrow = data.frame(recording = recording,
                               rBeta0 = chi_cont_rBeta0,
                               rBeta0Lower = chi_cont_rCI0[1,1],
                               rBeta0Upper = chi_cont_rCI0[1,2],
                               rP0 = chi_cont_rP0,
                               rBeta1 = chi_cont_rBeta1,
                               rBeta1Lower = chi_cont_rCI1[1,1],
                               rBeta1Upper = chi_cont_rCI1[1,2],
                               rP1 = chi_cont_rP1,
                               rBeta3 = chi_cont_rBeta3,
                               rBeta3Lower = chi_cont_rCI3[1,1],
                               rBeta3Upper = chi_cont_rCI3[1,2],
                               rP3 = chi_cont_rP3)
  chi_cont_human_response_results = rbind(chi_cont_human_response_results,chi_cont_newrow)
  
}

# Analyze full data (across recordings) using IVI-based approaches (our original
# and its variants controlling for previous IVIs)
chi_ivi_models = analyze_ivis(chi_ivi_records,chi_ivi_r_records,chi_previvi_resids,chi_prev3ivi_resids,as.factor(chi_recordingIDs))
chi_uncontrolled_response_model = chi_ivi_models[[1]]
chi_residual_response_model = chi_ivi_models[[2]]
chi_prev3residual_response_model = chi_ivi_models[[3]]
adu_ivi_models = analyze_ivis(adu_ivi_records,adu_ivi_r_records,adu_previvi_resids,adu_prev3ivi_resids,as.factor(adu_recordingIDs))
adu_uncontrolled_response_model = adu_ivi_models[[1]]
adu_residual_response_model = adu_ivi_models[[2]]
adu_prev3residual_response_model = adu_ivi_models[[3]]
chi_cont_ivi_models = analyze_ivis(chi_cont_ivi_records,chi_cont_ivi_r_records,chi_cont_previvi_resids,chi_cont_prev3ivi_resids,as.factor(chi_cont_recordingIDs))
chi_cont_uncontrolled_response_model = chi_cont_ivi_models[[1]]
chi_cont_residual_response_model = chi_cont_ivi_models[[2]]
chi_cont_prev3residual_response_model = chi_cont_ivi_models[[3]]

chi_rSummary0 = summary(chi_uncontrolled_response_model)
chi_rBeta0 = chi_rSummary0$coefficients["ivi_response_records","Estimate"]
chi_rP0 = chi_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(chi_uncontrolled_response_model)=="lm"){
  chi_rCI0 = confint(chi_uncontrolled_response_model, "ivi_response_records", level = 0.99)
} else{
  chi_rCI0 = confint.merMod(chi_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
}
adu_rSummary0 = summary(adu_uncontrolled_response_model)
adu_rBeta0 = adu_rSummary0$coefficients["ivi_response_records","Estimate"]
adu_rP0 = adu_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(adu_uncontrolled_response_model)=="lm"){
  adu_rCI0 = confint(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99)
} else{
  adu_rCI0 = confint.merMod(adu_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
}
chi_cont_rSummary0 = summary(chi_cont_uncontrolled_response_model)
chi_cont_rBeta0 = chi_cont_rSummary0$coefficients["ivi_response_records","Estimate"]
chi_cont_rP0 = chi_cont_rSummary0$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(chi_cont_uncontrolled_response_model)=="lm"){
  chi_cont_rCI0 = confint(chi_cont_uncontrolled_response_model, "ivi_response_records", level = 0.99)
} else{
  chi_cont_rCI0 = confint.merMod(chi_cont_uncontrolled_response_model, "ivi_response_records", level = 0.99) 
}

chi_rSummary1 = summary(chi_residual_response_model)
chi_rBeta1 = chi_rSummary1$coefficients["ivi_response_records","Estimate"]
chi_rP1 = chi_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(chi_residual_response_model)=="lm"){
  chi_rCI1 = confint(chi_residual_response_model, "ivi_response_records", level = 0.99)
} else{
  chi_rCI1 = confint.merMod(chi_residual_response_model, "ivi_response_records", level = 0.99) 
}
adu_rSummary1 = summary(adu_residual_response_model)
adu_rBeta1 = adu_rSummary1$coefficients["ivi_response_records","Estimate"]
adu_rP1 = adu_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(adu_residual_response_model)=="lm"){
  adu_rCI1 = confint(adu_residual_response_model, "ivi_response_records", level = 0.99)
} else{
  adu_rCI1 = confint.merMod(adu_residual_response_model, "ivi_response_records", level = 0.99) 
}
chi_cont_rSummary1 = summary(chi_cont_residual_response_model)
chi_cont_rBeta1 = chi_cont_rSummary1$coefficients["ivi_response_records","Estimate"]
chi_cont_rP1 = chi_cont_rSummary1$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(chi_cont_residual_response_model)=="lm"){
  chi_cont_rCI1 = confint(chi_cont_residual_response_model, "ivi_response_records", level = 0.99)
} else{
  chi_cont_rCI1 = confint.merMod(chi_cont_residual_response_model, "ivi_response_records", level = 0.99) 
}

chi_rSummary3 = summary(chi_prev3residual_response_model)
chi_rBeta3 = chi_rSummary3$coefficients["ivi_response_records","Estimate"]
chi_rP3 = chi_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(chi_prev3residual_response_model)=="lm"){
  chi_rCI3 = confint(chi_prev3residual_response_model, "ivi_response_records", level = 0.99)
} else{
  chi_rCI3 = confint.merMod(chi_prev3residual_response_model, "ivi_response_records", level = 0.99) 
}
adu_rSummary3 = summary(adu_prev3residual_response_model)
adu_rBeta3 = adu_rSummary3$coefficients["ivi_response_records","Estimate"]
adu_rP3 = adu_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(adu_prev3residual_response_model)=="lm"){
  adu_rCI3 = confint(adu_prev3residual_response_model, "ivi_response_records", level = 0.99)
} else{
  adu_rCI3 = confint.merMod(adu_prev3residual_response_model, "ivi_response_records", level = 0.99) 
}
chi_cont_rSummary3 = summary(chi_cont_prev3residual_response_model)
chi_cont_rBeta3 = chi_cont_rSummary3$coefficients["ivi_response_records","Estimate"]
chi_cont_rP3 = chi_cont_rSummary3$coefficients["ivi_response_records","Pr(>|t|)"]
if (class(chi_cont_prev3residual_response_model)=="lm"){
  chi_cont_rCI3 = confint(chi_cont_prev3residual_response_model, "ivi_response_records", level = 0.99)
} else{
  chi_cont_rCI3 = confint.merMod(chi_cont_prev3residual_response_model, "ivi_response_records", level = 0.99) 
}

chi_newrow = data.frame(recording = "all",
                        rBeta0 = chi_rBeta0,
                        rBeta0Lower = chi_rCI0[1,1],
                        rBeta0Upper = chi_rCI0[1,2],
                        rP0 = chi_rP0,
                        rBeta1 = chi_rBeta1,
                        rBeta1Lower = chi_rCI1[1,1],
                        rBeta1Upper = chi_rCI1[1,2],
                        rP1 = chi_rP1,
                        rBeta3 = chi_rBeta3,
                        rBeta3Lower = chi_rCI3[1,1],
                        rBeta3Upper = chi_rCI3[1,2],
                        rP3 = chi_rP3)
chi_human_response_results = rbind(chi_human_response_results,chi_newrow)

adu_newrow = data.frame(recording = "all",
                        rBeta0 = adu_rBeta0,
                        rBeta0Lower = adu_rCI0[1,1],
                        rBeta0Upper = adu_rCI0[1,2],
                        rP0 = adu_rP0,
                        rBeta1 = adu_rBeta1,
                        rBeta1Lower = adu_rCI1[1,1],
                        rBeta1Upper = adu_rCI1[1,2],
                        rP1 = adu_rP1,
                        rBeta3 = adu_rBeta3,
                        rBeta3Lower = adu_rCI3[1,1],
                        rBeta3Upper = adu_rCI3[1,2],
                        rP3 = adu_rP3)
adu_human_response_results = rbind(adu_human_response_results,adu_newrow)

chi_cont_newrow = data.frame(recording = "all",
                        rBeta0 = chi_cont_rBeta0,
                        rBeta0Lower = chi_cont_rCI0[1,1],
                        rBeta0Upper = chi_cont_rCI0[1,2],
                        rP0 = chi_cont_rP0,
                        rBeta1 = chi_cont_rBeta1,
                        rBeta1Lower = chi_cont_rCI1[1,1],
                        rBeta1Upper = chi_cont_rCI1[1,2],
                        rP1 = chi_cont_rP1,
                        rBeta3 = chi_cont_rBeta3,
                        rBeta3Lower = chi_cont_rCI3[1,1],
                        rBeta3Upper = chi_cont_rCI3[1,2],
                        rP3 = chi_cont_rP3)
chi_cont_human_response_results = rbind(chi_cont_human_response_results,chi_cont_newrow)

write.csv(chi_human_response_results, file = "data/chi_human_response_results.csv")
write.csv(adu_human_response_results, file = "data/adu_human_response_results.csv")
write.csv(chi_cont_human_response_results, file = "data/chi_cont_human_response_results.csv")