library(lme4)
library(lmerTest)

analyze_ivis <- function(ivi_records,ivi_response_records,previvi_resids,prev3ivi_resids,simIDs){
  
  # Compare IVI with vs. without response, without controlling for previous IVI
  uncontrolled_response_model = lmer(scale(log(ivi_records))~ivi_response_records+(ivi_response_records|simIDs))
  if (isSingular(uncontrolled_response_model)){
    uncontrolled_response_model = lmer(scale(log(ivi_records))~ivi_response_records+(1|simIDs)) 
    if (isSingular(uncontrolled_response_model)){
      uncontrolled_response_model = lm(scale(log(ivi_records))~ivi_response_records)
    }
  }
  
  # Compare with vs. without response residuals of the correlation between current IVI and previous IVI
  residual_response_model = lmer(scale(previvi_resids)~ivi_response_records+(ivi_response_records|simIDs))
  if (isSingular(residual_response_model)){
    residual_response_model = lmer(scale(previvi_resids)~ivi_response_records+(1|simIDs)) 
    if (isSingular(residual_response_model)){
      residual_response_model = lm(scale(previvi_resids)~ivi_response_records)
    }
  }
  
  # Compare with vs. without response residuals of the correlation between current IVI and previous 3 IVIs
  prev3residual_response_model = lmer(scale(prev3ivi_resids)~ivi_response_records+(ivi_response_records|simIDs))
  if (isSingular(prev3residual_response_model)){
    prev3residual_response_model = lmer(scale(prev3ivi_resids)~ivi_response_records+(1|simIDs))
    if (isSingular(prev3residual_response_model)){
      prev3residual_response_model = lm(scale(prev3ivi_resids)~ivi_response_records)
    }
  }
  
  return(list(uncontrolled_response_model,residual_response_model,prev3residual_response_model))
  
}