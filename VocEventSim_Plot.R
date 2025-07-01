library(ggplot2)
library(tidyverse)
library(patchwork)

setwd('~/Documents/GitHub/vocal-response-analysis-simulation/')

fitRank = 1
load('data/0344_000913/nonInteractive/VocEventSim_Optimize.RData')
chi_voc = chn_voc_record
adu_voc = adu_voc_record
datasets <- list(data.frame(chi_voc,adu_voc,t = seq(1,length(chi_voc)),type = rep("human",length(chi_voc))))
chi_voc = sims_chn_voc_records[[fitOrder[fitRank]]]
adu_voc = sims_adu_voc_records[[fitOrder[fitRank]]]
datasets = append(datasets,list(data.frame(chi_voc,adu_voc,t = seq(1,length(chi_voc)),type = rep("nonInteractive",length(chi_voc)))))
load('data/0344_000913/a2Interactive/VocEventSim_Optimize.RData')
chi_voc = sims_chn_voc_records[[fitOrder[fitRank]]]
adu_voc = sims_adu_voc_records[[fitOrder[fitRank]]]
datasets = append(datasets,list(data.frame(chi_voc,adu_voc,t = seq(1,length(chi_voc)),type = rep("a2Interactive",length(chi_voc)))))
load('data/0344_000913/bidirectional/VocEventSim_Optimize.RData')
chi_voc = sims_chn_voc_records[[fitOrder[fitRank]]]
adu_voc = sims_adu_voc_records[[fitOrder[fitRank]]]
datasets = append(datasets,list(data.frame(chi_voc,adu_voc,t = seq(1,length(chi_voc)),type = rep("bidirectional",length(chi_voc)))))

plot_quantile = 0.9

################################################################################
# Make raster plots of chi and adu vocs at full day, 1 hour, and 5 minute levels
# for the human and best-matched simulations of each type
################################################################################

make_quadrant_plot <- function(data,quadrant_title = "") {
  
  window_size_1hr = 60*60
  n_1hr <- nrow(data)
  window_sums_1hr <- sapply(1:(n_1hr - window_size_1hr + 1), function(i) {
    sum(data$chi_voc[i:(i + window_size_1hr - 1)])
  })
  target_count_1hr <- quantile(window_sums_1hr,probs=plot_quantile)
  matching_indices_1hr <- seq(which(window_sums_1hr == target_count_1hr)[1],which(window_sums_1hr == target_count_1hr)[1]+window_size_1hr-1)
  df_1hr = data[matching_indices_1hr,]
  
  window_size_5min = 5*60
  chi_voc_1hr = df_1hr$chi_voc
  n_5min <- length(chi_voc_1hr)
  window_sums_5min <- sapply(1:(n_5min - window_size_5min + 1), function(i) {
    sum(chi_voc_1hr[i:(i + window_size_5min - 1)])
  })
  target_count_5min <- quantile(window_sums_5min,probs=plot_quantile)
  matching_indices_5min <- seq(which(window_sums_5min == target_count_5min)[1],which(window_sums_5min == target_count_5min)[1]+window_size_5min-1)
  df_5min = df_1hr[matching_indices_5min,]
  
  p1 <- data %>%
    mutate(f = case_when(
      adu_voc == 1 ~ "red",
      chi_voc == 1 ~ "blue",
      t %in% df_1hr$t ~ "lightgray",
      TRUE ~ "white"
    )) %>%
    ggplot(aes(t,1,fill=f)) +
    geom_tile() +
    scale_fill_identity() +
    theme_classic() +
    ylab("full day") +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  
  p2 <- df_1hr %>%
    mutate(f = case_when(
      adu_voc == 1 ~ "red",
      chi_voc == 1 ~ "blue",
      t %in% df_5min$t ~ "lightgray",
      TRUE ~ "white"
    )) %>%
    ggplot(aes(t,1,fill=f)) +
    geom_tile() +
    scale_fill_identity() +
    theme_classic() +
    ylab("1 hour") +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  
  p3 <- df_5min %>%
    mutate(f = case_when(
      adu_voc == 1 ~ "red",
      chi_voc == 1 ~ "blue",
      TRUE ~ "white"
    )) %>%
    ggplot(aes(t,1,fill=f)) +
    geom_tile() +
    scale_fill_identity() +
    theme_classic() +
    ylab("5 minutes") +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )
  
  wrapped_quad <- wrap_elements(full = (p1 / p2 / p3) + plot_annotation(title = quadrant_title)) & theme(plot.title = element_text(hjust = 0.5,face="bold"))
  
}

quadrants <- mapply(
  make_quadrant_plot,
  datasets,
  quadrant_title = c("Human","No Interaction","Unidirectional Interaction","Bidirectional Interaction"),
  SIMPLIFY = FALSE
)

final_grid <- (quadrants[[1]] | quadrants[[2]]) /
                 (quadrants[[3]] | quadrants[[4]])

ggsave("grid_plot.pdf",plot=final_grid,width=8.5,height=11,units="in")

fitRank = 1
load('data/0344_000913/nonInteractive/VocEventSim_Optimize.RData')
chi_ivi = chn_ivi_record
adu_ivi = adu_ivi_record
chi_ivi_datasets <- list(chi_ivi)
adu_ivi_datasets <- list(adu_ivi)
chi_ivi = sims_chn_ivi_records[[fitOrder[fitRank]]]
adu_voc = sims_adu_voc_records[[fitOrder[fitRank]]]
chi_ivi_datasets = append(chi_ivi_datasets,list(chi_ivi))
adu_ivi_datasets = append(adu_ivi_datasets,list(adu_ivi))
load('data/0344_000913/a2Interactive/VocEventSim_Optimize.RData')
chi_ivi = sims_chn_ivi_records[[fitOrder[fitRank]]]
adu_voc = sims_adu_voc_records[[fitOrder[fitRank]]]
chi_ivi_datasets = append(chi_ivi_datasets,list(chi_ivi))
adu_ivi_datasets = append(adu_ivi_datasets,list(adu_ivi))
load('data/0344_000913/bidirectional/VocEventSim_Optimize.RData')
chi_ivi = sims_chn_ivi_records[[fitOrder[fitRank]]]
adu_voc = sims_adu_voc_records[[fitOrder[fitRank]]]
chi_ivi_datasets = append(chi_ivi_datasets,list(chi_ivi))
adu_ivi_datasets = append(adu_ivi_datasets,list(adu_ivi))
dataTypes = c("Human","No Interaction","Unidirectional Interaction","Bidirectional Interaction")

################################################################################
# Plot log-log *chi* ivi distributions and obtain power law and lognormal fits
# for the human and best-matched simulations of each type
################################################################################

library(poweRlaw)

pdf("chi_ivi.pdf", width = 8.5, height = 11)
par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))

for (n in 1:4){
  chi_ivi = chi_ivi_datasets[[n]]
  adu_ivi = adu_ivi_datasets[[n]]
  dataType = dataTypes[n]
  
  # power law model of human infant ivis
  pl_m = displ$new(chi_ivi)
  est_pl = estimate_xmin(pl_m)
  pl_m$setXmin(est_pl[[2]])
  pl_m$setPars(est_pl[[3]])
  plot(pl_m, xlab = "IVI (s)", ylab = "Count", main = dataTypes[n])
  lines(pl_m, col="purple",lwd = 2)
  
  # lognormal model of human infant ivis
  ln_m = dislnorm$new(chi_ivi)
  ln_m$setXmin(est_pl[[2]])
  est_ln = estimate_xmin(ln_m)
  ln_m$setPars(est_ln[[3]])
  lines(ln_m, col="orange",lwd = 2)
  
  if(n==1){
    legend(x="bottomleft",legend = c("power law fit","lognormal fit"), col=c("purple","orange"),lty=1,lwd = 2)
  }
  
  comp = compare_distributions(pl_m, ln_m)
  print(comp$test_statistic)
  print(comp$p_two_sided)
}

mtext("Child IVI distributions and fits", side = 3, outer = TRUE, cex = 1.5)
dev.off()

################################################################################
# Plot log-log *adu* ivi distributions and obtain power law and lognormal fits
# for the human and best-matched simulations of each type
################################################################################
