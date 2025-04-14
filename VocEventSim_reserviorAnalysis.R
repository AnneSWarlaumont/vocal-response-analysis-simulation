source('VocEventSim.R')

sim_length = 10*60*60
rthresh = 1

a2_minp = .001
a2_respsensitivity = 1
a2_othersensitivity = 1000

a1_minp = .000001
a1_respsensitivity = 1
a1_othersensitivity = 1

a1_meanlog = 0
a2_meanlog = 0
a1_sdlog = .2
a2_sdlog = .2
a1_maxp = .4
a2_maxp = .4
a1_p_voc = runif(1,min=a1_minp,max=a1_maxp)
a2_p_voc = runif(1,min=a2_minp,max=a2_maxp)

voc_records = two_agent_vocal_sim(sim_length,a1_p_voc,a2_p_voc,a1_meanlog,a2_meanlog,a1_sdlog,a2_sdlog,a1_minp,a2_minp,a1_maxp,a2_maxp,a1_othersensitivity,a2_othersensitivity,a1_respsensitivity,a2_respsensitivity,rthresh)
t = seq(from = 1, to = sim_length)
a1_voc_record = voc_records[[1]]
a2_voc_record = voc_records[[2]]
a2toa1_r_record = c(voc_records[[3]],NA)
a1toa2_r_record = c(voc_records[[4]],NA)
voc_dat = data.frame(t,a1_voc_record,a2_voc_record,a2toa1_r_record,a1toa2_r_record)

library(reservoirnet)
library(ggplot2)
library(dplyr)
library(slider)
library(ggplot2)

a1_dist_forecast = 60
traintest_sec = 5*60*60

voc_predict_dat <- voc_dat %>%
  # outcome at 14 seconds later
  mutate(a1outcome = lead(x = a1_voc_record, a1_dist_forecast),
         a1outcome_t = t + a1_dist_forecast) %>%
  # rolling average for a1_voc_record and a2toa1_r_record
  mutate_at(.vars = c("a1_voc_record","a2_voc_record","a2toa1_r_record"),
            .funs = function(x) slider::slide_dbl(.x = x,
                                                  .before = 10*60,
                                                  .f = mean,
                                                  na.rm = TRUE))


voc_predict_dat %>%
  tidyr::pivot_longer(cols = c("a1_voc_record","a2_voc_record","a2toa1_r_record")) %>%
  ggplot2::ggplot(mapping = aes(x = t, y = value)) + 
  geom_line() +
  facet_grid(name ~ ., scales = "free_y") +
  theme_bw() +
  geom_vline(mapping = aes(color = "train-test sets", xintercept = traintest_sec)) +
  labs(color = "") +
  theme(legend.position = "bottom")
